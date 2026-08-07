"""
0_pre_prep.py
=============
Reconstruct nominal investment weights by detailed asset type from BEA
Fixed Assets Table 2.7. Output mimics PQdat.txt: annual values applied
to every quarter, plus a 0/1 bonus-eligibility vector.

Outputs (in ./raw/):
    PQdat_BEA.csv     -- quarterly weights, 1958Q1 to current
    treat_BEA.csv     -- 0/1 eligibility vector
    asset_types.txt   -- numbered asset list with eligibility flags

Eligibility (Section 168(k) bonus depreciation):
    MACRS recovery period <= 20 years  -> 1 (eligible)
    Software (3-year)                   -> 1
    R&D and entertainment originals     -> 1 (special inclusions)
    Recovery period > 20 years          -> 0 (not eligible)
"""

#%%
# =============================================================================
# imports and config

import os
import time
import requests
import pandas as pd
import numpy as np

OUT_DIR = r"./raw"
os.makedirs(OUT_DIR, exist_ok=True)

START_YEAR = 1958
TABLE_INVESTMENT_BY_TYPE = "FAAt207"


def _load_keys(path):
    keys = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            keys[k.strip()] = v.strip().strip('"').strip("'")
    return keys

BEA_KEY = _load_keys(os.path.join(OUT_DIR, "api_keys.txt"))["BEA_KEY"]


#%%
# =============================================================================
# BEA Fixed Assets API helper

def bea_fa(table_name, freq="A"):
    url = "https://apps.bea.gov/api/data/"
    params = {"UserID": BEA_KEY, "method": "GetData",
              "DataSetName": "FixedAssets",
              "TableName": table_name, "Frequency": freq, "Year": "ALL",
              "ResultFormat": "JSON"}
    for attempt in range(4):
        try:
            r = requests.get(url, params=params, timeout=180)
            r.raise_for_status()
            results = r.json().get("BEAAPI", {}).get("Results", {})
            if isinstance(results, dict) and "Error" in results:
                raise RuntimeError(f"BEA API error on {table_name}: {results['Error']}")
            data = results.get("Data", [])
            if not data:
                raise RuntimeError(f"BEA returned no data for {table_name}")
            df = pd.DataFrame(data)
            df["DataValue"]  = pd.to_numeric(df["DataValue"].astype(str).str.replace(",", ""), errors="coerce")
            df["LineNumber"] = pd.to_numeric(df["LineNumber"], errors="coerce")
            df["YEAR"]       = pd.to_numeric(df["TimePeriod"], errors="coerce").astype(int)
            if "UNIT_MULT" in df.columns:
                um = pd.to_numeric(df["UNIT_MULT"], errors="coerce").fillna(0).astype(int)
                df["DataValue"] = df["DataValue"] * (10.0 ** (um - 9))
            return df
        except (requests.ConnectionError, requests.Timeout):
            if attempt < 3:
                time.sleep(2.0 ** attempt)
            else:
                raise


#%%
# =============================================================================
# pull Table 2.7

print(f"Downloading {TABLE_INVESTMENT_BY_TYPE}...")
fa_all = bea_fa(TABLE_INVESTMENT_BY_TYPE)
print(f"  rows: {len(fa_all)}")
print(f"  unique asset types: {fa_all['LineDescription'].nunique()}")
print(f"  year coverage: {fa_all['YEAR'].min()}-{fa_all['YEAR'].max()}")


#%%
# =============================================================================
# pivot to wide

line_meta = (fa_all[["LineNumber","LineDescription"]]
             .drop_duplicates().sort_values("LineNumber")
             .reset_index(drop=True))
wide = fa_all.pivot_table(index="LineNumber", columns="YEAR",
                          values="DataValue", aggfunc="first").sort_index()


#%%
# =============================================================================
# identify aggregate (parent) rows -- two passes
#
# Pass 1: numerical sum check. A row is aggregate if its value approximately
# equals the cumulative sum of subsequent rows up to some point.
# Pass 2: name-based blacklist for known parent categories that the
# numerical check might miss.

def is_aggregate(line_no, year, tol=0.01):
    if line_no not in wide.index or pd.isna(wide.loc[line_no, year]):
        return False
    val = wide.loc[line_no, year]
    if val == 0:
        return False
    later = wide.index[wide.index > line_no].tolist()
    cum = 0
    for k, ln in enumerate(later):
        v = wide.loc[ln, year]
        if pd.isna(v):
            continue
        cum += v
        if k >= 1 and abs(cum - val) / abs(val) < tol:
            return True
        if cum > val * (1 + tol) and k >= 1:
            break
    return False

# Try multiple recent years -- a row is aggregate if it looks aggregate in
# any of them (some aggregate rows happen to equal sub-sums in only some
# years due to rounding, but consistently across multiple years if true)
ref_years = [y for y in [wide.columns.max(), wide.columns.max()-1,
                          wide.columns.max()-5] if y in wide.columns]
aggregates = set()
for ln in wide.index:
    if any(is_aggregate(ln, yr) for yr in ref_years):
        aggregates.add(ln)

# Backstop: drop known parent categories by description
TOP_AGGREGATE_NAMES = {
    "private fixed assets",
    "private fixed investment",
    "equipment",
    "structures",
    "intellectual property products",
    "nonresidential",
    "residential",
    "nonresidential equipment",
    "residential equipment",
    "nonresidential structures",
    "residential structures",
    "nonresidential intellectual property products",
    "residential intellectual property products",
    "transportation equipment",
    "software",                       # parent for prepackaged/custom/own-account
    "research and development",       # parent for industry breakdown
    "entertainment, literary, and artistic originals",  # parent for movies/books/etc
    "other structures",
}
for ln, desc in line_meta[["LineNumber","LineDescription"]].itertuples(index=False):
    if desc.strip().lower() in TOP_AGGREGATE_NAMES:
        aggregates.add(ln)

leaf_lines = sorted(set(wide.index) - aggregates)
print(f"\n{len(leaf_lines)} leaf rows kept; {len(aggregates)} aggregates dropped.")


#%%
# =============================================================================
# build wide DataFrame with unique column names

desc_lookup = dict(zip(line_meta["LineNumber"], line_meta["LineDescription"]))
combined = wide.loc[leaf_lines].T

desc_counts = pd.Series([desc_lookup[ln] for ln in combined.columns]).value_counts()
new_cols = []
for ln in combined.columns:
    desc = desc_lookup[ln]
    if desc_counts[desc] > 1:
        new_cols.append(f"{desc} (L{int(ln)})")
    else:
        new_cols.append(desc)
combined.columns = new_cols
assert combined.columns.is_unique


#%%
# =============================================================================
# eligibility flags
#
# Equipment and IP -> eligible (5- or 7-year, software 3-year, R&D 5-year)
# Quasi-structures (utilities, telecom, oil/gas, farm, railroad) -> eligible
# Commercial / industrial / residential structures -> not eligible

ELIGIBILITY_RULES = [
    # === Equipment (5- or 7-year MACRS) -- eligible ===
    ("computer",                            1),
    ("communication equipment",             1),
    ("medical equipment",                   1),
    ("nonmedical instruments",              1),
    ("photocopy",                           1),
    ("office and accounting",               1),
    ("furniture and fixtures",              1),
    ("fabricated metal",                    1),
    ("engines",                             1),
    ("metalworking",                        1),
    ("special industry",                    1),
    ("general industrial",                  1),
    ("electrical transmission",             1),
    ("electrical equipment",                1),
    ("light trucks",                        1),
    ("trucks, buses",                       1),
    ("autos",                               1),
    ("aircraft",                            1),
    ("ships and boats",                     1),
    ("railroad equipment",                  1),
    ("agricultural machinery",              1),
    ("construction machinery",              1),
    ("mining and oilfield",                 1),
    ("service industry machinery",          1),
    ("other equipment",                     1),
    ("other nonresidential equipment",      1),
    ("other computer",                      1),

    # === IP products -- eligible ===
    # Software sub-categories (3-year MACRS):
    ("prepackaged",                         1),
    ("custom",                              1),
    ("own account",                         1),
    # R&D by industry (5-year amortization, eligible per HS):
    ("pharmaceutical",                      1),
    ("chemical manufacturing",              1),
    ("semiconductor",                       1),
    ("computer and electronic product",     1),
    ("motor vehicles",                      1),
    ("aerospace",                           1),
    ("other manufacturing",                 1),
    ("scientific research",                 1),
    ("all other nonmanufacturing",          1),
    ("universities",                        1),
    ("nonprofit institutions",              1),
    # Entertainment originals (eligible per HS treatment):
    ("theatrical movies",                   1),
    ("long-lived television",               1),
    ("books",                               1),
    ("music",                               1),

    # === Quasi-structures (15-20 year MACRS) -- eligible ===
    ("electric light and power",            1),
    ("petroleum and natural gas",           1),
    ("telecommunications",                  1),
    ("mining exploration",                  1),
    ("electric",                            1),  # utility electric structures
    ("other power",                         1),  # utility power structures
    ("communication",                       1),  # communication structures (15yr)
    ("farm",                                1),  # farm structures (20yr)
    ("railroads",                           1),  # railroad structures (15yr)

    # === Non-residential structures (39-year MACRS) -- not eligible ===
    ("commercial",                          0),
    ("hospital",                            0),
    ("medical buildings",                   0),
    ("manufacturing",                       0),
    ("industrial",                          0),
    ("religious",                           0),
    ("educational",                         0),
    ("amusement",                           0),
    ("lodging",                             0),
    ("warehouses",                          0),
    ("food and beverage",                   0),
    ("multimerchandise",                    0),
    ("office",                              0),  # office structures specifically
    ("health care",                         0),  # health care buildings
    ("special care",                        0),
    ("mining",                              0),  # mining structures (when not "machinery")
    ("air",                                 0),  # air transport structures
    ("land",                                0),

    # === Residential (27.5-year MACRS) -- not eligible ===
    ("housing units",                       0),
    ("manufactured homes",                  0),
    ("permanent site",                      0),
    ("1 to 4 unit",                         0),
    ("5-or more-unit",                      0),
    ("other residential",                   0),

    # === Misc -- not eligible ===
    ("brokers' commissions",                0),
    ("improvements",                        0),
    ("business",                            0),  # generic "Business" R&D parent if not caught
    ("other",                               0),  # catchall for "Other (L...)"
]

def eligibility_for(col_name):
    cn = col_name.lower()
    for needle, val in ELIGIBILITY_RULES:
        if needle in cn:
            return val
    return 0

treat = pd.Series({c: eligibility_for(c) for c in combined.columns})
print(f"\nEligibility: {int(treat.sum())} eligible / {len(treat)} total.")
unmatched = [c for c in combined.columns
             if not any(needle in c.lower()
                        for needle, _ in ELIGIBILITY_RULES)]
if unmatched:
    print(f"Unmatched (defaulted to 0):")
    for c in unmatched:
        print(f"  {c}")


#%%
# =============================================================================
# annual -> quarterly

end_year = combined.dropna(how="all").index.max()
qper = pd.period_range(f"{START_YEAR}Q1", f"{end_year}Q4", freq="Q")

def annual_to_quarterly(annual_series, qper):
    out = pd.Series(np.nan, index=qper)
    for q in qper:
        y = q.year
        if y in annual_series.index and not pd.isna(annual_series.loc[y]):
            out.loc[q] = annual_series.loc[y]
    return out

pq_q = pd.DataFrame({c: annual_to_quarterly(combined[c], qper)
                     for c in combined.columns}, index=qper)
pq_q = pq_q.dropna(axis=1, how="all")
treat = treat.loc[pq_q.columns]


#%%
# =============================================================================
# write outputs

pq_path    = os.path.join(OUT_DIR, "PQdat_BEA.csv")
treat_path = os.path.join(OUT_DIR, "treat_BEA.csv")
types_path = os.path.join(OUT_DIR, "asset_types.txt")

pq_q.to_csv(pq_path, index_label="QPER")
treat.to_frame("treat").T.to_csv(treat_path, index=False)
with open(types_path, "w", encoding="utf-8") as f:
    for i, c in enumerate(pq_q.columns, 1):
        f.write(f"{i:>3d}  {int(treat.loc[c])}  {c}\n")

print(f"\nWrote:")
print(f"  {pq_path}     ({len(pq_q)} rows, {pq_q.shape[1]} types)")
print(f"  {treat_path}  ({pq_q.shape[1]} eligibility flags)")
print(f"  {types_path}  (numbered asset list with treat flags)")


#%%
# =============================================================================
# sanity check: aggregate e_tau using BEA-derived weights

bonus = pd.Series(0.0, index=qper)
bonus.loc[pd.period_range("2001Q4", "2003Q2", freq="Q")] = 0.30
bonus.loc[pd.period_range("2003Q3", "2004Q4", freq="Q")] = 0.50
bonus.loc[pd.period_range("2008Q1", "2010Q4", freq="Q")] = 0.50
bonus.loc[pd.period_range("2011Q1", "2011Q4", freq="Q")] = 1.00
bonus.loc[pd.period_range("2012Q1", "2013Q4", freq="Q")] = 0.50
bonus.loc[pd.period_range("2014Q1", "2017Q3", freq="Q")] = 0.50
bonus.loc[pd.period_range("2017Q4", "2022Q4", freq="Q")] = 1.00
bonus.loc[pd.period_range("2023Q1", "2023Q4", freq="Q")] = 0.80
bonus.loc[pd.period_range("2024Q1", "2024Q4", freq="Q")] = 0.60
last_q = pd.Period(f"{end_year}Q4", "Q")
if pd.Period("2025Q1", "Q") <= last_q:
    bonus.loc[pd.Period("2025Q1", "Q")] = 0.40
if pd.Period("2025Q2", "Q") <= last_q:
    bonus.loc[pd.period_range("2025Q2", last_q, freq="Q")] = 1.00

shares = pq_q.div(pq_q.sum(axis=1), axis=0)
eligible_share = shares.mul(treat, axis=1).sum(axis=1, min_count=1)
e_tau = bonus * eligible_share

print(f"\ne_tau sanity check:")
print(f"  mean over 2001Q4+ when bonus > 0:    {e_tau.loc[bonus > 0].mean():.4f}")
print(f"  eligible_share at 2007Q1:            {eligible_share.loc[pd.Period('2007Q1','Q')]:.4f}")
print(f"  eligible_share, latest quarter:      {eligible_share.iloc[-1]:.4f}")
print(f"  eligible_share, 1965Q1:              {eligible_share.loc[pd.Period('1965Q1','Q')]:.4f}")
print(f"  eligible_share, 1985Q1:              {eligible_share.loc[pd.Period('1985Q1','Q')]:.4f}")
