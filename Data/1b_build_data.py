"""
build_data.py
=============
Construct all observables for the extended Zubairy (2014) DSGE estimation:

  Zubairy baseline (8): chat, ihat, ghat, bhat, pi, r, tau_l_hat, tau_k_hat
  Extensions      (4): lhat, trhat, tau_c_hat, tau_d_hat
  Tax-policy: tau_itc_hat from House-Shapiro files; e_tau_hat from BEA Fixed Assets weights and statutory bonus calendar.
  Additional tax rates: tau_i_hat and tau_div_hat from McGrattan-TAXSIM, 1960–2013.

All raw data is pulled from FRED API and BEA NIPA API. The HS replication
files PQdat_BEA.csv and treat_BEA.csv (produced by 0_pre_prep.py) are
read locally for tau_itc and e_tau weights.

Conventions (matching US_Data_forMatlab_old.xlsx):
- Real per-capita (chat, ihat, ghat, bhat, lhat, trhat): 100 * log-linear-detrend
- r        = FEDFUNDS / 400         (quarterly rate, decimal)
- pi       = GDPDEF[t] / GDPDEF[t-1]   (gross inflation, ratio)
- tau_*_hat = (tau - mean) / mean    (fractional deviation); tau_l also detrended
- tau_ic_hat, e_tau_hat            = level - mean (rates near zero, no detrend)
"""

#%%
# =============================================================================
# imports and configuration

import os
import time
import requests
import numpy as np
import pandas as pd

def _load_api_keys(path):
    """Read 'KEY=value' pairs from a simple text file. Returns a dict."""
    keys = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            keys[k.strip()] = v.strip().strip('"').strip("'")
    return keys

_keys = _load_api_keys("./raw/api_keys.txt")
FRED_KEY = _keys["FRED_KEY"]
BEA_KEY  = _keys["BEA_KEY"]

# Sample window
START_Q = "1958Q1"
END_Q   = "2026Q1"

# Folders (edit to your environment)
RAW_DIR = r"./raw"        # contains PQdat_BEA.csv, treat_BEA.csv (from 0_pre_prep.py)
OUT_DIR = r"./output"
os.makedirs(OUT_DIR, exist_ok=True)

# Population: "FRED" -> CNP16OV (monthly, average to quarterly)
#             "CSV"  -> read local LNU00000000Q.csv (Zubairy's exact series)
POP_SOURCE   = "FRED"
POP_CSV_PATH = r"./raw/LNU00000000Q.csv"


#%%
# =============================================================================
# helpers: FRED API, BEA NIPA API, transforms

def _request_with_retry(method, url, *, params=None, timeout=120, max_retries=4, backoff=2.0):
    """
    Wrap requests.get/post with exponential-backoff retry. Retries on:
    - Connection timeouts and read timeouts
    - ConnectionError (DNS, refused, reset)
    - 5xx server errors and 429 rate-limit
    Does NOT retry on 4xx client errors (except 429), since those indicate
    the request itself is wrong and won't succeed on retry.
    """
    last_exc = None
    for attempt in range(max_retries):
        try:
            r = requests.request(method, url, params=params, timeout=timeout)
            if r.status_code < 400:
                return r
            if r.status_code == 429 or r.status_code >= 500:
                last_exc = requests.HTTPError(f"HTTP {r.status_code}")
            else:
                r.raise_for_status()  # 4xx -> raise immediately, no retry
        except (requests.ConnectionError, requests.Timeout) as e:
            last_exc = e
        if attempt < max_retries - 1:
            wait = backoff ** attempt
            print(f"    [retry {attempt+1}/{max_retries-1}] {type(last_exc).__name__}; waiting {wait:.1f}s")
            time.sleep(wait)
    raise last_exc


def fred_get(series_id, start="1947-01-01", end=None):
    """Pull one FRED series by ID. Returns a pd.Series with DatetimeIndex."""
    url = "https://api.stlouisfed.org/fred/series/observations"
    params = {"series_id": series_id, "api_key": FRED_KEY, "file_type": "json",
              "observation_start": start}
    if end:
        params["observation_end"] = end
    r = _request_with_retry("GET", url, params=params, timeout=120, max_retries=4)
    obs = r.json().get("observations", [])
    if not obs:
        raise RuntimeError(f"FRED returned no data for {series_id}")
    df = pd.DataFrame(obs)
    df["date"]  = pd.to_datetime(df["date"])
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    s = df.set_index("date")["value"].dropna()
    s.name = series_id
    return s


def bea_nipa(table_name, freq="Q"):
    """Pull a full NIPA table from BEA. Returns long-format DataFrame."""
    url = "https://apps.bea.gov/api/data/"
    params = {"UserID": BEA_KEY, "method": "GetData", "DataSetName": "NIPA",
              "TableName": table_name, "Frequency": freq, "Year": "ALL",
              "ResultFormat": "JSON"}
    r = _request_with_retry("GET", url, params=params, timeout=180, max_retries=4)
    js = r.json()
    results = js.get("BEAAPI", {}).get("Results", {})
    if isinstance(results, dict) and "Error" in results:
        raise RuntimeError(f"BEA API error on {table_name}: {results['Error']}")
    data = results.get("Data", [])
    if not data:
        raise RuntimeError(f"BEA returned no data for {table_name}")
    df = pd.DataFrame(data)
    df["DataValue"]  = pd.to_numeric(df["DataValue"].astype(str).str.replace(",", ""), errors="coerce")
    df["LineNumber"] = pd.to_numeric(df["LineNumber"], errors="coerce")
    df["QPER"]       = pd.PeriodIndex(df["TimePeriod"], freq="Q")

    # BEA returns DataValue in raw units; UNIT_MULT tells us the power of 10.
    # E.g., UNIT_MULT=6 -> values are in millions; we want billions to match FRED.
    if "UNIT_MULT" in df.columns:
        unit_mult = pd.to_numeric(df["UNIT_MULT"], errors="coerce").fillna(0).astype(int)
        # Convert to billions of dollars: multiply by 10^(UNIT_MULT - 9)
        df["DataValue"] = df["DataValue"] * (10.0 ** (unit_mult - 9))
        mults = unit_mult.unique()
        print(f"    [{table_name}] UNIT_MULT values seen: {sorted(mults)} -> rescaled to $bn")
    else:
        print(f"    [{table_name}] WARNING: no UNIT_MULT field; assuming values are already in $bn")
    return df


def nipa_line(df_table, line_description, line_number=None):
    """Extract one line from a NIPA table by description (substring match)."""
    sub = df_table[df_table["LineDescription"].str.contains(
        line_description, case=False, regex=False, na=False)]
    if line_number is not None:
        sub = sub[sub["LineNumber"] == line_number]
    line_nums = sub["LineNumber"].unique()
    if len(line_nums) == 0:
        raise ValueError(f"No NIPA line matched '{line_description}' (line {line_number})")
    if len(line_nums) > 1:
        opts = sub[["LineNumber","LineDescription"]].drop_duplicates().sort_values("LineNumber")
        raise ValueError(
            f"Multiple lines matched '{line_description}'. "
            f"Pass line_number= to disambiguate. Options:\n{opts.to_string(index=False)}")
    s = sub.set_index("QPER")["DataValue"].sort_index()
    s.name = f"{line_description} (L{int(line_nums[0])})"
    return s


def to_qper(s_dt):
    out = s_dt.copy()
    out.index = out.index.to_period("Q")
    return out

def monthly_to_quarterly(s_m, how="mean"):
    g = s_m.groupby(s_m.index.to_period("Q"))
    return g.mean() if how == "mean" else g.sum()

def log_linear_detrend(x, scale=100.0):
    log_x = np.log(np.asarray(x, dtype=float))
    t = np.arange(len(log_x))
    X = np.column_stack([np.ones_like(t, dtype=float), t.astype(float)])
    beta = np.linalg.lstsq(X, log_x, rcond=None)[0]
    return scale * (log_x - X @ beta)

def linear_detrend_level(x, recenter_to_mean=True):
    arr = np.asarray(x, dtype=float)
    t = np.arange(len(arr))
    X = np.column_stack([np.ones_like(t, dtype=float), t.astype(float)])
    beta = np.linalg.lstsq(X, arr, rcond=None)[0]
    resid = arr - X @ beta
    return resid + arr.mean() if recenter_to_mean else resid

def to_yq_frame(period_index, **cols):
    out = {"year": period_index.year.astype(int),
           "quarter": period_index.quarter.astype(int)}
    out.update(cols)
    return pd.DataFrame(out).reset_index(drop=True)

def clip_window(s, start=START_Q, end=END_Q):
    return s.loc[pd.Period(start, "Q"):pd.Period(end, "Q")]


#%%
# =============================================================================
# download all FRED inputs

print("Downloading FRED series...")
fred = {}
fred_ids = ["PCND", "PCESV", "PCDG", "GPDI", "GCE",
            "GDPDEF", "FEDFUNDS",
            "MVPHGFD027MNFRBDAL",
            "HOANBS", "B087RC1Q027SBEA"]
if POP_SOURCE == "FRED":
    fred_ids.append("CNP16OV")

for sid in fred_ids:
    fred[sid] = fred_get(sid)
    print(f"  {sid:25s}  {fred[sid].index.min().date()} -> {fred[sid].index.max().date()}  ({len(fred[sid])} obs)")
    time.sleep(0.1)


#%%
# =============================================================================
# convert FRED inputs to quarterly PeriodIndex

pcnd_q     = to_qper(fred["PCND"])
pcesv_q    = to_qper(fred["PCESV"])
pcdg_q     = to_qper(fred["PCDG"])
gpdi_q     = to_qper(fred["GPDI"])
gce_q      = to_qper(fred["GCE"])
gdpdef_q   = to_qper(fred["GDPDEF"])
hoanbs_q   = to_qper(fred["HOANBS"])
govtran_q  = to_qper(fred["B087RC1Q027SBEA"])

# Monthly -> quarterly. All three are stocks/rates -> use MEAN consistently.
# Note: for debt, sum and mean give numerically identical output after the
# log-linear-detrend (log(3*mean) = log(3) + log(mean); the log(3) is absorbed
# by the trend intercept). Mean is the conventional choice for stock and
# rate variables, and matches the population aggregation.
ffr_q   = monthly_to_quarterly(fred["FEDFUNDS"],            how="mean")  # rate
debt_q  = monthly_to_quarterly(fred["MVPHGFD027MNFRBDAL"],  how="mean")  # stock

if POP_SOURCE == "FRED":
    pop_q = monthly_to_quarterly(fred["CNP16OV"], how="mean").rename("POP_THOUS")
elif POP_SOURCE == "CSV":
    raw = pd.read_csv(POP_CSV_PATH)
    raw = raw[raw["Value"].astype(str).ne("-")].copy()
    raw["Value"]   = raw["Value"].astype(float)
    raw["quarter"] = raw["Period"].astype(str).str.extract(r"Q0?([1-4])").astype(int)
    raw["QPER"]    = pd.PeriodIndex.from_fields(year=raw["Year"].astype(int),
                                                quarter=raw["quarter"], freq="Q")
    pop_q = raw.set_index("QPER")["Value"].sort_index().rename("POP_THOUS")
else:
    raise ValueError(f"Unknown POP_SOURCE: {POP_SOURCE}")


#%%
# =============================================================================
# download BEA NIPA tables (1.12, 3.1, 3.2, 3.3)

print("\nDownloading BEA NIPA tables...")
nipa = {}
for tname in ["T11200", "T30100", "T30200", "T30300"]:
    nipa[tname] = bea_nipa(tname, freq="Q")
    rows = nipa[tname]["LineDescription"].nunique()
    span = (nipa[tname]["QPER"].min(), nipa[tname]["QPER"].max())
    print(f"  {tname}  {rows} lines  {span[0]} -> {span[1]}")
    time.sleep(0.5)


#%%
# =============================================================================
# build: chat   (consumption = PCND + PCESV)

def build_chat(start=START_Q, end=END_Q):
    nominal  = pcnd_q + pcesv_q
    real_pc  = (nominal / (gdpdef_q / 100.0)) / (pop_q * 1000.0)
    real_pc  = clip_window(real_pc.dropna(), start, end)
    return to_yq_frame(real_pc.index, chat=log_linear_detrend(real_pc.values))

df_c = build_chat()
print(df_c[["chat"]].describe().T)


#%%
# =============================================================================
# build: ihat   (investment = GPDI + PCDG)

def build_ihat(start=START_Q, end=END_Q):
    nominal  = gpdi_q + pcdg_q
    real_pc  = (nominal / (gdpdef_q / 100.0)) / (pop_q * 1000.0)
    real_pc  = clip_window(real_pc.dropna(), start, end)
    return to_yq_frame(real_pc.index, ihat=log_linear_detrend(real_pc.values))

df_i = build_ihat()
print(df_i[["ihat"]].describe().T)


#%%
# =============================================================================
# build: ghat   (government consumption + investment = GCE)

def build_ghat(start=START_Q, end=END_Q):
    real_pc = (gce_q / (gdpdef_q / 100.0)) / (pop_q * 1000.0)
    real_pc = clip_window(real_pc.dropna(), start, end)
    return to_yq_frame(real_pc.index, ghat=log_linear_detrend(real_pc.values))

df_g = build_ghat()
print(df_g[["ghat"]].describe().T)


#%%
# =============================================================================
# build: bhat   (federal debt held by public, market value, Dallas Fed)

def build_bhat(start=START_Q, end=END_Q):
    real_pc = (debt_q / (gdpdef_q / 100.0)) / (pop_q * 1000.0)
    real_pc = clip_window(real_pc.dropna(), start, end)
    return to_yq_frame(real_pc.index, bhat=log_linear_detrend(real_pc.values))

df_b = build_bhat()
print(df_b[["bhat"]].describe().T)


#%%
# =============================================================================
# build: pi   (gross inflation = GDPDEF[t] / GDPDEF[t-1])

def build_pi(start=START_Q, end=END_Q):
    p  = gdpdef_q.dropna()
    pi = clip_window((p / p.shift(1)).dropna(), start, end)
    return to_yq_frame(pi.index, pi=pi.values)

df_pi = build_pi()
print(df_pi[["pi"]].describe().T)


#%%
# =============================================================================
# build: r   (federal funds rate / 400, quarterly decimal)

def build_r(start=START_Q, end=END_Q):
    rr = clip_window(ffr_q.dropna(), start, end) / 400.0
    return to_yq_frame(rr.index, r=rr.values)

df_r = build_r()
print(df_r[["r"]].describe().T)


#%%
# =============================================================================
# build: lhat   (per-capita hours, HOANBS / population, log-detrended)

def build_lhat(start=START_Q, end=END_Q):
    hpc = hoanbs_q / (pop_q * 1000.0)
    hpc = clip_window(hpc.dropna(), start, end)
    return to_yq_frame(hpc.index, lhat=log_linear_detrend(hpc.values))

df_l = build_lhat()
print(df_l[["lhat"]].describe().T)


#%%
# =============================================================================
# build: trhat   (real per-capita government transfers)

def build_trhat(start=START_Q, end=END_Q):
    real_pc = (govtran_q / (gdpdef_q / 100.0)) / (pop_q * 1000.0)
    real_pc = clip_window(real_pc.dropna(), start, end)
    return to_yq_frame(real_pc.index, trhat=log_linear_detrend(real_pc.values))

df_tr = build_trhat()
print(df_tr[["trhat"]].describe().T)


#%%
# =============================================================================
# extract NIPA components for tax-rate construction (Jones 2002 method)
#
# Description-based matching is robust to BEA renumbering tables.

T112 = nipa["T11200"]
T31  = nipa["T30100"]
T32  = nipa["T30200"]
T33  = nipa["T30300"]


def find_line(df_table, pattern, prefer_line=None, exclude=()):
    """
    Find a NIPA line by description pattern. If multiple lines match, prefer
    one whose description matches `pattern` exactly (case-insensitive); if
    that fails, prefer the smallest line number; respect prefer_line if given.
    `exclude` is a tuple of substrings: any matching description is dropped.
    """
    cands = df_table[df_table["LineDescription"].str.contains(
        pattern, case=False, regex=False, na=False)]
    if exclude:
        for ex in exclude:
            cands = cands[~cands["LineDescription"].str.contains(
                ex, case=False, regex=False, na=False)]
    cands = cands[["LineNumber","LineDescription"]].drop_duplicates().sort_values("LineNumber")
    if cands.empty:
        raise ValueError(f"No NIPA line matched '{pattern}'")
    if prefer_line is not None and prefer_line in cands["LineNumber"].values:
        return int(prefer_line)
    # Prefer exact (case-insensitive) match on description
    exact = cands[cands["LineDescription"].str.lower() == pattern.lower()]
    if not exact.empty:
        return int(exact.iloc[0]["LineNumber"])
    # Otherwise smallest line number among matches (typically the headline aggregate)
    return int(cands.iloc[0]["LineNumber"])


def get_line(df_table, line_number, table_label=""):
    """Pull line by explicit number; print description for verification."""
    sub = df_table[df_table["LineNumber"] == line_number]
    if sub.empty:
        raise ValueError(f"{table_label} line {line_number} not found")
    desc = sub["LineDescription"].iloc[0]
    s = sub.set_index("QPER")["DataValue"].sort_index()
    s.name = f"{table_label}L{line_number}"
    print(f"  {table_label} L{line_number:>2d}: {desc}")
    return s


print("\nDiscovering NIPA line numbers in current BEA vintage...")
print("(Pattern -> auto-detected line number)")

# T11200 -- National Income by Type of Income
n_EC      = find_line(T112, "Compensation of employees", prefer_line=2)
n_W       = find_line(T112, "Wages and salaries",        prefer_line=3)
n_PRI     = find_line(T112, "Proprietors' income")
n_Rental  = find_line(T112, "Rental income of persons")
n_CP      = find_line(T112, "Corporate profits")
n_NI      = find_line(T112, "Net interest")

# T30100 -- Government Total Receipts and Expenditures
n_CSI     = find_line(T31,  "Contributions for government social insurance")
n_TPI     = find_line(T31,  "Taxes on production and imports")

# T30200 -- Federal Government
# "Personal current taxes" can also appear as a memo or sub-line; pick top-level
n_FIT     = find_line(T32,  "Personal current taxes",
                      exclude=("less:", "addenda", "less"))
n_CT_fed  = find_line(T32,  "Taxes on corporate income",
                      exclude=("less:", "addenda"))

# T30300 -- State and Local Government
n_SIT     = find_line(T33,  "Personal current taxes",
                      exclude=("less:", "addenda"))
n_CT_sl   = find_line(T33,  "Taxes on corporate income",
                      exclude=("less:", "addenda"))
n_PT      = find_line(T33,  "Property taxes",
                      exclude=("less:", "addenda"))

print("\nExtracting series with discovered line numbers:")
print("  T11200 (National Income by Type):")
EC      = get_line(T112, n_EC,      "T11200")
W       = get_line(T112, n_W,       "T11200")
PRI     = get_line(T112, n_PRI,     "T11200")
Rental  = get_line(T112, n_Rental,  "T11200")
CP      = get_line(T112, n_CP,      "T11200")
NI      = get_line(T112, n_NI,      "T11200")

print("  T30100 (Government Total):")
CSI     = get_line(T31,  n_CSI,     "T30100")
TPI     = get_line(T31,  n_TPI,     "T30100")

print("  T30200 (Federal):")
FIT     = get_line(T32,  n_FIT,     "T30200")
CT_fed  = get_line(T32,  n_CT_fed,  "T30200")

print("  T30300 (State and Local):")
SIT     = get_line(T33,  n_SIT,     "T30300")
CT_sl   = get_line(T33,  n_CT_sl,   "T30300")
PT      = get_line(T33,  n_PT,      "T30300")

# Sanity check: most recent quarter values should be in expected ranges.
# All values are SAAR in $bn.
expected = {
    "TPI":    (1500, 2300, "broad indirect taxes incl. sales/excise/customs"),
    "PT":     (500,   900, "state/local property taxes"),
    "FIT":    (1500, 3000, "federal personal income taxes"),
    "SIT":    (300,   700, "state/local personal income taxes"),
    "CT_fed": (200,   700, "federal corporate income taxes"),
    "CT_sl":  (40,    150, "state/local corporate income taxes"),
    "CSI":    (1300, 2100, "social insurance contributions"),
    "EC":     (10000, 17500, "compensation of employees"),
    "W":      (8000, 13500, "wages and salaries"),
    "PRI":    (1500, 2700, "proprietors' income"),
    "CP":     (1500, 4500, "corporate profits"),
    "NI":     (200,  1200, "net interest"),
    "Rental": (200,  1200, "rental income"),
}
print("\nSanity check (most recent quarter values, SAAR $bn):")
print(f"  {'series':8s} {'value':>10s}  {'expected':>15s}  status  description")
locals_map = {"TPI":TPI,"PT":PT,"FIT":FIT,"SIT":SIT,"CT_fed":CT_fed,
              "CT_sl":CT_sl,"CSI":CSI,"EC":EC,"W":W,"PRI":PRI,"CP":CP,
              "NI":NI,"Rental":Rental}
for name, s in locals_map.items():
    last = s.dropna()
    if len(last) == 0:
        print(f"  {name:8s} {'no data':>10s}")
        continue
    val = last.iloc[-1]
    lo, hi, desc = expected[name]
    in_range = lo <= val <= hi
    flag = "  OK  " if in_range else "  !!  "
    print(f"  {name:8s} {val:>10,.1f}  [{lo:>5,.0f}, {hi:>5,.0f}]  {flag}  {desc}")

# personal income tax rate
CI    = Rental + CP + NI + PRI / 2.0
tau_p = (FIT + SIT) / (W + PRI / 2.0 + CI)
tau_p = clip_window(tau_p.dropna(), START_Q, END_Q)
print(f"\ntau_p:   mean={tau_p.mean():.4f}  std={tau_p.std():.4f}  (Jones 2002 reports ~0.10-0.15)")


#%%
# =============================================================================
# build: tau_l hats (labor tax rate, Jones 2002)
#
# Produces both detrended and demean-only fractional-deviation variants.
# 3_clean_for_matlab.py picks one as the canonical tau_l_hat for the .mod file.

def build_tau_l(start=START_Q, end=END_Q):
    tl = (tau_p * (W + PRI / 2.0) + CSI) / (EC + PRI / 2.0)
    tl = clip_window(tl.dropna(), start, end)
    mu = tl.mean()
    # Linear-detrended fractional deviation
    tl_dt = pd.Series(linear_detrend_level(tl.values), index=tl.index)
    tl_hat_dt = (tl_dt - mu) / mu
    # Demean-only fractional deviation
    tl_hat_nd = (tl - mu) / mu
    return to_yq_frame(tl.index,
                       tau_l=tl.values,
                       tau_l_detrend=tl_dt.values,
                       tau_l_hat_detrend=tl_hat_dt.values,
                       tau_l_hat_nodetrend=tl_hat_nd.values)

df_tl = build_tau_l()
print(df_tl[["tau_l","tau_l_hat_detrend","tau_l_hat_nodetrend"]].describe().T)


#%%
# =============================================================================
# build: tau_k hats (capital tax rate, Jones 2002)
#
# Produces both detrended and demean-only fractional-deviation variants.
# 3_clean_for_matlab.py picks one as the canonical tau_k_hat for the .mod file.

def build_tau_k(start=START_Q, end=END_Q):
    CT = CT_fed + CT_sl
    # Mendoza-Razin-Tesar (1994):
    #   numerator   = total capital tax revenue
    #   denominator = capital income plus corporate and property taxes
    # tk = (tau_p * CI + CT + PT) / (CI + CT + PT)  # lower alternative
    # Jones (2002):
    #   numerator   = total capital tax revenue
    #   denominator = measured capital income
    tk = (tau_p * CI + CT + PT) / CI
    tk = clip_window(tk.dropna(), start, end)
    mu = tk.mean()
    # Linear-detrended fractional deviation
    tk_dt = pd.Series(linear_detrend_level(tk.values), index=tk.index)
    tk_hat_dt = (tk_dt - mu) / mu
    # Demean-only fractional deviation
    tk_hat_nd = (tk - mu) / mu
    return to_yq_frame(tk.index,
                       tau_k=tk.values,
                       tau_k_detrend=tk_dt.values,
                       tau_k_hat_detrend=tk_hat_dt.values,
                       tau_k_hat_nodetrend=tk_hat_nd.values)

df_tk = build_tau_k()
print(df_tk[["tau_k","tau_k_hat_detrend","tau_k_hat_nodetrend"]].describe().T)


#%%
# =============================================================================
# build: tau_c hats (consumption tax rate, Mendoza-Razin-Tesar 1994)
#
# Produces both detrended and demean-only fractional-deviation variants.
# 3_clean_for_matlab.py picks one as the canonical tau_c_hat for the .mod file.

def build_tau_c(start=START_Q, end=END_Q):
    pce_total = pcnd_q + pcesv_q + pcdg_q
    indirect  = TPI - PT
    tc = indirect / (pce_total + gce_q - indirect)
    tc = clip_window(tc.dropna(), start, end)
    mu = tc.mean()
    # Linear-detrended fractional deviation
    tc_dt = pd.Series(linear_detrend_level(tc.values), index=tc.index)
    tc_hat_dt = (tc_dt - mu) / mu
    # Demean-only fractional deviation (MRT 1994 standard)
    tc_hat_nd = (tc - mu) / mu
    return to_yq_frame(tc.index,
                       tau_c=tc.values,
                       tau_c_detrend=tc_dt.values,
                       tau_c_hat_detrend=tc_hat_dt.values,
                       tau_c_hat_nodetrend=tc_hat_nd.values)

df_tc = build_tau_c()
print(df_tc[["tau_c","tau_c_hat_detrend","tau_c_hat_nodetrend"]].describe().T)


#%%
# =============================================================================
# build: tau_d hats (dividend tax rate, proxied by personal income tax rate)
#
# No clean quarterly source. Pre-2003 dividends taxed as ordinary income (so
# tau_p is exact); post-JGTRRA preferential rates make tau_p an upper bound.
# Swap in McGrattan-Prescott or Gourio-Miao series here when available.
#
# Produces both detrended and demean-only fractional-deviation variants.
# 3_clean_for_matlab.py picks one as the canonical tau_d_hat for the .mod file.

def build_tau_d(start=START_Q, end=END_Q):
    td = clip_window(tau_p, start, end)
    mu = td.mean()
    # Linear-detrended fractional deviation
    td_dt = pd.Series(linear_detrend_level(td.values), index=td.index)
    td_hat_dt = (td_dt - mu) / mu
    # Demean-only fractional deviation
    td_hat_nd = (td - mu) / mu
    return to_yq_frame(td.index,
                       tau_d=td.values,
                       tau_d_detrend=td_dt.values,
                       tau_d_hat_detrend=td_hat_dt.values,
                       tau_d_hat_nodetrend=td_hat_nd.values)

df_td = build_tau_d()
print(df_td[["tau_d","tau_d_hat_detrend","tau_d_hat_nodetrend"]].describe().T)


#%%
# =============================================================================
# build: tau_itc_hat (aggregate Investment Tax Credit rate)
#
# The federal ITC existed 1962-1986. House-Shapiro's ITCdat.txt provides
# quarterly ITC rates by asset type for 1959Q1-2007Q1; rates are zero
# pre-1962 (ITC didn't exist) and post-1986 (repealed by TRA 1986).
#
# Coverage handling:
#   1958Q1-1958Q4:    ITC = 0 (pre-existence, backfilled)
#   1959Q1-2007Q1:    HS data (zero in periods when ITC didn't apply)
#   2007Q2-current:   ITC = 0 (statutory, repealed)
#
# Aggregation uses HS's own PQdat.txt for weights, since HS's 38-type ITC
# rates align with HS's 38-type investment shares (BEA's current taxonomy
# in PQdat_BEA.csv has ~70 different categories).

def build_tau_ic(itc_path, pq_path, start_q=START_Q, end_q=END_Q):
    df_itc = pd.read_csv(itc_path, sep=r"\s+", header=None)
    df_itc.columns = [f"itc_{i+1:02d}" for i in range(df_itc.shape[1])]
    df_inv = pd.read_csv(pq_path, sep=r"\s+", header=None)
    df_inv.columns = [f"inv_{i+1:02d}" for i in range(df_inv.shape[1])]

    q0 = pd.Period("1959Q1", freq="Q")
    df_itc["QPER"] = [q0 + i for i in range(len(df_itc))]
    df_inv["QPER"] = [q0 + i for i in range(len(df_inv))]

    raw_start = max(df_itc["QPER"].min(), df_inv["QPER"].min())  # 1959Q1
    raw_end   = min(df_itc["QPER"].max(), df_inv["QPER"].max())  # 2007Q1
    idx = pd.period_range(start_q, end_q, freq="Q")

    df_itc = df_itc.set_index("QPER").reindex(idx)
    df_inv = df_inv.set_index("QPER").reindex(idx)

    # Pre-1959Q1 (1958): ITC = 0 (didn't exist), use first observed weights
    pre_mask = idx < raw_start
    if pre_mask.any():
        df_itc.loc[pre_mask] = 0.0
        # Backfill weights from 1959Q1 (first observed) -- proxy for 1958
        first_obs = df_inv.loc[raw_start]
        df_inv.loc[pre_mask] = first_obs.values

    # Post-2007Q1: ITC = 0 (repealed), forward-fill weights
    post_mask = idx > raw_end
    df_itc.loc[post_mask] = df_itc.loc[post_mask].fillna(0.0)
    df_inv = df_inv.ffill()

    shares = df_inv.div(df_inv.sum(axis=1), axis=0)
    shares.columns = df_itc.columns
    tau_ic = (df_itc * shares).sum(axis=1, min_count=1)

    mu_ic = tau_ic.dropna().mean()
    tau_ic_hat = tau_ic.to_numpy() - mu_ic

    return to_yq_frame(idx, tau_itc=tau_ic.to_numpy(),
                       tau_itc_hat=tau_ic_hat)

df_itc = build_tau_ic(itc_path=os.path.join(RAW_DIR, "ITCdat.txt"),
                      pq_path =os.path.join(RAW_DIR, "PQdat.txt"))
print(df_itc[["tau_itc","tau_itc_hat"]].describe().T)


#%%
# =============================================================================
# build: e_tau_hat   (aggregate bonus-depreciation expensing share)
#
# Bonus depreciation eligibility: Section 168(k) covers MACRS recovery
# period <= 20 years (equipment, software, R&D, entertainment originals,
# utility quasi-structures). 39-year structures and 27.5-year residential
# are not eligible.
#
# Weights from raw/PQdat_BEA.csv (BEA Fixed Assets Table 2.7, annual values
# applied to all four quarters). Eligibility flags from raw/treat_BEA.csv.
# Both produced by 0_pre_prep.py.
#
# Statutory bonus calendar:
#   Pre-2001Q4:    0     (no bonus regime)
#   2001Q4-2003Q2: 30%   (JCWAA 2002, retroactive to 9/11/2001)
#   2003Q3-2004Q4: 50%   (JGTRRA 2003)
#   2005Q1-2007Q4: 0     (bonus expired)
#   2008Q1-2010Q4: 50%   (Economic Stimulus Act 2008 + ARRA 2009 + SBJA 2010)
#   2011Q1-2011Q4: 100%  (Tax Relief Act 2010)
#   2012Q1-2013Q4: 50%   (American Taxpayer Relief Act 2012)
#   2014Q1-2014Q4: 50%   (Tax Increase Prevention Act 2014, retroactive)
#   2015Q1-2017Q3: 50%   (PATH Act 2015)
#   2017Q4-2022Q4: 100%  (TCJA, property placed in service after 9/27/2017)
#   2023Q1-2023Q4: 80%   (TCJA phase-down)
#   2024Q1-2024Q4: 60%
#   2025Q1:        40%   (most of Q1 pre-1/19/2025 OBBBA effective date)
#   2025Q2-2026+:  100%  (OBBBA, July 2025, retroactive to 1/20/2025, permanent)

def build_e_tau(pq_path, treat_path, start_q=START_Q, end_q=END_Q):
    # PQdat_BEA.csv: index=QPER (quarterly period strings), columns=asset types
    df_inv = pd.read_csv(pq_path, index_col="QPER")
    df_inv.index = pd.PeriodIndex(df_inv.index, freq="Q")

    # treat_BEA.csv: single row, columns=asset types matching df_inv
    treat = pd.read_csv(treat_path).iloc[0]
    treat = treat.reindex(df_inv.columns).fillna(0).astype(int)

    raw_start = df_inv.index.min()
    raw_end   = df_inv.index.max()
    idx = pd.period_range(start_q, end_q, freq="Q")

    # Reindex to full period; forward-fill if BEA data is shorter than END_Q
    df_inv = df_inv.reindex(idx)
    df_inv = df_inv.ffill()
    if raw_start > idx.min():
        df_inv.loc[idx < raw_start] = np.nan

    shares = df_inv.div(df_inv.sum(axis=1), axis=0)
    eligible_share = shares.mul(treat, axis=1).sum(axis=1, min_count=1)

    bonus = pd.Series(0.0, index=idx)
    bonus.loc[pd.period_range("2001Q4", "2003Q2", freq="Q")] = 0.30
    bonus.loc[pd.period_range("2003Q3", "2004Q4", freq="Q")] = 0.50
    bonus.loc[pd.period_range("2008Q1", "2010Q4", freq="Q")] = 0.50
    bonus.loc[pd.period_range("2011Q1", "2011Q4", freq="Q")] = 1.00
    bonus.loc[pd.period_range("2012Q1", "2013Q4", freq="Q")] = 0.50
    bonus.loc[pd.period_range("2014Q1", "2017Q3", freq="Q")] = 0.50
    bonus.loc[pd.period_range("2017Q4", "2022Q4", freq="Q")] = 1.00
    bonus.loc[pd.period_range("2023Q1", "2023Q4", freq="Q")] = 0.80
    bonus.loc[pd.period_range("2024Q1", "2024Q4", freq="Q")] = 0.60
    bonus.loc[pd.Period("2025Q1", "Q")]                       = 0.40
    bonus.loc[pd.period_range("2025Q2", end_q,    freq="Q")] = 1.00

    e_tau = bonus * eligible_share
    mu_e  = e_tau.dropna().mean()
    e_tau_hat = e_tau.to_numpy() - mu_e

    return to_yq_frame(idx, e_tau=e_tau.to_numpy(), e_tau_hat=e_tau_hat)

df_e = build_e_tau(pq_path   =os.path.join(RAW_DIR, "PQdat_BEA.csv"),
                   treat_path=os.path.join(RAW_DIR, "treat_BEA.csv"))
print(df_e[["e_tau","e_tau_hat"]].describe().T)


#%%
# =============================================================================
# build: tau_i_hat (interest income tax) and tau_div_hat (dividend tax)
#
# Source: McGrattan-style TAXSIM federal marginal tax rates, computed by
# Feenberg-Coutts NBER TAXSIM model, distributed by McGrattan as
# raw/taxsim_fed_McGrattan.txt. Annual coverage 1960-2013.
#
# Columns used:
#   "Interest"           -> tau_i  (interest income tax rate)
#   "Qualified Dividends"-> tau_div (dividend income tax rate; pre-2003
#                                    this equals the ordinary dividend
#                                    column, post-2003 it captures the
#                                    JGTRRA preferential rate)
#
# Coverage outside 1960-2013:
#   1958Q1-1959Q4: NaN (TAXSIM data starts 1960)
#   2014Q1-current: NaN (McGrattan file ends 2013)
#
# Annual values are applied to all four quarters of each year. Both
# fractional-deviation variants are produced (parallel to tau_l/tau_k/
# tau_c/tau_d): linear-detrended (_hat_detrend) and demean-only
# (_hat_nodetrend). The linear trend is fit on observed quarters only,
# leaving NaN outside the data window.
#
# These series are NOT used in the .mod file's varobs (because of high
# collinearity with tau_l_hat), but they are exported for data-completeness
# and any downstream robustness analysis.

def build_taxsim_rates(path, start_q=START_Q, end_q=END_Q):
    rows = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            year = int(parts[0])
            interest      = float(parts[2]) if parts[2] != "NaN" else np.nan
            qualified_div = float(parts[4]) if parts[4] != "NaN" else np.nan
            rows.append((year, interest, qualified_div))

    df_ann = pd.DataFrame(rows, columns=["YEAR", "tau_i", "tau_div"])
    df_ann = df_ann.set_index("YEAR") / 100.0  # percent -> decimal

    idx = pd.period_range(start_q, end_q, freq="Q")
    tau_i_q   = pd.Series(np.nan, index=idx)
    tau_div_q = pd.Series(np.nan, index=idx)
    for q in idx:
        y = q.year
        if y in df_ann.index:
            tau_i_q.loc[q]   = df_ann.loc[y, "tau_i"]
            tau_div_q.loc[q] = df_ann.loc[y, "tau_div"]

    def detrend_with_gaps(s):
        """Linear-detrend; fits the trend on observed (non-NaN) quarters
        only and returns NaN outside that window. Recentered to observed
        mean to match the convention used by linear_detrend_level."""
        out = pd.Series(np.nan, index=s.index)
        obs_mask = ~s.isna()
        if obs_mask.sum() < 2:
            return out
        obs_pos = np.where(obs_mask.values)[0].astype(float)
        obs_val = s.values[obs_mask.values]
        X = np.column_stack([np.ones_like(obs_pos), obs_pos])
        beta = np.linalg.lstsq(X, obs_val, rcond=None)[0]
        resid = obs_val - X @ beta
        out.iloc[obs_mask.values] = resid + obs_val.mean()
        return out

    # Demean using only observed quarters
    mu_i = tau_i_q.dropna().mean()
    mu_d = tau_div_q.dropna().mean()

    # Two variants each (parallel to tau_l/tau_k/tau_c/tau_d):
    # detrend = linear-detrended fractional deviation
    # nodetrend = demean-only fractional deviation
    tau_i_dt   = detrend_with_gaps(tau_i_q)
    tau_div_dt = detrend_with_gaps(tau_div_q)
    tau_i_hat_detrend     = (tau_i_dt   - mu_i) / mu_i
    tau_i_hat_nodetrend   = (tau_i_q    - mu_i) / mu_i
    tau_div_hat_detrend   = (tau_div_dt - mu_d) / mu_d
    tau_div_hat_nodetrend = (tau_div_q  - mu_d) / mu_d

    return to_yq_frame(idx,
                       tau_i                 = tau_i_q.to_numpy(),
                       tau_i_detrend         = tau_i_dt.to_numpy(),
                       tau_i_hat_detrend     = tau_i_hat_detrend.to_numpy(),
                       tau_i_hat_nodetrend   = tau_i_hat_nodetrend.to_numpy(),
                       tau_div               = tau_div_q.to_numpy(),
                       tau_div_detrend       = tau_div_dt.to_numpy(),
                       tau_div_hat_detrend   = tau_div_hat_detrend.to_numpy(),
                       tau_div_hat_nodetrend = tau_div_hat_nodetrend.to_numpy())

df_taxsim = build_taxsim_rates(os.path.join(RAW_DIR, "taxsim_fed_McGrattan.txt"))
print(df_taxsim[["tau_i","tau_i_hat_detrend","tau_i_hat_nodetrend",
                 "tau_div","tau_div_hat_detrend","tau_div_hat_nodetrend"]].describe().T)


#%%
# =============================================================================
# merge all series and export US_Data_Full.xlsx

frames = [df_c, df_i, df_g, df_b, df_pi, df_r,
          df_l, df_tr,
          df_tl, df_tk, df_tc, df_td,
          df_itc, df_e, df_taxsim]

df_out = frames[0]
for f in frames[1:]:
    df_out = df_out.merge(f, on=["year", "quarter"], how="outer")
df_out = df_out.sort_values(["year", "quarter"]).reset_index(drop=True)

# Final column order: identifiers first
id_cols  = ["year", "quarter"]
obs_cols = ["chat", "ihat", "ghat", "bhat", "pi", "r",
            "lhat", "trhat",
            "tau_l", "tau_l_detrend", "tau_l_hat_detrend", "tau_l_hat_nodetrend",
            "tau_k", "tau_k_detrend", "tau_k_hat_detrend", "tau_k_hat_nodetrend",
            "tau_c", "tau_c_detrend", "tau_c_hat_detrend", "tau_c_hat_nodetrend",
            "tau_d", "tau_d_detrend", "tau_d_hat_detrend", "tau_d_hat_nodetrend",
            "tau_itc", "tau_itc_hat",
            "e_tau", "e_tau_hat",
            "tau_i", "tau_i_detrend", "tau_i_hat_detrend", "tau_i_hat_nodetrend",
            "tau_div", "tau_div_detrend", "tau_div_hat_detrend", "tau_div_hat_nodetrend"]
df_out = df_out[id_cols + [c for c in obs_cols if c in df_out.columns]]

print("\nFinal series summary:")
print(df_out.drop(columns=id_cols).agg(["mean","std","min","max"]).T)

out_path = os.path.join(OUT_DIR, "US_Data_Full.xlsx")
df_out.to_excel(out_path, index=False)
print(f"\nWrote {out_path}  ({len(df_out)} rows, {df_out.shape[1]} cols)")
