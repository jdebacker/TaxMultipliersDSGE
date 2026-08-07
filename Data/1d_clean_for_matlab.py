"""
3_clean_for_matlab.py
=====================
Prepare two MATLAB-ready datasets from US_Data_Full.xlsx:

  US_Data_Matlab.xlsx          -- 1958Q1 to 2019Q4 (extended sample)
  US_Data_Matlab_Zubairy.xlsx  -- 1958Q1 to 2008Q4 (Zubairy's sample window)

For each detrendable tax rate (tau_l, tau_k, tau_c, tau_d), 1_build_data.py
produces both `tau_x_hat_detrend` and `tau_x_hat_nodetrend`. The HAT_CHOICE
dict below picks one variant per tax to be renamed to the canonical
`tau_x_hat` for the .mod file. Switching specifications is a one-line edit.
"""

#%%
# =============================================================================
# imports and HAT_CHOICE selector

import os
import pandas as pd
import numpy as np
from scipy.io import savemat

OUT_DIR = r"./output"
SRC = os.path.join(OUT_DIR, "US_Data_Full.xlsx")

# Pick which variant becomes the canonical tau_x_hat for the .mod file.
# Set each value to "detrend" or "nodetrend".
HAT_CHOICE = {
    "tau_l":   "detrend",     # linear-detrended (strong, smooth trend)
    "tau_k":   "nodetrend",   # demean-only (cyclical/policy variation dominates)
    "tau_c":   "detrend",     # linear-detrended (strong, smooth trend; ADF + magnitude agree)
    "tau_d":   "detrend",     # linear-detrended (strong, smooth trend)
    "tau_i":   "nodetrend",   # demean-only weak/cyclical
    "tau_div": "detrend",     # linear-detrended strong trend
}


#%%
# =============================================================================
# load full dataset and build PeriodIndex (compatible with older pandas)

df = pd.read_excel(SRC)

df["QPER"] = pd.PeriodIndex(
    [pd.Period(year=int(y), quarter=int(q), freq="Q")
     for y, q in zip(df["year"], df["quarter"])]
)
# Drop the original year/quarter columns; they get re-added from the index
# at the end of each window-export cell.
df = df.drop(columns=["year", "quarter"]).set_index("QPER").sort_index()

print(f"Loaded {SRC}")
print(f"  rows: {len(df)}  cols: {df.shape[1]}")
print(f"  range: {df.index.min()} -> {df.index.max()}")
print(f"\nHAT_CHOICE:")
for tax, choice in HAT_CHOICE.items():
    print(f"  {tax}_hat <- {tax}_hat_{choice}")


#%%
# =============================================================================
# build column-drop list and rename map from HAT_CHOICE

# Always-drop: level series and the unused hat variant for each tax.
# Year and quarter are kept (set as a date index for slicing, then
# re-added as columns at the front).
drop_cols = []
for tax, choice in HAT_CHOICE.items():
    drop_cols += [tax, f"{tax}_detrend"]
    other = "nodetrend" if choice == "detrend" else "detrend"
    drop_cols.append(f"{tax}_hat_{other}")
drop_cols += ["tau_itc", "e_tau"]

# Rename the chosen variant for each tax to the canonical tau_x_hat,
# and rename pi -> pi_obs and r -> r_obs to match the .mod file's varobs.
rename_map = {f"{tax}_hat_{choice}": f"{tax}_hat"
              for tax, choice in HAT_CHOICE.items()}
rename_map.update({"pi": "pi_obs", "r": "r_obs"})


#%%
# =============================================================================
# build US_Data_Matlab.xlsx -- 1958Q1 to 2019Q4

WINDOW_FULL = (pd.Period("1958Q1","Q"), pd.Period("2019Q4","Q"))

df_matlab = df.loc[WINDOW_FULL[0]:WINDOW_FULL[1]].copy()
df_matlab = df_matlab.drop(columns=[c for c in drop_cols if c in df_matlab.columns])
df_matlab = df_matlab.rename(columns=rename_map)

# Add year/quarter back as plain columns at the front so MATLAB can read them
df_matlab.insert(0, "quarter", df_matlab.index.quarter)
df_matlab.insert(0, "year",    df_matlab.index.year)
df_matlab = df_matlab.reset_index(drop=True)

out1 = os.path.join(OUT_DIR, "US_Data_Matlab.xlsx")
df_matlab.to_excel(out1, index=False)
print(f"\nWrote {out1}")
print(f"  rows: {len(df_matlab)}  cols: {df_matlab.shape[1]}")
print(f"  range: {WINDOW_FULL[0]} -> {WINDOW_FULL[1]}")
print(f"  columns: {df_matlab.columns.tolist()}")

# Also write a .mat file: each column becomes a column-vector variable.
# MATLAB-friendly variable names (no special chars). Loaded via load(),
# producing workspace variables matching the column names.
out1_mat = os.path.join(OUT_DIR, "US_Data_Matlab.mat")
mat_dict_full = {col: df_matlab[col].to_numpy().reshape(-1, 1)
                 for col in df_matlab.columns}
savemat(out1_mat, mat_dict_full, do_compression=True)
print(f"Wrote {out1_mat}")


#%%
# =============================================================================
# build US_Data_Matlab_Zubairy.xlsx -- 1958Q1 to 2008Q4

WINDOW_ZUB = (pd.Period("1958Q1","Q"), pd.Period("2008Q4","Q"))

df_zub = df.loc[WINDOW_ZUB[0]:WINDOW_ZUB[1]].copy()
df_zub = df_zub.drop(columns=[c for c in drop_cols if c in df_zub.columns])
df_zub = df_zub.rename(columns=rename_map)

df_zub.insert(0, "quarter", df_zub.index.quarter)
df_zub.insert(0, "year",    df_zub.index.year)
df_zub = df_zub.reset_index(drop=True)

out2 = os.path.join(OUT_DIR, "US_Data_Matlab_Zubairy.xlsx")
df_zub.to_excel(out2, index=False)
print(f"\nWrote {out2}")
print(f"  rows: {len(df_zub)}  cols: {df_zub.shape[1]}")
print(f"  range: {WINDOW_ZUB[0]} -> {WINDOW_ZUB[1]}")
print(f"  columns: {df_zub.columns.tolist()}")

out2_mat = os.path.join(OUT_DIR, "US_Data_Matlab_Zubairy.mat")
mat_dict_zub = {col: df_zub[col].to_numpy().reshape(-1, 1)
                for col in df_zub.columns}
savemat(out2_mat, mat_dict_zub, do_compression=True)
print(f"Wrote {out2_mat}")


#%%
# =============================================================================
# build US_Data_Matlab_1961-2013.xlsx -- TAXSIM-based dividend tax variant
#
# Coverage matches the McGrattan-TAXSIM file (1960-2013), trimmed to start
# at 1961Q1 because 1960's Interest column is NaN in the source file.
# The proxy tau_d (= tau_p) is dropped and tau_div (TAXSIM-based) is
# renamed to tau_d so the .mod file's existing varobs reference works
# unchanged. tau_i is kept in the dataset for completeness.

WINDOW_TAXSIM = (pd.Period("1961Q1","Q"), pd.Period("2013Q4","Q"))

df_tx = df.loc[WINDOW_TAXSIM[0]:WINDOW_TAXSIM[1]].copy()
df_tx = df_tx.drop(columns=[c for c in drop_cols if c in df_tx.columns])
df_tx = df_tx.rename(columns=rename_map)

# Replace tau_d (proxy) with tau_div (TAXSIM-based) under the canonical name
if "tau_d_hat" in df_tx.columns:
    df_tx = df_tx.drop(columns=["tau_d_hat"])
if "tau_div_hat" in df_tx.columns:
    df_tx = df_tx.rename(columns={"tau_div_hat": "tau_d_hat"})

df_tx.insert(0, "quarter", df_tx.index.quarter)
df_tx.insert(0, "year",    df_tx.index.year)
df_tx = df_tx.reset_index(drop=True)

out3 = os.path.join(OUT_DIR, "US_Data_Matlab_1961_2013.xlsx")
df_tx.to_excel(out3, index=False)
print(f"\nWrote {out3}")
print(f"  rows: {len(df_tx)}  cols: {df_tx.shape[1]}")
print(f"  range: {WINDOW_TAXSIM[0]} -> {WINDOW_TAXSIM[1]}")
print(f"  columns: {df_tx.columns.tolist()}")

out3_mat = os.path.join(OUT_DIR, "US_Data_Matlab_1961_2013.mat")
mat_dict_tx = {col: df_tx[col].to_numpy().reshape(-1, 1)
               for col in df_tx.columns}
savemat(out3_mat, mat_dict_tx, do_compression=True)
print(f"Wrote {out3_mat}")


#%%
# =============================================================================
# Print all calibration anchors over both windows

import requests, time


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


def fred_q(series_id, fred_key):
    """Pull a FRED series and convert to PeriodIndex (Q). Monthly -> qtly mean."""
    url = "https://api.stlouisfed.org/fred/series/observations"
    params = {"series_id": series_id, "api_key": fred_key,
              "file_type": "json", "observation_start": "1947-01-01"}
    for attempt in range(4):
        try:
            r = requests.get(url, params=params, timeout=120)
            r.raise_for_status()
            obs = r.json().get("observations", [])
            d = pd.DataFrame(obs)
            d["date"]  = pd.to_datetime(d["date"])
            d["value"] = pd.to_numeric(d["value"], errors="coerce")
            s = d.set_index("date")["value"].dropna()
            if (s.index[1] - s.index[0]).days < 60:
                return s.groupby(s.index.to_period("Q")).mean()
            s.index = s.index.to_period("Q")
            return s
        except (requests.ConnectionError, requests.Timeout):
            if attempt < 3:
                time.sleep(2.0 ** attempt)
            else:
                raise


def report_window(label, window, src_full, b_y, g_y):
    a, b = window
    sub = src_full.loc[a:b]
    print(f"\n{label}:")
    for v in ["tau_l", "tau_k", "tau_c", "tau_d", "tau_itc", "e_tau",
              "tau_i", "tau_div"]:
        if v in sub.columns:
            print(f"  {v}_bar = {sub[v].mean():.6f}")
    by = b_y.loc[a:b].dropna().mean()
    gy = g_y.loc[a:b].dropna().mean()
    print(f"  b_y      = {by:.4f}")
    print(f"  g_y      = {gy:.4f}")


# --- pull FRED inputs
FRED_KEY = _load_keys("./raw/api_keys.txt")["FRED_KEY"]
print("\nFetching GDP, GCE, debt from FRED for b/y and g/y...")
gdp  = fred_q("GDP",                FRED_KEY)   # billions $, SAAR
gce  = fred_q("GCE",                FRED_KEY)   # billions $, SAAR
debt = fred_q("MVPHGFD027MNFRBDAL", FRED_KEY)   # billions $, monthly -> qtly mean

# Print raw magnitudes so the unit is verifiable, not guessed
sample_q = pd.Period("2007Q4", "Q")
print(f"\nUnit check at {sample_q}:")
print(f"  GDP   = {gdp.loc[sample_q]:>12,.1f}")
print(f"  GCE   = {gce.loc[sample_q]:>12,.1f}")
print(f"  debt  = {debt.loc[sample_q]:>12,.1f}")

# Debt and GDP are both in billions, so this is directly debt/GDP
b_y = debt / gdp
g_y = gce / gdp

# --- reload level data
src_full = pd.read_excel(SRC)
src_full["QPER"] = pd.PeriodIndex(
    [pd.Period(year=int(y), quarter=int(q), freq="Q")
     for y, q in zip(src_full["year"], src_full["quarter"])]
)
src_full = src_full.set_index("QPER")

print("\n" + "="*70)
print("Steady-state calibration values for tax_dsge_est.mod")
print("="*70)

report_window("Full sample (1958Q1-2019Q4)",    WINDOW_FULL,   src_full, b_y, g_y)
report_window("Zubairy window (1958Q1-2008Q4)", WINDOW_ZUB,    src_full, b_y, g_y)
report_window("TAXSIM window (1961Q1-2013Q4)",  WINDOW_TAXSIM, src_full, b_y, g_y)