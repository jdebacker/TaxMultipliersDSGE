"""
2_diagnose_data.py
==================
Sanity checks on the constructed series in US_Data_Full.xlsx:

  1. Verify tax-rate formulas produce values in [0, 1]
  2. Sample means of tax rates over Zubairy's window (literature reference)
  3. ADF stationarity tests on the hat series
  4. Plot real-quantity hats so they look stationary around 0
  5. Plot tax-rate levels with linear-trend overlays
  6. Linear-trend slope-significance test
  6b. Economic-magnitude trend assessment + structural-break flag
  7. Diagnostic context: how to use these results to inform HAT_CHOICE
  8. Debt-to-GDP convention check (b/y as years of annual output)

This script DIAGNOSES whether each tax rate has a meaningful trend; it
does NOT decide. The detrend-vs-demean choice for each tax variable is
made in 3_clean_for_matlab.py via HAT_CHOICE.

Run AFTER 1_build_data.py finishes. Reads US_Data_Full.xlsx from OUT_DIR.
"""

#%%
# =============================================================================
# imports

import os
import time
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import adfuller
import scipy.stats as stats

# Paths/keys must match 1_build_data.py
OUT_DIR  = r"./output"
RAW_DIR  = r"./raw"
START_Q  = "1958Q1"
END_Q    = "2026Q1"

def _load_api_keys(path):
    keys = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            keys[k.strip()] = v.strip().strip('"').strip("'")
    return keys

_keys = _load_api_keys(os.path.join(RAW_DIR, "api_keys.txt"))
FRED_KEY = _keys["FRED_KEY"]


df = pd.read_excel(os.path.join(OUT_DIR, "US_Data_Full.xlsx"))

# Build PeriodIndex compatibly across pandas versions
# (pd.PeriodIndex.from_fields requires pandas >= 2.2)
df["QPER"] = pd.PeriodIndex(
    [pd.Period(year=int(y), quarter=int(q), freq="Q")
     for y, q in zip(df["year"], df["quarter"])]
)
df = df.set_index("QPER")

# Restrict to Zubairy's window for direct comparison to her paper
zub = df.loc[pd.Period("1958Q1","Q"):pd.Period("2008Q4","Q")]


#%%
# =============================================================================
# CHECK 1: tax rates must be in [0, 1]
#
# If any are outside [0, 1], the formula is wrong or the inputs are mis-scaled.

print("="*70)
print("CHECK 1: Tax rates should be in [0, 1]")
print("="*70)
tax_vars = ["tau_l", "tau_k", "tau_c", "tau_d", "tau_itc", "e_tau",
            "tau_i", "tau_div"]
for v in tax_vars:
    if v not in df.columns: continue
    s = df[v].dropna()
    out_of_range = ((s < 0) | (s > 1)).sum()
    flag = "  <-- BROKEN" if out_of_range > 0 else "  OK"
    print(f"  {v:10s}  min={s.min():>8.4f}  max={s.max():>8.4f}  "
          f"mean={s.mean():>7.4f}  out_of_[0,1]={out_of_range}{flag}")


#%%
# =============================================================================
# CHECK 2: Sample means over Zubairy window
#
# Zubairy (2014), p.178 ("Calibration and Priors"):
#   - g/y = 0.18                       (govt share of GDP)
#   - b/y = 0.33                       (annual debt/GDP) -- see CHECK 8
#   - tau_k mean ~ 0.41                (capital tax rate)
#   - tau_l mean ~ 0.23                (labor tax rate)

print()
print("="*70)
print("CHECK 2: Sample means over Zubairy window (1958Q1-2008Q4)")
print("="*70)
print(f"  {'series':12s}  {'ours':>10s}  {'Zubairy':>10s}  {'diff':>8s}")
print(f"  {'-'*12}  {'-'*10}  {'-'*10}  {'-'*8}")
targets = {"tau_l": 0.23, "tau_k": 0.41}
for v, target in targets.items():
    if v not in zub.columns: continue
    ours = zub[v].mean()
    print(f"  {v:12s}  {ours:>10.4f}  {target:>10.4f}  {ours-target:>+8.4f}")


#%%
# =============================================================================
# CHECK 3: Stationarity (Augmented Dickey-Fuller)
#
# Hat versions should reject the unit-root null (small p-value < 0.05).

print()
print("="*70)
print("CHECK 3: ADF test on hat series (p < 0.05 means stationary)")
print("="*70)
hat_vars = ["chat", "ihat", "ghat", "bhat", "lhat", "trhat",
            "tau_l_hat_detrend", "tau_l_hat_nodetrend",
            "tau_k_hat_detrend", "tau_k_hat_nodetrend",
            "tau_c_hat_detrend", "tau_c_hat_nodetrend",
            "tau_d_hat_detrend", "tau_d_hat_nodetrend",
            "tau_i_hat_detrend", "tau_i_hat_nodetrend",
            "tau_div_hat_detrend", "tau_div_hat_nodetrend",
            "tau_itc_hat", "e_tau_hat", "pi", "r"]
for v in hat_vars:
    if v not in zub.columns: continue
    s = zub[v].dropna().values
    if len(s) < 10:
        continue
    try:
        adf, pval = adfuller(s, autolag="AIC")[:2]
        flag = "stationary" if pval < 0.05 else "UNIT ROOT (concerning)"
        print(f"  {v:22s}  ADF={adf:>7.3f}  p={pval:>6.4f}  {flag}")
    except Exception as e:
        print(f"  {v:22s}  ADF failed: {e}")


#%%
# =============================================================================
# CHECK 4: Plot real-quantity hat series + inflation
#
# Should look stationary around 0 (or 1, for pi).

fig, axes = plt.subplots(7, 2, figsize=(13, 18))

panels = [
    ("Consumption", None, "chat"),
    ("Investment",  None, "ihat"),
    ("Government",  None, "ghat"),
    ("Debt",        None, "bhat"),
    ("Hours",       None, "lhat"),
    ("Transfers",   None, "trhat"),
    ("Inflation",   "pi", None),
]
for i, (label, lvl, hat) in enumerate(panels):
    ax_l = axes[i, 0]
    ax_r = axes[i, 1]
    series_to_plot = hat or lvl
    df[series_to_plot].plot(ax=ax_l, color="navy", lw=0.8)
    ax_l.set_title(f"{label} ({series_to_plot})")
    if hat:
        ax_l.axhline(0, color="k", lw=0.5)
    df[series_to_plot].diff().plot(ax=ax_r, color="darkred", lw=0.6)
    ax_r.set_title("first difference (residual diagnostic)")
    ax_r.axhline(0, color="k", lw=0.5)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "diagnostic_real_series.png"), dpi=110)
plt.show()


#%%
# =============================================================================
# CHECK 5: Plot tax rate levels (untransformed) with linear-trend overlay
#
# Visible upward/downward slopes are evidence FOR linear detrending; weak
# slopes or break-driven movement are evidence for demean-only. Use this
# alongside CHECK 6 and CHECK 6b to inform HAT_CHOICE in 3_clean_for_matlab.py.

fig, axes = plt.subplots(4, 2, figsize=(13, 12))

tax_panels = [
    ("Labor tax (tau_l)",          "tau_l"),
    ("Capital tax (tau_k)",        "tau_k"),
    ("Consumption tax (tau_c)",    "tau_c"),
    ("Dividend tax proxy (tau_d)", "tau_d"),
    ("Interest tax (tau_i)",       "tau_i"),
    ("Dividend tax TAXSIM (tau_div)", "tau_div"),
    ("ITC (tau_itc)",              "tau_itc"),
    ("Bonus expensing (e_tau)",    "e_tau"),
]
for i, (label, lvl) in enumerate(tax_panels):
    ax_l = axes[i // 2, i % 2]
    df[lvl].plot(ax=ax_l, color="navy", lw=1.0, label="level")
    s = df[lvl].dropna()
    if len(s) < 2:
        ax_l.set_title(label + "  (no data)")
        continue
    t = np.arange(len(s))
    coef = np.polyfit(t, s.values, 1)
    trend = pd.Series(np.polyval(coef, t), index=s.index)
    trend.plot(ax=ax_l, color="orange", ls="--", lw=1.0,
               label=f"linear trend (slope={coef[0]:.5f})")
    ax_l.axhline(s.mean(), color="gray", ls=":", lw=0.7, label="mean")
    ax_l.set_title(label)
    ax_l.legend(loc="best", fontsize=8)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "diagnostic_tax_levels.png"), dpi=110)
plt.show()


#%%
# =============================================================================
# CHECK 6: Slope-significance test
#
# Formal version of CHECK 5. Tests whether the linear-trend slope differs
# from zero. NOTE: with N > 200, statistical significance becomes automatic
# even for tiny or break-driven movements -- always read this alongside
# CHECK 5 (visual) and CHECK 6b (economic magnitude). Do not use this
# check alone to decide HAT_CHOICE.

print()
print("="*70)
print("CHECK 6: Linear-trend slope significance for each tax rate")
print("        (caveat: significance does not imply economically meaningful drift)")
print("="*70)
print(f"  {'series':10s}  {'slope':>11s}  {'t-stat':>8s}  {'p-val':>8s}  flag")
print(f"  {'-'*10}  {'-'*11}  {'-'*8}  {'-'*8}  {'-'*30}")
for v in ["tau_l", "tau_k", "tau_c", "tau_d", "tau_itc", "e_tau",
          "tau_i", "tau_div"]:
    if v not in zub.columns: continue
    s = zub[v].dropna().values
    if len(s) < 10:
        continue
    t = np.arange(len(s))
    res = stats.linregress(t, s)
    flag = "slope significantly != 0" if res.pvalue < 0.05 else "slope not distinguishable from 0"
    print(f"  {v:10s}  {res.slope:>+11.6f}  {res.slope/res.stderr:>+8.2f}  "
          f"{res.pvalue:>8.4f}  {flag}")


#%%
# =============================================================================
# CHECK 6b: Economic magnitude of trend (not just statistical significance)
#
# CHECK 6 flags any non-zero slope as significant once N is large, and can't
# tell a true trend from a structural break (ITC repeal, bonus episodes).
# This check reports three magnitude metrics:
#
#   - drift_rel    = (slope * N) / mean         (total drift relative to level)
#   - R^2          = fraction of variance the trend explains
#   - break_flag   = YES if residual variance differs by 4x between halves
#                    (suggests regime change rather than smooth drift)

print()
print("="*70)
print("CHECK 6b: Economic magnitude of trend (not just stat. significance)")
print("="*70)
print(f"  {'series':10s}  {'drift_rel':>10s}  {'R^2':>6s}  {'break?':>7s}  evidence pattern")
print(f"  {'-'*10}  {'-'*10}  {'-'*6}  {'-'*7}  {'-'*40}")

for v in ["tau_l", "tau_k", "tau_c", "tau_d", "tau_itc", "e_tau",
          "tau_i", "tau_div"]:
    if v not in zub.columns: continue
    s = zub[v].dropna().values
    n = len(s)
    if n < 10:
        continue
    t = np.arange(n)
    res = stats.linregress(t, s)
    fit = res.intercept + res.slope * t
    resid = s - fit

    # Drift relative to mean (absolute value)
    drift_rel = abs(res.slope * n / s.mean()) if s.mean() != 0 else float('inf')

    # R^2 of linear trend
    ss_tot = ((s - s.mean())**2).sum()
    ss_res = (resid**2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    # Break flag: residuals concentrated in one half of sample suggest
    # structural break rather than smooth drift.
    half = n // 2
    var_ratio = resid[:half].var() / resid[half:].var() if resid[half:].var() > 0 else float('inf')
    is_break = (var_ratio > 4) or (var_ratio < 0.25)

    # Categorize evidence pattern (descriptive, not prescriptive)
    if is_break:
        pattern = "structural-break pattern"
    elif drift_rel >= 0.20 and r2 >= 0.30:
        pattern = "strong, smooth trend"
    elif drift_rel >= 0.10 and r2 >= 0.15:
        pattern = "moderate trend"
    else:
        pattern = "weak / mostly cyclical"

    print(f"  {v:10s}  {drift_rel:>10.3f}  {r2:>6.3f}  {'YES' if is_break else 'no':>7s}  {pattern}")

print("""
  Reading guide:
    drift_rel  = total drift over sample as fraction of mean.
                   Above ~0.20 indicates economically meaningful drift.
    R^2        = how well a linear trend fits the data.
                   Below ~0.20 indicates movement is mostly cyclical or
                   break-driven, not smooth trend.
    break?     = YES if residual variance differs by 4x between halves,
                   suggesting a regime change rather than a smooth trend.
                   In that case a linear detrend is the wrong fix; level
                   deviation captures the structure better.

  Interpreting patterns:
    "strong, smooth trend"     -- detrending is well-supported
    "moderate trend"           -- detrending defensible, also fine to demean only
    "weak / mostly cyclical"   -- demean-only preserves more signal
    "structural-break pattern" -- linear detrend would mis-fit; consider
                                  level deviation or break dummies

  These are diagnostic categories, not prescriptions. Combine with CHECK 5
  (visual) to decide each entry of HAT_CHOICE in 3_clean_for_matlab.py.
""")

print()
print("="*70)
print("CHECK 7: Diagnostic context for tax-variable transformation choices")
print("="*70)
print("""
  This script DIAGNOSES whether each tax rate needs detrending; it does NOT
  decide. The decision is made in 3_clean_for_matlab.py via HAT_CHOICE.

  1_build_data.py produces TWO hat variants for each detrendable tax rate:
    tau_l_hat_detrend   tau_l_hat_nodetrend
    tau_k_hat_detrend   tau_k_hat_nodetrend
    tau_c_hat_detrend   tau_c_hat_nodetrend
    tau_d_hat_detrend   tau_d_hat_nodetrend

  Both are produced symmetrically; neither is privileged at build time.
  3_clean_for_matlab.py reads HAT_CHOICE and renames the chosen variant to
  the canonical tau_x_hat for the .mod file.

  tau_itc and e_tau use level deviation (tau - mean) without a trend
  because both are exactly zero for ~75% of the sample; (tau-mean)/mean
  would explode. Their step-function shapes are structural breaks
  (ITC repeal 1986; bonus depreciation episodes from 2001 onward), not
  smooth trends.

  Use the diagnostic evidence above to inform HAT_CHOICE:
    CHECK 5 -- visual inspection of each tax level with linear-trend overlay
    CHECK 6 -- statistical significance of the linear-trend slope
    CHECK 6b -- economic magnitude (drift_rel, R^2, structural-break flag)

  When CHECK 6 and CHECK 6b disagree (e.g. CHECK 6 says "DETREND" but
  CHECK 6b flags the residuals as a structural break), trust CHECK 6b plus
  the level plot. Slope significance becomes automatic with N > 200 even
  for tiny or break-driven movements.

  Useful reference points from the literature:
    Jones (2002): detrend tau_l (rising payroll/SS contributions are not
                  in standard DSGE models); do not detrend tau_k.
    Mendoza-Razin-Tesar (1994): demean-only tau_c.
""")


#%%
# =============================================================================
# CHECK 8: Debt-to-GDP convention check
#
# Verifies that average-debt-as-stock + SAAR-GDP-as-flow gives the right units
# (years of annual GDP), matching Zubairy's calibration target b_bar/y_bar = 0.33.
#
# Pulls FRED GDP (SAAR) and the debt series from the existing FRED MVPHGFD,
# averaging the monthly debt to a quarterly stock.

def _fred_get(series_id, start="1947-01-01"):
    url = "https://api.stlouisfed.org/fred/series/observations"
    params = {"series_id": series_id, "api_key": FRED_KEY, "file_type": "json",
              "observation_start": start}
    for attempt in range(4):
        try:
            r = requests.get(url, params=params, timeout=120)
            r.raise_for_status()
            obs = r.json().get("observations", [])
            df_ = pd.DataFrame(obs)
            df_["date"]  = pd.to_datetime(df_["date"])
            df_["value"] = pd.to_numeric(df_["value"], errors="coerce")
            return df_.set_index("date")["value"].dropna()
        except (requests.ConnectionError, requests.Timeout):
            if attempt < 3:
                wait = 2.0 ** attempt
                print(f"    [retry {attempt+1}/3] timeout; waiting {wait:.1f}s")
                time.sleep(wait)
            else:
                raise

print("="*70)
print("CHECK 8: Debt-to-GDP convention check")
print("="*70)
print("Fetching GDP and debt from FRED for the ratio...")
gdp_m = _fred_get("GDP")
debt_m = _fred_get("MVPHGFD027MNFRBDAL")

gdp_q = gdp_m.copy(); gdp_q.index = gdp_q.index.to_period("Q")
debt_q = debt_m.groupby(debt_m.index.to_period("Q")).mean()

debt_to_gdp = (debt_q / gdp_q).dropna()
zub_window = debt_to_gdp.loc[pd.Period("1958Q1","Q"):pd.Period("2008Q4","Q")]
full_window = debt_to_gdp.loc[pd.Period(START_Q,"Q"):pd.Period(END_Q,"Q")]
print(f"\nDebt/GDP sanity check:")
print(f"  1958Q1-2008Q4 (Zubairy sample) mean: {zub_window.mean():.3f}  (target ~0.33)")
print(f"  {START_Q}-{END_Q}    (full sample)   mean: {full_window.mean():.3f}")
print(f"\n  ratio has units of 'years of annual GDP' because:")
print(f"    numerator: average debt outstanding during the quarter (stock, $bn)")
print(f"    denominator: GDP at annualized rate (flow, $bn)")
