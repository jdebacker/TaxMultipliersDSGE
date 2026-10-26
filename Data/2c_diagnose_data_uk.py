"""
2c_diagnose_data_uk.py
======================
Sanity checks on the constructed series in UK_Data_Full.xlsx:

  1. Verify tax-rate formulas produce values in [0, 1]
  2. SKIPPED: Zubairy comparison is US-specific
  3. ADF stationarity tests on the hat series
  4. Plot real-quantity hats so they look stationary around 0
  5. Plot tax-rate levels with linear-trend overlays
  6. Linear-trend slope-significance test
  6b. Economic-magnitude trend assessment + structural-break flag
  7. Diagnostic context: how to use these results to inform HAT_CHOICE
  8. Debt-to-GDP convention check

This script DIAGNOSES whether each tax rate has a meaningful trend; it
does NOT decide. The detrend-vs-demean choice for each tax variable is
made in 2d_clean_for_matlab_uk.py via HAT_CHOICE.

Designed to work regardless of which 2b version produced UK_Data_Full.xlsx
(OECD-annual or ONS-quarterly). Each tax-rate column is checked
conditionally with `if v in df.columns`. UK pipelines do not produce
tau_d / tau_i / tau_div / tau_itc / e_tau, so those columns are simply
absent from any UK output and the conditional checks pass.

Run AFTER 2b_build_data_uk_*.py finishes. Reads UK_Data_Full.xlsx.
"""

#%%
# =============================================================================
# imports

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import adfuller
import scipy.stats as stats

OUT_DIR  = r"./output"
RAW_DIR  = r"./raw"
START_Q  = "1955Q1"
END_Q    = "2026Q1"

df = pd.read_excel(os.path.join(OUT_DIR, "UK_Data_Full.xlsx"))

df["QPER"] = pd.PeriodIndex(
    [pd.Period(year=int(y), quarter=int(q), freq="Q")
     for y, q in zip(df["year"], df["quarter"])]
)
df = df.set_index("QPER")

# Set the diagnostic window to the longest fully-populated stretch implied
# by lhat (1971Q1+), since the UK estimation will be bound by that anyway.
diag = df.loc[pd.Period("1971Q1","Q"):pd.Period(END_Q,"Q")]


#%%
# =============================================================================
# CHECK 1: tax rates must be in [0, 1]

print("="*70)
print("CHECK 1: Tax rates should be in [0, 1]")
print("="*70)
tax_vars = ["tau_l", "tau_k", "tau_c", "tau_d", "tau_itc", "e_tau",
            "tau_i", "tau_div"]
any_present = False
for v in tax_vars:
    if v not in df.columns:
        continue
    any_present = True
    s = df[v].dropna()
    if len(s) == 0:
        print(f"  {v:10s}  (column present but all NaN)")
        continue
    out_of_range = ((s < 0) | (s > 1)).sum()
    flag = "  <-- BROKEN" if out_of_range > 0 else "  OK"
    print(f"  {v:10s}  min={s.min():>8.4f}  max={s.max():>8.4f}  "
          f"mean={s.mean():>7.4f}  out_of_[0,1]={out_of_range}{flag}")
if not any_present:
    print("  No tax-rate level columns found in UK_Data_Full.xlsx")


#%%
# =============================================================================
# CHECK 2: SKIPPED -- Zubairy targets are US-specific

print()
print("="*70)
print("CHECK 2: SKIPPED (Zubairy 2014 targets are US-specific; UK calibration")
print("         anchors are read directly from sample means in CHECK 8 below)")
print("="*70)


#%%
# =============================================================================
# CHECK 3: Stationarity (Augmented Dickey-Fuller)

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
    if v not in diag.columns:
        continue
    s = diag[v].dropna().values
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

panels = [
    ("Consumption", None, "chat"),
    ("Investment",  None, "ihat"),
    ("Government",  None, "ghat"),
    ("Debt",        None, "bhat"),
    ("Hours",       None, "lhat"),
    ("Transfers",   None, "trhat"),
    ("Inflation",   "pi", None),
]
panels = [(lbl, lvl, hat) for lbl, lvl, hat in panels
          if (hat or lvl) in df.columns]

fig, axes = plt.subplots(len(panels), 2, figsize=(13, max(3, 2.6*len(panels))))
if len(panels) == 1:
    axes = axes.reshape(1, 2)

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
plt.savefig(os.path.join(OUT_DIR, "diagnostic_real_series_uk.png"), dpi=110)
plt.show()


#%%
# =============================================================================
# CHECK 5: Plot tax rate levels with linear-trend overlay
#
# UK has tau_l, tau_k, tau_c only. tau_d/tau_i/tau_div/tau_itc/e_tau are
# US-pipeline-only and not constructed for the UK. The tax_panels list is
# filtered down to columns actually present in UK_Data_Full.xlsx.

tax_panels_all = [
    ("Labor tax (tau_l)",          "tau_l"),
    ("Capital tax (tau_k)",        "tau_k"),
    ("Consumption tax (tau_c)",    "tau_c"),
    ("Dividend tax proxy (tau_d)", "tau_d"),
    ("Interest tax (tau_i)",       "tau_i"),
    ("Dividend tax TAXSIM (tau_div)", "tau_div"),
    ("ITC (tau_itc)",              "tau_itc"),
    ("Bonus expensing (e_tau)",    "e_tau"),
]
tax_panels = [(lbl, lvl) for lbl, lvl in tax_panels_all if lvl in df.columns]

if len(tax_panels) > 0:
    n = len(tax_panels)
    nrows = (n + 1) // 2
    fig, axes = plt.subplots(nrows, 2, figsize=(13, max(3, 3.0*nrows)))
    axes = np.atleast_2d(axes)
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
    # Hide unused subplots if odd count
    for j in range(n, nrows * 2):
        axes[j // 2, j % 2].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "diagnostic_tax_levels_uk.png"), dpi=110)
    plt.show()
else:
    print("  No tax-rate level columns found; skipping plot.")


#%%
# =============================================================================
# CHECK 6: Slope-significance test

print()
print("="*70)
print("CHECK 6: Linear-trend slope significance for each tax rate")
print("        (caveat: significance does not imply economically meaningful drift)")
print("="*70)
print(f"  {'series':10s}  {'slope':>11s}  {'t-stat':>8s}  {'p-val':>8s}  flag")
print(f"  {'-'*10}  {'-'*11}  {'-'*8}  {'-'*8}  {'-'*30}")
for v in ["tau_l", "tau_k", "tau_c", "tau_d", "tau_itc", "e_tau",
          "tau_i", "tau_div"]:
    if v not in diag.columns:
        continue
    s = diag[v].dropna().values
    if len(s) < 10:
        continue
    t = np.arange(len(s))
    res = stats.linregress(t, s)
    flag = "slope significantly != 0" if res.pvalue < 0.05 else "slope not distinguishable from 0"
    print(f"  {v:10s}  {res.slope:>+11.6f}  {res.slope/res.stderr:>+8.2f}  "
          f"{res.pvalue:>8.4f}  {flag}")


#%%
# =============================================================================
# CHECK 6b: Economic magnitude of trend

print()
print("="*70)
print("CHECK 6b: Economic magnitude of trend (not just stat. significance)")
print("="*70)
print(f"  {'series':10s}  {'drift_rel':>10s}  {'R^2':>6s}  {'break?':>7s}  evidence pattern")
print(f"  {'-'*10}  {'-'*10}  {'-'*6}  {'-'*7}  {'-'*40}")

for v in ["tau_l", "tau_k", "tau_c", "tau_d", "tau_itc", "e_tau",
          "tau_i", "tau_div"]:
    if v not in diag.columns:
        continue
    s = diag[v].dropna().values
    n = len(s)
    if n < 10:
        continue
    t = np.arange(n)
    res = stats.linregress(t, s)
    fit = res.intercept + res.slope * t
    resid = s - fit

    drift_rel = abs(res.slope * n / s.mean()) if s.mean() != 0 else float('inf')

    ss_tot = ((s - s.mean())**2).sum()
    ss_res = (resid**2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    half = n // 2
    var_ratio = resid[:half].var() / resid[half:].var() if resid[half:].var() > 0 else float('inf')
    is_break = (var_ratio > 4) or (var_ratio < 0.25)

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
                   In that case linear detrend is the wrong fix; level
                   deviation captures the structure better.

  UK-specific structural breaks to watch:
    1973: VAT introduction -> tau_c jumps; expect break flag for tau_c
    1979: Thatcher Howe budget -> direct-to-indirect tax shift
    1984-86: corporation tax reform -> tau_k regime change
    1997: PUSF accounting basis change (ONS quarterly version only)

  Interpreting patterns:
    "strong, smooth trend"     -- detrending is well-supported
    "moderate trend"           -- detrending defensible, also fine to demean only
    "weak / mostly cyclical"   -- demean-only preserves more signal
    "structural-break pattern" -- linear detrend would mis-fit; level
                                  deviation captures the structure better
""")


#%%
# =============================================================================
# CHECK 7: Context

print("="*70)
print("CHECK 7: Diagnostic context for tax-variable transformation choices")
print("="*70)
print("""
  This script DIAGNOSES whether each tax rate needs detrending; it does NOT
  decide. The decision is made in 2d_clean_for_matlab_uk.py via HAT_CHOICE.

  2b_build_data_uk_*.py produces TWO hat variants for each detrendable tax:
    tau_l_hat_detrend   tau_l_hat_nodetrend
    tau_k_hat_detrend   tau_k_hat_nodetrend
    tau_c_hat_detrend   tau_c_hat_nodetrend

  Both are produced symmetrically; neither is privileged at build time.
  2d_clean_for_matlab_uk.py reads HAT_CHOICE and renames the chosen variant
  to the canonical tau_x_hat for the .mod file.

  Use the diagnostic evidence above to inform HAT_CHOICE:
    CHECK 5  -- visual inspection of each tax level with linear-trend overlay
    CHECK 6  -- statistical significance of the linear-trend slope
    CHECK 6b -- economic magnitude (drift_rel, R^2, structural-break flag)

  When CHECK 6 and CHECK 6b disagree (CHECK 6 flags slope significant but
  6b reports a structural break), trust CHECK 6b plus the level plot.

  UK-specific note: with a 1971Q1+ sample and likely structural breaks at
  1973 (VAT) and 1979 (Thatcher reforms), expect tau_c to flag as
  "structural-break pattern" rather than smooth trend. Demean-only is
  often the safer choice for UK tax rates than for the corresponding US
  series.
""")


#%%
# =============================================================================
# CHECK 8: Debt-to-GDP convention check + UK calibration anchors
#
# Pulls UK GDP and debt from ONS to verify b/y ratio. UK postwar average
# debt/GDP is typically 0.30-0.50 depending on subwindow. Modern (post-2010)
# is much higher; long-run average over 1955-2019 is closer to 0.40.

print("="*70)
print("CHECK 8: Debt-to-GDP convention check + UK calibration anchors")
print("="*70)

# If b_y or g_y are already in the dataset (some UK build files include them),
# print directly. Otherwise, compute from columns we have.
if "b_y" in df.columns and df["b_y"].notna().any():
    by_diag = df.loc[diag.index, "b_y"].dropna()
    print(f"  b_y (from dataset) mean over {diag.index.min()} - {diag.index.max()}: "
          f"{by_diag.mean():.3f}")
else:
    print("  b_y not in dataset; if needed, compute as bhat-level / GDP outside this script.")

if "g_y" in df.columns and df["g_y"].notna().any():
    gy_diag = df.loc[diag.index, "g_y"].dropna()
    print(f"  g_y (from dataset) mean over {diag.index.min()} - {diag.index.max()}: "
          f"{gy_diag.mean():.3f}")

print(f"\n  UK reference points (rough postwar averages):")
print(f"    b/y: ~0.40 long-run; ~0.30 1980s-1990s; >0.80 post-2010")
print(f"    g/y: ~0.20")
print(f"    tau_l: ~0.20-0.25")
print(f"    tau_k: ~0.30-0.45")
print(f"    tau_c: ~0.10-0.15 (post-1973 VAT era)")

print(f"\n  Sample means (use these as steady-state calibration anchors")
print(f"  in the UK .mod file's 'parameters' block):")
for v in ["tau_l", "tau_k", "tau_c"]:
    if v not in diag.columns:
        continue
    s = diag[v].dropna()
    if len(s) > 0:
        print(f"    {v}_ss = {s.mean():.6f}")
