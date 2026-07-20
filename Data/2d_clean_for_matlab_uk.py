"""
2d_clean_for_matlab_uk.py
=========================
Prepare UK MATLAB-ready datasets for DSGE estimation.

Input:
    output/UK_Data_Full.xlsx

Outputs:
    output/UK_Data_Matlab_tax.xlsx    + .mat   1997Q2 - 2016Q4
        Balanced over the tax window. tau_l, tau_k, tau_c all present with no
        NaN (tau_k is the binding constraint, ending 2016Q4). Use this for the
        estimation that includes tax-rate observables in varobs.

    output/UK_Data_Matlab_macro.xlsx  + .mat   1971Q1 - 2019Q4
        Long macro sample. chat, ihat, ghat, lhat, pi_obs, r_obs run from
        1971Q1; bhat and trhat have leading NaN (debt ~1990+, transfers ~1988+).
        Tax-rate columns are dropped from this file (they only start 1997Q2 and
        are not used as observables in the long-sample spec). Use this for the
        macro-only estimation; choose in the .mod's varobs whether to include
        bhat/trhat given their later start.

For tau_l, tau_k, tau_c, 2b_build_data_uk.py produces both detrended and
nodetrend hat variants. HAT_CHOICE selects which becomes the canonical
tau_x_hat in the tax file.
"""

#%%
# =============================================================================
# imports and config

import os
import pandas as pd
from scipy.io import savemat

OUT_DIR = r"./output"
SRC = os.path.join(OUT_DIR, "UK_Data_Full.xlsx")

WINDOW_TAX   = (pd.Period("1997Q2", "Q"), pd.Period("2016Q4", "Q"))
WINDOW_MACRO = (pd.Period("1971Q1", "Q"), pd.Period("2019Q4", "Q"))

# Pick which variant becomes the canonical tau_x_hat for the .mod file.
# Set each value to "detrend" or "nodetrend".
HAT_CHOICE = {
    "tau_l": "detrend",     # R^2=0.49, clean upward trend
    "tau_k": "nodetrend",   # R^2=0.15, mostly cyclical, linear fit poor
    "tau_c": "detrend",     # R^2=0.36, clean upward trend; 2009 dip is a level event
}

# Macro-window observables kept in UK_Data_Matlab_macro.xlsx. bhat and trhat
# are retained even though they have leading NaN (debt ~1990+, transfers ~1988+);
# the .mod varobs decides whether to use them.
MACRO_COLS = ["chat", "ihat", "ghat", "bhat", "pi_obs", "r_obs", "lhat", "trhat"]


#%%
# =============================================================================
# load full dataset

df = pd.read_excel(SRC)
df["QPER"] = pd.PeriodIndex(
    [pd.Period(year=int(y), quarter=int(q), freq="Q")
     for y, q in zip(df["year"], df["quarter"])]
)
df = df.drop(columns=["year", "quarter"]).set_index("QPER").sort_index()

print(f"Loaded {SRC}")
print(f"  rows: {len(df)}  cols: {df.shape[1]}")
print(f"  range: {df.index.min()} -> {df.index.max()}")
print("\nHAT_CHOICE:")
for tax, choice in HAT_CHOICE.items():
    print(f"  {tax}_hat <- {tax}_hat_{choice}")


#%%
# =============================================================================
# helper: write an xlsx + mat pair from a windowed, renamed frame

def export_window(frame, window, out_basename, label):
    a, b = window
    sub = frame.loc[a:b].copy()
    sub.insert(0, "quarter", sub.index.quarter)
    sub.insert(0, "year", sub.index.year)
    sub = sub.reset_index(drop=True)

    out_xlsx = os.path.join(OUT_DIR, out_basename + ".xlsx")
    sub.to_excel(out_xlsx, index=False)

    out_mat = os.path.join(OUT_DIR, out_basename + ".mat")
    savemat(out_mat,
            {col: sub[col].to_numpy().reshape(-1, 1) for col in sub.columns},
            do_compression=True)

    print(f"\n[{label}] wrote:")
    print(f"  {out_xlsx}")
    print(f"  {out_mat}")
    print(f"  window: {a} -> {b}   rows: {len(sub)}")
    print(f"  columns: {sub.columns.tolist()}")
    print(f"  non-NaN by series:")
    for col in sub.columns:
        if col in ("year", "quarter"):
            continue
        n = sub[col].notna().sum()
        first = sub.loc[sub[col].notna(), col].index.min() if n else None
        first_yq = (f"{sub.loc[first,'year']}Q{sub.loc[first,'quarter']}"
                    if first is not None else "—")
        tag = "" if n == len(sub) else f"   (starts {first_yq}; {len(sub)-n} leading NaN)"
        print(f"    {col:12s}  {n:>4d} obs{tag}")
    return sub


#%%
# =============================================================================
# TAX FILE: 1997Q2 - 2016Q4, balanced, with canonical tau_x_hat

# Drop level series, detrended-level series, and the unused hat variant.
drop_cols = []
for tax, choice in HAT_CHOICE.items():
    drop_cols += [tax, f"{tax}_detrend"]
    other = "nodetrend" if choice == "detrend" else "detrend"
    drop_cols.append(f"{tax}_hat_{other}")

rename_map = {f"{tax}_hat_{choice}": f"{tax}_hat"
              for tax, choice in HAT_CHOICE.items()}
rename_map.update({"pi": "pi_obs", "r": "r_obs"})

df_tax = df.drop(columns=[c for c in drop_cols if c in df.columns])
df_tax = df_tax.rename(columns=rename_map)

export_window(df_tax, WINDOW_TAX, "UK_Data_Matlab_tax", "TAX 1997Q2-2016Q4")


#%%
# =============================================================================
# MACRO FILE: 1971Q1 - 2019Q4, macro observables only (tax columns dropped)

# Keep only the macro observables; rename pi/r. Drop every tax-related column.
df_macro = df.rename(columns={"pi": "pi_obs", "r": "r_obs"})
df_macro = df_macro[[c for c in MACRO_COLS if c in df_macro.columns]]

export_window(df_macro, WINDOW_MACRO, "UK_Data_Matlab_macro", "MACRO 1971Q1-2019Q4")
