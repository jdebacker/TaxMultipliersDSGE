# Series construction summary

This document describes how every observable is constructed from raw FRED, BEA NIPA, and BEA Fixed Assets inputs. The pipeline runs in four scripts:

| Script | Role | Output |
|---|---|---|
| `1a_pre_prep.py` | Pull BEA Fixed Assets Table 2.7, build investment weights and bonus-eligibility flags by detailed asset type | `raw/PQdat_BEA.csv`, `raw/treat_BEA.csv`, `raw/asset_types.txt` |
| `1b_build_data.py` | Pull raw data from APIs and construct all observables | `US_Data_Full.xlsx` |
| `1c_diagnose_data.py` | Sanity checks (range, stationarity, trend tests, plots, debt/GDP ratio) | console + 2 PNGs |
| `1d_clean_for_matlab.py` | Pick hat variants and slice windows for estimation | `US_Data_Matlab.xlsx`, `US_Data_Matlab_Zubairy.xlsx` |

---

## 1. Overview

The dataset contains observables in three groups:

| Group | Series | Source | Role in model |
|---|---|---|---|
| **A. Macro baseline** | `chat`, `ihat`, `ghat`, `bhat`, `pi`, `r`, `tau_l_hat`, `tau_k_hat` | FRED + BEA NIPA | Standard NK DSGE observables |
| **B. Extensions** | `lhat`, `trhat`, `tau_c_hat`, `tau_d_hat` | FRED + BEA NIPA | Hours, transfers, consumption tax, dividend tax |
| **C. Tax policy** | `tau_itc_hat`, `e_tau_hat` | BEA Fixed Assets Table 2.7 + statutory bonus calendar | Investment tax credit, bonus depreciation expensing share |

For each detrendable tax rate (`tau_l`, `tau_k`, `tau_c`, `tau_d`), `1b_build_data.py` produces **both** `tau_x_hat_detrend` and `tau_x_hat_nodetrend` — symmetrically named, neither privileged in the build. The selection of which variant becomes the canonical `tau_x_hat` is made in `1d_clean_for_matlab.py` via the `HAT_CHOICE` dict, so robustness comparisons are a one-line edit. See Section 3.

Default sample window in `1b_build_data.py`: **1958Q1 – 2026Q1**. The MATLAB-export script trims to two windows for estimation:

- `US_Data_Matlab.xlsx`: 1958Q1 – 2019Q4 (extended sample)
- `US_Data_Matlab_Zubairy.xlsx`: 1958Q1 – 2008Q4 (Zubairy's sample window, useful for direct comparison)

All inputs are pulled programmatically:
- FRED via the `fred/series/observations` endpoint
- BEA NIPA via the `apps.bea.gov/api` endpoint (NIPA dataset)
- BEA Fixed Assets Table 2.7 via the `apps.bea.gov/api` endpoint (FixedAssets dataset; pulled by `1a_pre_prep.py`) — produces `PQdat_BEA.csv` and `treat_BEA.csv` in `raw/`
- API keys read from `raw/api_keys.txt` (`FRED_KEY=...` and `BEA_KEY=...`, one per line)

The script handles two BEA-specific quirks automatically:

1. **NIPA line-number drift.** BEA renumbers its tables periodically. `find_line()` discovers each line by description substring and prints the line number it found. Sanity-check ranges flag any mismatches.

2. **Unit scaling.** BEA's API returns NIPA flows in *millions* of dollars (`UNIT_MULT=6`), while FRED returns most aggregates in *billions*. The `bea_nipa()` wrapper reads `UNIT_MULT` from the response and rescales to billions, so all downstream arithmetic is in consistent units.

---

## 2. Aggregation conventions

The dataset combines three semantically different types of variables, each with its own appropriate aggregation rule:

| Type | Examples | Convention | Why |
|---|---|---|---|
| **Flow** (accumulates over time) | GDP, C, I, G, transfers, tax revenues | BEA SAAR (seasonally adjusted annual rate, = quarterly flow × 4) | Allows direct comparison to annual aggregates. All NIPA flows share this convention, so the GDP identity `Y = C + I + G + NX` holds exactly in the published data. |
| **Stock** (level at a point in time) | Federal debt, population | Quarterly mean of monthly snapshots | A stock has no "amount per quarter" — summing monthly values produces the meaningless unit "dollar-months." The mean represents "average level outstanding during the quarter." |
| **Rate / index** | Federal funds rate, GDP deflator, hours index | Quarterly mean | Rates have no quantity to accumulate. |

This mixing is **intentional and economically correct**, not an inconsistency. The ratio `b_t / y_t` (debt-to-GDP) pairs an average-stock numerator with a SAAR-flow denominator and yields units of *years of annual output*. CHECK 8 in `1c_diagnose_data.py` verifies this ratio sits in the expected ~0.30 range over the 1958–2008 window.

The three monthly→quarterly aggregations done in-script all use **mean**:
- `FEDFUNDS` (rate)
- `CNP16OV` (population stock)
- `MVPHGFD027MNFRBDAL` (debt stock)

---

## 3. Detrending conventions

The DSGE has a stationary steady state. Real per-capita aggregates trend over time due to forces outside the model (TFP growth, capital deepening, demographic shifts), so they're log-linearly detrended before estimation. Tax *rates* are bounded and can't grow unboundedly, but several show slow drift driven by mechanisms the model doesn't capture (rising payroll/social-insurance contributions, services-share shift, etc.) — so detrending may make sense for those too.

| Treatment | Applied to | Rationale |
|---|---|---|
| `100 × (log x − linear trend)` | `chat`, `ihat`, `ghat`, `bhat`, `lhat`, `trhat` | Real per-capita levels grow due to forces not in the model. The detrended residual is the cyclical component the model speaks to. |
| Both fractional-deviation variants produced (`_hat_detrend` AND `_hat_nodetrend`) | `tau_l`, `tau_k`, `tau_c`, `tau_d` | Each tax rate gets two columns; the choice between them is made downstream via `HAT_CHOICE` in `1d_clean_for_matlab.py`. See 3.1. |
| Demean as level: `τ − μ` | `tau_itc_hat`, `e_tau_hat` | Both are policy series exactly zero for ~75% of the sample. Fractional deviation `(τ−μ)/μ` would explode; their step-function shapes are structural breaks (ITC repeal 1986; bonus depreciation episodes from 2001), not smooth trends. |
| None (kept as level) | `pi`, `r` | Inflation and the federal funds rate are bounded and naturally stationary. |

### 3.1 Symmetric hat variants for the tax rates

`1b_build_data.py` produces both variants for every detrendable tax rate, with no asymmetry in the build output:

| Series | Build output columns |
|---|---|
| `tau_l` | `tau_l`, `tau_l_detrend`, `tau_l_hat_detrend`, `tau_l_hat_nodetrend` |
| `tau_k` | `tau_k`, `tau_k_detrend`, `tau_k_hat_detrend`, `tau_k_hat_nodetrend` |
| `tau_c` | `tau_c`, `tau_c_detrend`, `tau_c_hat_detrend`, `tau_c_hat_nodetrend` |
| `tau_d` | `tau_d`, `tau_d_detrend`, `tau_d_hat_detrend`, `tau_d_hat_nodetrend` |

Where `*_detrend` is the linearly-detrended level series and `*_hat_detrend` / `*_hat_nodetrend` are the two fractional-deviation alternatives.

**`1d_clean_for_matlab.py` picks the canonical version** via the `HAT_CHOICE` dict at the top of the script:

```python
The default choices reflect the diagnostic evidence in `1c_diagnose_data.py`:

- `tau_l`, `tau_c`, and `tau_d` show economically meaningful trends, so the default specification uses the linear-detrended variants.
- `tau_k` is dominated by cyclical and policy variation, including the 1980s decline and the TCJA cut, so the default specification uses the demean-only variant.
```

For each tax, the chosen variant is renamed to the canonical `tau_x_hat` in the MATLAB exports; the other variant and the level series are dropped. The `.mod` file always references `tau_l_hat`, `tau_k_hat`, etc., regardless of which underlying transformation was selected. **Switching specifications for a robustness check is a one-line edit.**

The default choices reflect the diagnostic evidence in `1c_diagnose_data.py` (CHECK 6b: economic-magnitude trend assessment):

- `tau_l`, `tau_d` trend strongly upward (rising payroll/social-insurance contributions for `tau_l`; rising income-tax progressivity for `tau_d` via the `tau_p` proxy). Detrending matches Jones (2002).
- `tau_k` movement is dominated by cyclical and policy variation (1980s drop, 2017 TCJA cut), not smooth drift. Demean-only preserves cyclical signal and matches Jones (2002).
- `tau_c` declines slowly over the postwar period as the services share of consumption rises. Since the model has no mechanism for this slow drift, the default specification uses the linear-detrended variant. The demean-only variant is kept for robustness.

`1c_diagnose_data.py` runs both a slope-significance test (CHECK 6) and an economic-magnitude assessment (CHECK 6b) on each tax level. CHECK 6b reports `drift_rel` (total drift as fraction of mean), `R²` (how well a linear trend fits), and a structural-break flag (residual variance asymmetry between halves of the sample). The combination distinguishes real trends from breaks and from significant-but-tiny drift.

---

## 4. Group A — Macro baseline observables

### 4.1 `chat` — real per-capita consumption

**Inputs**
- `PCND` — personal consumption expenditures, nondurable goods (FRED, BEA NIPA Table 1.1.5, SAAR)
- `PCESV` — personal consumption expenditures, services (FRED, BEA NIPA Table 1.1.5, SAAR)
- `GDPDEF` — GDP implicit price deflator (FRED, index)
- `CNP16OV` — civilian noninstitutional population age 16+ (FRED, monthly→quarterly mean, in thousands)

**Construction**
1. Nominal consumption = `PCND + PCESV`
2. Real consumption = nominal / (`GDPDEF` / 100)
3. Real per-capita = real / (population × 1000)
4. `chat` = `100 × (log(real per-capita) − linear trend)`

### 4.2 `ihat` — real per-capita investment

**Inputs**
- `GPDI` — gross private domestic investment (FRED, BEA NIPA Table 1.1.5, SAAR)
- `PCDG` — personal consumption expenditures, durable goods (FRED, BEA NIPA Table 1.1.5, SAAR)
- `GDPDEF`, `CNP16OV` (as above)

**Construction**
1. Nominal investment = `GPDI + PCDG` (consumer durables treated as investment, standard RBC/DSGE convention)
2. Same deflator, per-capita, log-linear-detrend pipeline as `chat`.

### 4.3 `ghat` — real per-capita government spending

**Inputs**
- `GCE` — government consumption expenditures and gross investment (FRED, BEA NIPA Table 1.1.5 Line 22, SAAR)
- `GDPDEF`, `CNP16OV`

**Construction**
Same pipeline as `chat`, with `GCE` as the nominal input.

### 4.4 `bhat` — real per-capita federal debt

**Inputs**
- `MVPHGFD027MNFRBDAL` — market value of privately held gross federal debt, Dallas Fed (FRED, monthly, $M)
- `GDPDEF`, `CNP16OV`

**Construction**
1. Monthly debt → quarterly **mean** (debt is a stock; mean represents average outstanding)
2. Real per-capita debt, log-linear-detrend ×100

Because debt is a stock and GDP is a SAAR flow, the ratio `bhat → b/y` has units of years of annual output, in line with how DSGEs typically express debt-to-GDP calibrations.

### 4.5 `pi` — gross inflation

`pi[t] = GDPDEF[t] / GDPDEF[t−1]`. Mean ≈ 1.008–1.010 in postwar US. No detrend.

### 4.6 `r` — federal funds rate

`r[t] = FEDFUNDS[t] / 400`, where `FEDFUNDS` is the monthly federal funds effective rate (annualized, percent), aggregated monthly→quarterly by mean. The `/400` converts annual percent to a quarterly decimal rate, matching the DSGE measurement equation `r_obs = R/400`.

### 4.7 `tau_l_hat` — labor tax rate (Jones 2002)

**Inputs (BEA NIPA, all SAAR, auto-rescaled to $bn)**
- `FIT` — federal personal current taxes (Table 3.2)
- `SIT` — state and local personal current taxes (Table 3.3)
- `W` — wages and salaries (Table 1.12)
- `PRI` — proprietors' income (Table 1.12)
- `Rental` — rental income of persons (Table 1.12)
- `CP` — corporate profits (Table 1.12)
- `NI` — net interest (Table 1.12)
- `CSI` — contributions for government social insurance (Table 3.1)
- `EC` — compensation of employees (Table 1.12)

**Construction**
1. Capital income: `CI = Rental + CP + NI + PRI/2`
2. Personal income tax rate: `tau_p = (FIT + SIT) / (W + PRI/2 + CI)`
3. Labor tax rate: `tau_l = (tau_p × (W + PRI/2) + CSI) / (EC + PRI/2)`

**Hat variants produced by `1b_build_data.py`**
- `tau_l_hat_detrend`: linear-detrend then fractional deviation `(tau_l_detrend − μ)/μ`
- `tau_l_hat_nodetrend`: demean-only fractional deviation `(tau_l − μ)/μ`

The labor tax rate has a documented upward trend from rising payroll/social-insurance contributions (`CSI` grew from ~3% of GDP in the 1950s to ~7% by the 2000s). The DSGE has no mechanism for this drift. Default `HAT_CHOICE["tau_l"] = "detrend"` follows Jones (2002).

### 4.8 `tau_k_hat` — capital tax rate (Jones 2002)

**Additional inputs**
- `CT_fed` — federal taxes on corporate income (Table 3.2)
- `CT_sl` — state and local taxes on corporate income (Table 3.3)
- `PT` — property taxes (Table 3.3, state and local)

**Construction**
1. `CT = CT_fed + CT_sl`
2. `tau_k = (tau_p × CI + CT + PT) / CI` — capital tax revenue over pre-tax capital income, following Jones (2002). The denominator is `CI` alone (which already includes corporate profits gross of corporate taxes), not `CI + CT + PT`.

**Hat variants produced by `1b_build_data.py`**
- `tau_k_hat_detrend`: linear-detrend then fractional deviation
- `tau_k_hat_nodetrend`: demean-only fractional deviation

Movement is dominated by cyclical variation (1980s drop, 2000s recovery, 2010s decline, 2020s rebound) and discrete policy events (most notably the 2017 TCJA cut from 35% to 21% federal rate). Default `HAT_CHOICE["tau_k"] = "nodetrend"` preserves cyclical signal that the model is supposed to capture.

Property tax data starts 1958Q1, which fixes the earliest possible sample start.

---

## 5. Group B — Extension observables

### 5.1 `lhat` — per-capita hours worked

**Input**
- `HOANBS` — Nonfarm Business Sector: Hours Worked, BLS index (2017=100), already quarterly

**Construction**
1. Hours per capita: `HOANBS / (population × 1000)` (units arbitrary because HOANBS is an index; only log-changes matter)
2. `lhat = 100 × (log(per-capita hours) − linear trend)`

The base year of the index (2017=100) is irrelevant because log-linear-detrend absorbs any constant scaling. This follows the Smets-Wouters (2007) convention.

### 5.2 `trhat` — real per-capita government transfers

**Input**
- `B087RC1Q027SBEA` — government social benefits to persons (FRED, BEA NIPA Table 3.1, SAAR, current $bn)

This is the broad transfer aggregate: Social Security, Medicare, Medicaid, unemployment insurance, veterans' benefits, etc.

**Construction**
Identical pipeline to `chat`/`ihat`/`ghat`/`bhat`: deflate by `GDPDEF`, divide by population, 100 × log-linear-detrend. Because transfers are a flow (SAAR), they sit on the same footing as `G` for any government-budget identity calculations.

### 5.3 `tau_c_hat` — consumption tax rate (Mendoza-Razin-Tesar 1994)

**Inputs**
- `TPI` — taxes on production and imports (BEA NIPA Table 3.1). Pre-2003 this category was called "indirect business taxes." Includes sales tax, excise tax, customs duties, business property taxes.
- `PT` — property taxes (BEA NIPA Table 3.3, state and local)
- `PCE_total` = `PCND + PCESV + PCDG`
- `GCE`

**Construction**

```
tau_c = (TPI − PT) / (PCE_total + GCE − (TPI − PT))
```

- Numerator isolates indirect taxes that fall on consumption (subtracting `PT` removes property tax, which falls on the capital base).
- Denominator is the pre-tax consumption base: total nominal consumption (private + government) minus the indirect tax we just isolated.

This formula crosses sources (BEA NIPA flows in the numerator, FRED PCE/GCE in the denominator), so unit consistency is critical. The `bea_nipa()` wrapper rescales NIPA values from millions to billions automatically.

**Hat variants produced by `1b_build_data.py`**
- `tau_c_hat_detrend`: linear-detrend then fractional deviation
- `tau_c_hat_nodetrend`: demean-only fractional deviation

The consumption tax rate declined ~25% over the postwar period (from ~6.7% to ~5.0%) as the services share of consumption rose and services are taxed less heavily than goods in most US states. Default `HAT_CHOICE["tau_c"] = "detrend"` follows the standard MRT (1994) operationalization.

### 5.4 `tau_d_hat` — dividend tax rate (proxy)

There is no clean quarterly time series of effective dividend tax rates. The script uses `tau_p` (the personal income tax rate from Section 4.7) as a proxy:

- **Exact pre-2003**: dividends were taxed as ordinary income, so `tau_p` is the effective dividend rate.
- **Upper bound post-2003 (JGTRRA)**: qualified dividends were taxed at preferential rates (15%, then 20% top bracket).

**Hat variants produced by `1b_build_data.py`**
- `tau_d_hat_detrend`: linear-detrend then fractional deviation
- `tau_d_hat_nodetrend`: demean-only fractional deviation

Default `HAT_CHOICE["tau_d"] = "detrend"`. `tau_d` (proxied by `tau_p`) inherits a clear upward drift from rising income-tax progressivity, which the DSGE has no mechanism for.

The `build_tau_d` function in `1b_build_data.py` is isolated, so swapping in a McGrattan-Prescott (2005), Gourio-Miao (2010), or constant-calibration alternative for the underlying level series is a one-line change. **Flag this proxy in the paper's data appendix.**

---

## 6. Group C — Tax policy series

### 6.1 `tau_itc_hat` — aggregate Investment Tax Credit rate

The federal ITC existed 1962–1986 (repealed by the Tax Reform Act of 1986). House-Shapiro's `ITCdat.txt` provides quarterly ITC rates by asset type for 1959Q1–2007Q1; rates are zero before 1962 and after 1986 within that file. Outside the file:

- 1958Q1–1958Q4: backfilled with `tau_itc = 0` (ITC didn't exist)
- 2007Q2–current: filled with `tau_itc = 0` (statutory, ITC repealed long since)

**Inputs**
- `raw/ITCdat.txt` — quarterly ITC rate by 38 asset types, 1959Q1–2007Q1 (HS replication)
- `raw/PQdat.txt` — quarterly nominal investment by 38 asset types (HS replication; used for type weights since the 38-type ITC rates align with the 38-type investment shares)

**Construction**
1. Type shares: `s_HS[t,m] = PQ_HS[t,m] / Σ_j PQ_HS[t,j]`
2. Aggregate ITC: `tau_itc[t] = Σ_m s_HS[t,m] × ITC[t,m]`
3. `tau_itc_hat = tau_itc − mean(tau_itc)` (level deviation)

**Calibration note:** set `tauic_ss = mean(tau_itc)` in the `.mod` file (printed by `1d_clean_for_matlab.py`). Because `tau_itc` is zero for ~75% of the sample, level deviation is the right form — fractional `(τ − μ)/μ` would explode.

**Why HS files for ITC but BEA files for `e_tau`?** ITC rates by asset type only exist in HS's compilation (the IRS published rates by IRS asset class, which HS crosswalked to their 38-type taxonomy). BEA Fixed Assets uses a different ~70-type taxonomy with no ITC information. So we use HS's matched (rates, weights) pair for ITC, and BEA's much longer time series for e_tau weights.

### 6.2 `e_tau_hat` — aggregate bonus depreciation expensing share

The DSGE has one aggregate capital good, but bonus depreciation eligibility is per asset type (depending on MACRS recovery period). We aggregate to a single quarterly series by weighting eligibility by each type's share of nominal investment:

```
eligible_share[t] = Σ_m (PQ_BEA[t,m] / Σ_j PQ_BEA[t,j]) × treat[m]
e_tau[t]          = bonus[t] × eligible_share[t]
e_tau_hat[t]      = e_tau[t] − mean(e_tau)        (level deviation)
```

**Inputs (produced by `1a_pre_prep.py`):**
- `raw/PQdat_BEA.csv` — quarterly nominal investment by detailed asset type, 1958Q1 to current. Built from BEA Fixed Assets Table 2.7 (annual values applied to all four quarters of each year). After dropping aggregate parent rows, ~70 leaf asset types remain.
- `raw/treat_BEA.csv` — 0/1 bonus-eligibility flag for each asset type, assigned by MACRS recovery period:
  - ≤20 years → 1 (equipment, software, R&D, entertainment originals, utility/telecom/farm/oil-gas quasi-structures)
  - >20 years → 0 (commercial structures, residential, land)
- A statutory quarterly bonus depreciation calendar `bonus[t]`, hard-coded in `1b_build_data.py`.

**Statutory bonus depreciation calendar**

| Period | `bonus[t]` | Legislation |
|---|---|---|
| Before 2001Q4 | 0 | (no bonus regime) |
| 2001Q4 – 2003Q2 | 0.30 | JCWAA 2002 (retroactive to 9/11/2001) |
| 2003Q3 – 2004Q4 | 0.50 | JGTRRA 2003 |
| 2005Q1 – 2007Q4 | 0 | bonus expired |
| 2008Q1 – 2010Q4 | 0.50 | Economic Stimulus Act 2008 + ARRA 2009 + SBJA 2010 |
| 2011Q1 – 2011Q4 | 1.00 | Tax Relief Act 2010 |
| 2012Q1 – 2013Q4 | 0.50 | American Taxpayer Relief Act 2012 |
| 2014Q1 – 2014Q4 | 0.50 | Tax Increase Prevention Act 2014 (retroactive) |
| 2015Q1 – 2017Q3 | 0.50 | PATH Act 2015 |
| 2017Q4 – 2022Q4 | 1.00 | TCJA (property placed in service after 9/27/2017) |
| 2023Q1 – 2023Q4 | 0.80 | TCJA phase-down |
| 2024Q1 – 2024Q4 | 0.60 | TCJA phase-down |
| 2025Q1 | 0.40 | pre-OBBBA effective date (1/20/2025) |
| 2025Q2 – 2026+ | 1.00 | OBBBA (July 2025), retroactive to 1/20/2025, made permanent |

The 2025Q1 = 0.40 entry is a clean approximation; OBBBA was retroactive to January 20, so the first 19 days of Q1 are at 40% and the rest at 100%. A more precise weight (~0.52, weighted by days) is one line to swap in if needed.

**Calibration note:** set `etau_ss` in the `.mod` file equal to `mean(e_tau)` over the estimation window (printed at the end of `1d_clean_for_matlab.py`'s run). Because `e_tau` is exactly zero for ~75% of the sample, a fractional deviation `(e_tau − μ)/μ` would explode — that's why this is a level deviation.

### 6.3 `tau_i_hat` and `tau_div_hat` — interest and dividend tax rates

Annual federal marginal tax rates on interest income and qualified dividends from the **McGrattan-TAXSIM** file (`raw/taxsim_fed_McGrattan.txt`), distributed by Ellen McGrattan and originally computed using NBER TAXSIM (Feenberg-Coutts 1993). Coverage 1960–2013, annual.

**Inputs**
- `Interest` column → `tau_i` (interest income tax rate)
- `Qualified Dividends` column → `tau_div` (dividend tax rate). Pre-2003 this equals the ordinary dividend column; post-2003 it captures the preferential tax treatment of qualified dividends. Before qualified-dividend treatment is relevant, the series should be interpreted as the TAXSIM dividend-income marginal tax rate.
**Construction**
1. Read annual values, convert percent to decimal.
2. Apply each year's value to all four quarters (1958Q1–1959Q4 and 2014Q1–current are NaN).
3. Demean using only observed (non-NaN) quarters; produce two variants parallel to the other tax rates:
   - `tau_x_hat_detrend`: linear-detrended fractional deviation. The trend is fit on observed quarters only and the recentering uses the observed mean.
   - `tau_x_hat_nodetrend`: demean-only fractional deviation, `(tau_x − mean) / mean`.

The default `HAT_CHOICE` in `1d_clean_for_matlab.py` uses `nodetrend` for `tau_i` (no clear monotonic trend in interest tax rates) and `detrend` for `tau_div` (clear downward trend from TRA 1986 and JGTRRA 2003).


---

## 7. BEA Fixed Assets reconstruction (Table 2.7)

`1a_pre_prep.py` pulls **BEA Fixed Assets Table 2.7** ("Investment in Private Fixed Assets, Equipment, Structures, and Intellectual Property Products by Type"; API code `FAAt207`). The script:

1. Downloads the full annual table (1901–current) via the BEA API.
2. Detects and drops aggregate parent rows. Two passes: (a) numerical sum-check (a row whose value approximately equals the sum of the next contiguous block of rows is an aggregate); (b) a backstop name list for known parent categories like "Equipment", "Structures", "Software", "Research and development", "Entertainment, literary, and artistic originals".
3. Disambiguates duplicate column names (BEA reuses "Other" under multiple parents) by appending the line number.
4. Assigns 0/1 eligibility flags by substring rules on the asset description, mapping MACRS recovery periods to bonus eligibility.
5. Expands annual values to quarterly by applying each year's value to all four quarters.
6. Writes `raw/PQdat_BEA.csv`, `raw/treat_BEA.csv`, `raw/asset_types.txt`.

**Differences from House-Shapiro (2008):**

| | House-Shapiro 2008 | This pipeline |
|---|---|---|
| Source | Unpublished BEA detailed asset files + IRS asset-class crosswalk | BEA Fixed Assets Table 2.7 (public API) |
| Frequency | Quarterly nominal investment by type | Annual, applied to 4 quarters of each year |
| Asset taxonomy | 38 types, 2007 vintage | ~70 leaf types, current BEA taxonomy |
| Eligibility flag | IRS asset-class crosswalk by type | MACRS recovery period rule (≤20yr → eligible) |
| Aggregation | `bonus × eligible_share` (same as ours) and a PV-weighted variant | `bonus × eligible_share` only |
| Coverage | 1959Q1 – 2007Q1 | 1958Q1 – current |

**What's lost:** within-year compositional shifts (annual data smoothed to quarterly), and the PV-weighted refinement that distinguishes 3-year from 20-year MACRS classes. Both are second-order for an aggregate DSGE observable.

**What's gained:** continuous coverage from 1958 to today with no splice or forward-fill, automatic updates as BEA releases new data, programmatic reproducibility through the public API.

**Validation against HS:** at 2007Q1 (the HS endpoint), our `eligible_share ≈ 0.54`, in line with the HS-derived value (~0.55). The eligible share rises slowly over time (≈0.47 in 1965, 0.51 in 1985, 0.54 in 2007, 0.56 in 2024), reflecting the long-run shift from structures-heavy to equipment- and software-heavy investment.

---

## 9. Implementation notes

**API keys.** `1b_build_data.py` and `1c_diagnose_data.py` read `FRED_KEY` and `BEA_KEY` from `raw/api_keys.txt`. Format: one `KEY=value` pair per line; comments start with `#`. Don't commit the file to a public repo.

**Population.** Default is `CNP16OV` from FRED (monthly, civilian noninstitutional population 16+, in thousands), averaged to quarterly. Alternative: BLS series `LNU00000000Q` (quarterly NSA, same CPS source) — set `POP_SOURCE = "CSV"` at the top of `1b_build_data.py` and place the file at `POP_CSV_PATH`.

**NIPA line numbers.** The script auto-discovers each line by description substring (e.g., "Personal current taxes", "Taxes on production and imports") and prints what it found. This is robust to BEA's periodic table renumbering.

**BEA back-data revisions.** BEA does periodic comprehensive revisions (most recent: 2023, 2018, 2013), which revise back-data. Pre-2008 numbers from the current BEA vintage will not perfectly match values that were available in 2014 or earlier. This is irreducible without using an archived data vintage.

---

## 10. Observed sample means

Reproduced sample means from a representative run:

| Series | Full sample (1958Q1–2019Q4) | 1958Q1–2008Q4 (Zubairy window) |
|---|---|---|
| `b̄/ȳ` | 0.355 | 0.302 |
| `ḡ/ȳ` | 0.205 | 0.209 |
| `τ̄_l` | 0.216 | 0.213 |
| `τ̄_k` | 0.416 | 0.429 |
| `τ̄_c` | 0.058 | 0.060 |
| `τ̄_d` (proxy) | 0.139 | 0.138 |
| `τ̄_itc` | 0.015 | 0.018 |
| `ē_τ` | 0.089 | 0.018 |

`1d_clean_for_matlab.py` prints the empirical means for both windows at the end of its run, so they can be pasted directly into the `.mod` file calibration block (`tau_l_bar`, `tau_k_bar`, `tau_itc_bar`, `e_tau_bar`, etc.).

---

## 11. Pipeline checklist before estimation

| Item | Where to verify |
|---|---|
| `raw/api_keys.txt` exists with `FRED_KEY` and `BEA_KEY` | manual |
| `1a_pre_prep.py` has been run; `raw/PQdat_BEA.csv` and `raw/treat_BEA.csv` exist | manual |
| Eligibility flags in `raw/asset_types.txt` look reasonable | manual one-time check |
| All observables present in `US_Data_Full.xlsx` with `year`, `quarter` keys | `1b_build_data.py` end-of-run summary |
| Every NIPA component sanity-check value sits in expected range (no `!!` flags) | `1b_build_data.py` "Sanity check" printout |
| Every tax rate within [0, 1] | `1c_diagnose_data.py` CHECK 1 |
| All hat series reject unit-root null (low-power caveat applies) | `1c_diagnose_data.py` CHECK 3 |
| Real-quantity hats look stationary around 0 | `1c_diagnose_data.py` CHECK 4 plot |
| No undetected trends in tax rates | `1c_diagnose_data.py` CHECK 5 plot + CHECK 6 + CHECK 6b |
| Debt/GDP ≈ 0.30 over 1959–2008 window | `1c_diagnose_data.py` CHECK 8 |
| `τ_d` proxy flagged in paper data appendix | manual |
| `HAT_CHOICE` set to intended specification before running `1d_clean_for_matlab.py` | `1d_clean_for_matlab.py` top |
| `tau_*_bar` values set in `.mod` file from data means | `1d_clean_for_matlab.py` end-of-run printout |
| `tauic_ss = mean(tau_itc)` and `etau_ss = mean(e_tau)` set in `.mod` file | manual |
| Robustness check: flip one entry of `HAT_CHOICE`, re-run, re-estimate | downstream |
