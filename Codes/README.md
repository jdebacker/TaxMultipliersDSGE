# Codes

Estimation and multiplier code for "Tax Multipliers in a DSGE Model". Requires MATLAB and **Dynare 5.5** (other Dynare versions give different chains from the same files).

## Folders

| folder | model | estimated parameters | role in the paper |
|---|---|---|---|
| `tax_dsge` | fiscal rules with debt and output feedback | 49 | baseline |
| `tax_dsge_ar1` | AR(1) fiscal rules, only transfers respond to debt | 34 | Section on financing, appendix |
| `tax_dsge_kgoods` | baseline rules + capital-goods friction (epsK = 4) | 49 | Section on capital goods, appendix |
| `tax_dsge_kgoods_ar1` | AR(1) rules + capital-goods friction | 34 | appendix |
| `tax_dsge_uk` | UK, AR(1) rules, no investment tax credit | 32 | UK section, appendix |

`old/` holds earlier versions of the model file and is not used.

Each US folder contains the model (`tax_dsge_est.mod`), the data (`US_Data_Matlab_1961_2013.mat`, 1961Q1-2013Q4), the starting point for the mode search (`start_mode.mat`), the mode and proposal covariance used by the MCMC (`tax_dsge_est_mode_cal.mat`), the two post-estimation scripts, and the results reported in the paper (`table_multipliers_pvown.csv`, `table_self_financing.csv`, `post_pvown.csv`). The UK folder has the same structure with `tax_dsge_uk.mod`, `UK_Data_Matlab_tax.mat` (1997Q2-2016Q4), `start_uk.mat`, `tax_dsge_uk_mode_cal.mat`, and a driver `run_uk.m`.

The MCMC chains are not included (several hundred MB per model). They are recreated by step 1 below, or can be sent on request.

## Stages

The estimation steps are selected with a Dynare macro variable in the `.mod` file:

| command | what it does |
|---|---|
| `dynare tax_dsge_est noclearall` | steady state and model checks only |
| `dynare tax_dsge_est noclearall -DSTAGE=1` | posterior mode search (newrat), starting from `start_mode.mat`; only needed to redo the mode |
| `dynare tax_dsge_est noclearall -DSTAGE=2` | MCMC, 4 chains of 1,000,000 draws, from `tax_dsge_est_mode_cal.mat` (about 5 hours) |
| `dynare tax_dsge_est noclearall -DSTAGE=5` | reload existing chains without sampling |

## Reproducing the results

For each US model, from inside its folder:

1. `dynare tax_dsge_est noclearall -DSTAGE=2` (estimation; skip if the chains are already in `tax_dsge_est\metropolis\`)
2. `dynare tax_dsge_est noclearall -DSTAGE=5` (only in a new MATLAB session; after step 1 the workspace is already loaded)
3. `compute_multipliers` writes `table_multipliers_pvown.csv`: present-value own-revenue multipliers at H = 1, 4, 12, 20 (Tables 3, 4, 6), and `table_self_financing.csv`: the share of self-financing draws
4. `post_estimation` writes `post_pvown.csv`: re-solves on the same draws with sigma_u, b_g and, in the friction models, epsK changed (Table 5)

For the UK, from inside `tax_dsge_uk`: `run_uk('estimate')` (mode search and MCMC), then `run_uk('post')` (reload, `compute_multipliers`, `post_estimation`; Table 7).

## Figures

Once all five models have chains, run from this folder:

```
make_paper_figures
```

It reads the five sets of chains without modifying them, re-solves the model on the same 1,000 posterior draws used for the tables, and writes to `paper_figures\`: impulse responses for the baseline model, multipliers by horizon (baseline with bands, the four U.S. models, UK against U.S. AR(1)), and posterior summaries for every model. The log compares its medians with each folder's `table_multipliers_pvown.csv`. `make_paper_figures('plot')` redraws the figures from the saved results without re-solving.

## Expensing multiplier denominator

The expensing own-cost denominator measures net revenue lost relative to keeping `e_tau = etau_ss` on the same simulated investment and capital-tax-rate paths. The reporting state `k_tau_etau` tracks the difference in tax basis:

```
k_tau_etau = (1 - delta_tau) * k_tau_etau(-1) - (e_tau - etau_ss) * i
spend_etau = tau_k * ((e_tau - etau_ss) * i + delta_tau * k_tau_etau(-1))
```

In capital-goods models, investment expenditure is `PI * i`. A rise in expensing creates a positive cost today and a negative tax-basis difference; smaller future depreciation deductions then recover revenue and reduce the net cost. This comparison excludes changes in taxable rental income and does not attribute changes in deductions from investment or capital tax rates alone to expensing. Both reporting variables are zero in steady state. To first order the main denominator depends on steady-state investment and capital tax rates, the simulated expensing path, and its accumulated effect on tax basis. The mechanical denominator uses the same recovery accounting with the instrument's AR(1) path, omitting endogenous fiscal feedback.

These reporting equations do not change `k_tau`, `rev_k`, or the government budget, which already account for expensing and depreciation. The total-revenue robustness calculation retains full revenue accounting. Numerator and denominator are both cumulated through the reported horizon; recovery after that horizon is not included.

Existing CSV tables and saved figures use the previous expensing denominator and must be regenerated with `compute_multipliers`, `post_estimation`, and `make_paper_figures`. Existing parameter draws can be reused; re-estimation is not required.


## Self-financing tax cuts and incentives

A tax cut, ITC increase, or expensing increase is self-financing in a given posterior draw and horizon when its discounted cumulative GDP response is positive and its discounted cumulative **total tax revenue** response is nonnegative. Zero is evaluated with a numerical tolerance of `1e-10` in model units. All changes are oriented in the expansionary policy direction, including reversing the simulated positive rate innovation for tax cuts. Spending and transfers are not classified.

The main table, `table_multipliers_pvown.csv`, keeps the finite own-revenue ratio for every draw, since it measures output per dollar of the instrument's own cost. `compute_multipliers` also writes:

- `table_self_financing.csv`: for each tax instrument and horizon, the share of posterior draws that are self-financing.
- `table_multipliers_pvown_inf.csv`: the own-revenue table with the `+Inf` reporting override applied to self-financing draws, draw by draw, before medians and intervals are computed.

The mechanical and total-revenue tables keep the `+Inf` override. `post_estimation` reports finite ratios in `post_pvown.csv` and saves the total-revenue response (`NET_COST`) in `post_raw.mat`, so the override or the shares can be recomputed for the sensitivity settings. `make_paper_figures` plots the finite ratios.

Intentional infinities survive the finite-ratio clipping rule and are written as `Inf` in CSVs and console tables. Interpolated percentiles touching `+Inf` remain infinite. Paired percentage changes involving infinite levels are undefined and use finite pairs only.

These changes affect reporting only and require no re-estimation. Run `addpath('Codes/tests'); test_fiscal_reporting` in MATLAB or Octave to check the reporting rules.
