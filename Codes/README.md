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

Each US folder contains the model (`tax_dsge_est.mod`), the data (`US_Data_Matlab_1961_2013.mat`, 1961Q1-2013Q4), the starting point for the mode search (`start_mode.mat`), the mode and proposal covariance used by the MCMC (`tax_dsge_est_mode_cal.mat`), the two post-estimation scripts, and the results reported in the paper (`table_multipliers_pvown.csv`, `post_pvown.csv`). The UK folder has the same structure with `tax_dsge_uk.mod`, `UK_Data_Matlab_tax.mat` (1997Q2-2016Q4), `start_uk.mat`, `tax_dsge_uk_mode_cal.mat`, and a driver `run_uk.m`.

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
3. `compute_multipliers` writes `table_multipliers_pvown.csv`: present-value own-revenue multipliers at H = 1, 4, 12, 20 (Tables 3, 4, 6)
4. `post_estimation` writes `post_pvown.csv`: re-solves on the same draws with sigma_u, b_g and, in the friction models, epsK changed (Table 5)

For the UK, from inside `tax_dsge_uk`: `run_uk('estimate')` (mode search and MCMC), then `run_uk('post')` (reload, `compute_multipliers`, `post_estimation`; Table 7).

## Figures

Once all five models have chains, run from this folder:

```
make_paper_figures
```

It reads the five sets of chains without modifying them, re-solves the model on the same 1,000 posterior draws used for the tables, and writes to `paper_figures\`: impulse responses for the baseline model, multipliers by horizon (baseline with bands, the four U.S. models, UK against U.S. AR(1)), and posterior summaries for every model. The log compares its medians with each folder's `table_multipliers_pvown.csv`. `make_paper_figures('plot')` redraws the figures from the saved results without re-solving.
