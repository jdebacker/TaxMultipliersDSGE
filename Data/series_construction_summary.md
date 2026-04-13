# Constructing the extended-model observables and policy series

## Big picture

Because the DSGE has one aggregate capital good, but House and Shapiro provide type-specific objects (ITC by type and bonus eligibility by type), an aggregation step is required to obtain quarterly aggregate policy rates.

We construct aggregate policy rates using nominal investment weights:
\[
\tau_{ic,t}=\sum_m s_{t,m}\,ITC_{t,m},\qquad
e_{\tau,t}=bonus_t \cdot \sum_m s_{t,m}\,treat_m
\]
where \(s_{t,m}\) are nominal investment shares by type from PQ (House and Shapiro) or BEA. Without weights, each type would be treated as equally important, which is not appropriate.

---

## 1) Debt series: `bhat`

### Raw inputs
- Monthly Dallas Fed market value privately held gross federal debt: `MVPHGFD027MNFRBDAL`
- Quarterly GDP deflator: `GDPDEF` (index)
- Quarterly civilian non-institutional population age 16+: `LNU00000000Q` (thousands)

### Steps
1. Convert monthly debt to quarterly by averaging the three months within each quarter (sum also works — the intercept in the linear trend absorbs the log(3) difference — but averaging is more conventional).
2. Deflate to real terms by dividing by \((GDPDEF/100)\).
3. Convert to per-capita by dividing by \((POP\_THOUS \times 1000)\).
4. Make stationary by taking logs and removing a linear time trend on the target sample window (1959Q1–2013Q4), then multiplying residuals by 100.
5. Output a quarterly series `bhat` with `year` and `quarter`.

---

## 2) ITC series: `tau_ic` and `tau_ic_hat`

### Raw inputs
- `ITCdat.txt`: quarterly ITC rate by type (38 columns, House and Shapiro)
- `PQdat.txt`: quarterly nominal investment spending by type (38 columns, House and Shapiro), used for weights

### Steps
1. Attach a quarter index starting at 1959Q1 to both `ITCdat.txt` and `PQdat.txt`.
2. Extend to 2013Q4:
   - ITC rates are filled with 0 after the House and Shapiro sample ends (ITC was repealed in 1986, so ITCdat already has zeros from ~1987 onward; filling with 0 post-2007 is correct).
   - Investment weights are forward-filled after the House and Shapiro sample ends (placeholder until BEA weights are added post-2007; see Section 5 below).
3. Compute nominal type shares each quarter:
   \[
   s_{t,m} = \frac{PQ_{t,m}}{\sum_j PQ_{t,j}}.
   \]
4. Aggregate the ITC rate:
   \[
   \tau_{ic,t}=\sum_m s_{t,m}\,ITC_{t,m}.
   \]
5. Construct the hat series as a demeaned level (rate series, no time trend):
   - `tau_ic_hat = tau_ic - mean(tau_ic)`.

### Calibration note
When adding `tau_ic` as an observable in the extended mod file, set `tau_ic_bar = mean(tau_ic)` so the model steady state matches the data mean. Since `tau_ic_hat = tau_ic - mean(tau_ic)`, the measurement equation `tau_ic_hat_obs = tau_ic - tau_ic_bar` is consistent. This is especially important because `tau_ic` is near zero for most of the sample — using a fractional deviation `(tau_ic - bar)/bar` would blow up.

---

## 3) Depreciation policy series: `e_tau` and `e_tau_hat`

### Raw inputs
- `PQdat.txt`: quarterly nominal investment spending by type (weights)
- `treat.txt`: 0/1 indicator by type for bonus-eligible property (MACRS recovery period ≤ 20 years)
- A statutory quarterly bonus depreciation calendar, coded as `bonus_t` (values such as 0, 0.30, 0.50, 1.00)

### Bonus depreciation calendar

| Period | `bonus_t` | Legislation | Notes |
|--------|-----------|-------------|-------|
| Before 2001Q4 | 0 | — | No bonus depreciation |
| 2001Q4–2003Q2 | 0.30 | JCWAA (signed March 2002, retroactive to 9/11/2001) | House & Shapiro Figure 4 shows clear investment effects starting 2001Q4; retroactivity means property placed in service after 9/11/2001 qualified |
| 2003Q3–2004Q4 | 0.50 | JGTRRA (signed May 2003, expires Jan 1, 2005) | Property must be placed in service before January 1, 2005 |
| 2005Q1–2007Q4 | 0 | — | Bonus expired; Section 179 expensing extended separately but not modeled here |
| 2008Q1–2010Q4 | 0.50 | Economic Stimulus Act of 2008 + ARRA 2009 + SBJA 2010 | |
| 2011Q1–2011Q4 | 1.00 | Tax Relief, Unemployment Insurance Reauthorization, and Job Creation Act of 2010 | 100% first-year expensing |
| 2012Q1–2013Q4 | 0.50 | American Taxpayer Relief Act of 2012 | |

**In `clean_data.py`:**
```python
bonus = pd.Series(0.0, index=df_inv.index)
bonus.loc[pd.period_range("2001Q4", "2003Q2", freq="Q")] = 0.30
bonus.loc[pd.period_range("2003Q3", "2004Q4", freq="Q")] = 0.50
bonus.loc[pd.period_range("2008Q1", "2010Q4", freq="Q")] = 0.50
bonus.loc[pd.period_range("2011Q1", "2011Q4", freq="Q")] = 1.00
bonus.loc[pd.period_range("2012Q1", "2013Q4", freq="Q")] = 0.50
```

### Steps
1. Attach a quarter index starting at 1959Q1 to `PQdat.txt` and extend to 2013Q4 by forward-filling weights (placeholder until BEA weights are added post-2007; see Section 5 below).
2. Compute nominal type shares \(s_{t,m}\) from PQ.
3. Compute the eligible investment share each quarter:
   \[
   eligible\_share_t = \sum_m s_{t,m}\,treat_m.
   \]
4. Convert the statutory bonus rate into an aggregate expensing share:
   \[
   e_{\tau,t} = bonus_t \cdot eligible\_share_t.
   \]
5. Construct the hat series as a demeaned level:
   - `e_tau_hat = e_tau - mean(e_tau)`.

### Calibration note
When adding `e_tau` as an observable in the extended mod file, set `e_tau_bar = mean(e_tau)` so the model steady state matches the data mean. The measurement equation `e_tau_hat_obs = e_tau - e_tau_bar` is then consistent with the data construction. If `e_tau_bar` is left at 0 but the data was demeaned around a nonzero mean, the Kalman filter sees a permanent level mismatch. In practice: run `build_e_tau`, print `e_tau.mean()`, and use that value.

---

## 4) House and Shapiro type ordering (m = 1..38)

The 38 capital types and their ordering are taken from the column headers of `MACRO_ITC_aer.xls` (sheet `ITC`):
- The first two columns are `YEAR` and `QTR`.
- The next 38 columns define the type names and their order, which aligns with `ITCdat.txt`, `PQdat.txt`, and the 38-length vectors `treat.txt`, `dhat.txt`, and `edelta.txt`.

Ordered list (m = 1..38):
1. Computers and peripheral equipment
2. Software
3. Communication equipment
4. Instruments (med)
5. Instruments (nonmed)
6. Photocopy and related equipment
7. Other office equipment
8. Other fabricated metal products
9. Steam engines
10. Internal combustion engines
11. Metalworking machinery
12. Special industry machinery, n.e.c.
13. General industrial, including materials handling, equipment
14. Electrical transmission, distribution, and industrial apparatus
15. Trucks, buses, and truck trailers
16. Autos
17. Aircraft
18. Ships and boats
19. Railroad equipment
20. Farm tractors
21. Agricultural machinery, except tractors
22. Construction tractors
23. Construction machinery, except tractors
24. Mining and oilfield machinery
25. Service industry machinery
26. Commercial incl Office
27. Hospital and institutional buildings
28. Industrial buildings
29. Electric light and power (structures)
30. Gas (structures)
31. Telecommunications
32. Petroleum and natural gas
33. Other mining exploration
34. Religious buildings
35. Educational buildings
36. Railroads
37. Farm related buildings and structures
38. Single Family Res

**Bonus eligibility:** Types 1–25 are equipment (MACRS recovery period ≤ 20 years) and are generally eligible. Types 26–28 and 34–35 are structures with 39-year recovery (ineligible). Types 29–33 and 36–37 are "quasi-structures" with 15- or 20-year recovery (eligible). Type 38 (residential) is ineligible. The exact 0/1 classification is in `treat.txt`.

---

## 5) Updating investment weights post-2007

### Why post-2007 weights matter
`PQdat.txt` covers 1959Q1–2007Q1 (193 quarterly observations, per House and Shapiro README). Forward-filling PQ weights freezes investment composition after 2007Q1, which makes `eligible_share` constant and therefore produces piecewise-constant `e_tau` whenever `bonus_t` is constant.

Post-2007 weights are **less important for ITC** in practice because the aggregate ITC rate is essentially zero after the ITC era. Post-2007 weights **matter much more for `e_tau`** because the eligible investment share depends on composition, which shifts over time (e.g., software grows, structures shrink).

### What to download from BEA
Source: https://www.bea.gov/data/investment-fixed-assets/by-type

You need two tables, both in **current dollars** (nominal), annual:
- **Table 5.5.5**: Private Fixed Investment in Equipment by Type — covers types 1–25 (equipment)
- **Table 5.4.5**: Private Fixed Investment in Structures by Type — covers types 26–38 (structures, quasi-structures, residential)

Download years 2007–2013 (or wider if extending the sample).

### Mapping BEA types to the 38 House and Shapiro types

Most are direct matches. Known tricky cases:
- **Instruments (types 4–5):** BEA may not split medical vs. nonmedical instruments the same way. If BEA reports a single "instruments" category, split using the PQ ratio from the last available year (2007Q1).
- **Steam engines (type 9):** May be merged into a broader BEA category. Use PQ ratio to allocate.
- **Railroad structures (type 36):** BEA reclassified these post-1997 into "land improvements" (primarily railroads). Use that line.
- **Single Family Res (type 38):** Use BEA residential fixed investment.

### Steps
1. Download annual current-dollar investment by detailed asset type from BEA Tables 5.5.5 and 5.4.5 (2007–2013).
2. Map BEA asset labels to the 38 House and Shapiro types. For types that don't map 1:1, allocate using PQ shares from 2007Q1 as a bridge.
3. Compute annual shares \(s_{y,m}\) and apply the same shares to all four quarters in year \(y\).
4. Replace the forward-filled PQ weights with BEA-based weights for 2007Q2–2013Q4.
5. Recompute `e_tau` and `e_tau_hat` (and `tau_ic`/`tau_ic_hat`, though the effect is negligible since ITC is zero).

---

## 6) Notes on additional depreciation objects
- `dhat.txt` and `Dep_sched.txt` are only needed if the approach is expanded to construct a present-value depreciation allowance object using detailed depreciation schedules and discounting assumptions (i.e., the Hall-Jorgenson \(z^m\) from House and Shapiro eq. 10).
- Under the current DSGE implementation, the key observed depreciation-policy input is the expensing share `e_tau`, while the tax-basis evolution and depreciation term \(\delta_\tau k^\tau_{t-1}\) are model-implied given \(\delta_\tau\) and `e_tau`.
- `edelta.txt` contains 38 economic depreciation rates \(\delta^m\) by type (annual). These are not needed for the aggregate DSGE but would be needed for any type-level analysis following House and Shapiro's structural estimation.

---

## 7) Checklist before estimation

| Item | Check |
|------|-------|
| `bhat`: verify sample period matches other observables (1959Q1–2013Q4, 220 obs) | ☐ |
| `tau_ic_hat`: confirm ITC is zero post-1986 (sanity check) | ☐ |
| `e_tau_hat`: confirm bonus calendar matches Section 3 table above | ☐ |
| `e_tau_bar` in mod file = `mean(e_tau)` from data | ☐ |
| `tau_ic_bar` in mod file = `mean(tau_ic)` from data | ☐ |
| `US_Data_forMatlab_old.xlsx` has `year` and `quarter` columns for merge | ☐ |
| Post-2007 BEA weights downloaded and mapped (or accepted forward-fill as placeholder) | ☐ |
| Merged output has no NaN in any observable column | ☐ |
