import numpy as np
import pandas as pd

#%%
# =============================================================================
# clean debt  -> bhat = 100 * (log real per-capita debt detrended linearly)

def build_bhat(debt_path, gdpdef_path, pop_path, start_q="1959Q1", end_q="2008Q4"):
    # Load raw series (Dallas Fed debt is monthly, deflator/pop are quarterly)
    debt_m = pd.read_excel(debt_path, sheet_name="Monthly")
    gdpdef_q = pd.read_excel(gdpdef_path, sheet_name="Quarterly")
    pop_q = pd.read_csv(pop_path)

    debt_m["observation_date"] = pd.to_datetime(debt_m["observation_date"])
    gdpdef_q["observation_date"] = pd.to_datetime(gdpdef_q["observation_date"])

    # Monthly -> quarterly (Zubairy: sum monthly series within each quarter)
    debt_q = (debt_m
              .assign(QPER=debt_m["observation_date"].dt.to_period("Q"))
              .groupby("QPER", as_index=False)["MVPHGFD027MNFRBDAL"].sum()
              .rename(columns={"MVPHGFD027MNFRBDAL": "DEBT_BIL"}))
    
    # GDP deflator is an index, convert to quarter period
    gdpdef_q = (gdpdef_q
                .assign(QPER=gdpdef_q["observation_date"].dt.to_period("Q"))
                [["QPER", "GDPDEF"]])

    # Population is in thousands, convert to quarter period
    pop_q = pop_q[pop_q["Value"].ne("-")].copy()
    pop_q["Value"] = pop_q["Value"].astype(float)
    pop_q["quarter"] = pop_q["Period"].str.extract(r"Q0?([1-4])").astype(int)
    pop_q["QPER"] = pd.PeriodIndex.from_fields(year=pop_q["Year"].astype(int), quarter=pop_q["quarter"], freq="Q")
    pop_q = pop_q[["QPER", "Value"]].rename(columns={"Value": "POP_THOUS"})

    # Merge to quarterly panel and keep target window
    df = (debt_q.merge(gdpdef_q, on="QPER", how="inner")
               .merge(pop_q, on="QPER", how="inner")
               .sort_values("QPER")
               .reset_index(drop=True))
    df = df[(df["QPER"] >= pd.Period(start_q, "Q")) & (df["QPER"] <= pd.Period(end_q, "Q"))].reset_index(drop=True)

    # Real per-capita debt: divide by (GDPDEF/100) and population in persons
    real_debt = df["DEBT_BIL"] / (df["GDPDEF"] / 100.0)          # billions of real dollars
    pc_real_debt = real_debt / (df["POP_THOUS"] * 1000.0)        # billions of real dollars per person

    # Stationary series: log, remove time trend, scale by 100
    log_b = np.log(pc_real_debt.to_numpy())
    t = np.arange(len(df))
    X = np.column_stack([np.ones(len(df)), t])
    beta = np.linalg.lstsq(X, log_b, rcond=None)[0]
    bhat = 100.0 * (log_b - X @ beta)
    
    # from statsmodels.tsa.filters.hp_filter import hpfilter
    # log_b = np.log(pc_real_debt.to_numpy())
    # _, trend = hpfilter(log_b, lamb=1600)
    # bhat = 100.0 * (log_b - trend)
    
    out = pd.DataFrame({
        "year": df["QPER"].dt.year.astype(int),
        "quarter": df["QPER"].dt.quarter.astype(int),
        "bhat": bhat
    })
    
    return out


#%%
# =============================================================================
# clean ITC  -> tau_ic (level) and tau_ic_hat (demeaned level)

def build_tau_ic(itc_path, pq_path, start_q="1959Q1", end_q="2013Q4"):
    # ITC rates by type and nominal spending by type (weights), both 193x38 in HS replication
    df_itc = pd.read_csv(itc_path, sep=r"\s+", header=None)
    df_itc.columns = [f"itc_{i+1:02d}" for i in range(df_itc.shape[1])]

    df_inv = pd.read_csv(pq_path, sep=r"\s+", header=None)
    df_inv.columns = [f"inv_{i+1:02d}" for i in range(df_inv.shape[1])]

    # Attach quarter index (HS panel starts 1959Q1), extend to end_q
    q0 = pd.Period("1959Q1", freq="Q")
    df_itc["QPER"] = [q0 + i for i in range(len(df_itc))]
    df_inv["QPER"] = [q0 + i for i in range(len(df_inv))]

    raw_end_itc = df_itc["QPER"].max()
    raw_end_inv = df_inv["QPER"].max()
    raw_end = min(raw_end_itc, raw_end_inv)

    idx = pd.period_range(start_q, end_q, freq="Q")
    data_filled = (idx > raw_end).astype(int)

    df_itc = df_itc.set_index("QPER").reindex(idx).fillna(0.0)
    df_inv = df_inv.set_index("QPER").reindex(idx).ffill()

    # Spending shares by type, then aggregate ITC rate
    shares = df_inv.div(df_inv.sum(axis=1), axis=0)
    shares.columns = df_itc.columns
    tau_ic = (df_itc * shares).sum(axis=1)

    # Keep "hat" as deviation from mean (rate series, no artificial trend)
    tau_ic_hat = tau_ic.to_numpy() - tau_ic.mean()

    out = pd.DataFrame({
        "year": tau_ic.index.year.astype(int),
        "quarter": tau_ic.index.quarter.astype(int),
        "data_filled_tau_ic": data_filled,
        "tau_itc": tau_ic.to_numpy(),
        "tau_itc_hat": tau_ic_hat
    })
    return out


#%%
# =============================================================================
# clean depreciation  -> e_tau (level) and e_tau_hat (demeaned level)

def build_e_tau(pq_path, treat_path, start_q="1959Q1", end_q="2013Q4"):
    # Nominal spending by type gives weights, treat is 0/1 for bonus-eligible types
    df_inv = pd.read_csv(pq_path, sep=r"\s+", header=None)
    df_inv.columns = [f"inv_{i+1:02d}" for i in range(df_inv.shape[1])]

    treat = pd.read_csv(treat_path, sep=r"\s+", header=None).iloc[:, 0].to_numpy()
    treat = pd.Series(treat, index=df_inv.columns)

    q0 = pd.Period("1959Q1", freq="Q")
    df_inv["QPER"] = [q0 + i for i in range(len(df_inv))]

    raw_end = df_inv["QPER"].max()

    idx = pd.period_range(start_q, end_q, freq="Q")
    data_filled = (idx > raw_end).astype(int)

    df_inv = df_inv.set_index("QPER").reindex(idx).ffill()

    shares = df_inv.div(df_inv.sum(axis=1), axis=0)
    eligible_share = shares.mul(treat, axis=1).sum(axis=1)

    # Statutory bonus rate by quarter (edit windows if needed)
    bonus = pd.Series(0.0, index=df_inv.index)
    bonus.loc[pd.period_range("2001Q4", "2003Q2", freq="Q")] = 0.30
    bonus.loc[pd.period_range("2003Q3", "2004Q4", freq="Q")] = 0.50
    bonus.loc[pd.period_range("2008Q1", "2010Q4", freq="Q")] = 0.50
    bonus.loc[pd.period_range("2011Q1", "2011Q4", freq="Q")] = 1.00
    bonus.loc[pd.period_range("2012Q1", "2013Q4", freq="Q")] = 0.50

    # Aggregate expensing share: statutory rate times eligible investment share
    e_tau = (bonus * eligible_share).to_numpy()
    e_tau_hat = e_tau - e_tau.mean()

    out = pd.DataFrame({
        "year": idx.year.astype(int),
        "quarter": idx.quarter.astype(int),
        "data_filled_e_tau": data_filled,
        "e_tau": e_tau,
        "e_tau_hat": e_tau_hat
    })
    return out


#%%
# =============================================================================
# run and merge into main dataset

base = r"C:\Users\Owner\OneDrive - University of South Carolina\economics 1\5 tax multipliers in a DSGE model\Data\raw"
output = r"C:\Users\Owner\OneDrive - University of South Carolina\economics 1\5 tax multipliers in a DSGE model\Data"

df_b = build_bhat(
    debt_path=rf"{base}\MVPHGFD027MNFRBDAL.xlsx",
    gdpdef_path=rf"{base}\GDPDEF.xlsx",
    pop_path=rf"{base}\LNU00000000Q.csv",
)

df_itc = build_tau_ic(
    itc_path=rf"{base}\ITCdat.txt",
    pq_path=rf"{base}\PQdat.txt",
)

df_e = build_e_tau(
    pq_path=rf"{base}\PQdat.txt",
    treat_path=rf"{base}\treat.txt",
)

df = pd.read_excel(rf"{base}\US_Data_forMatlab_old.xlsx")
df["year"] = df["year"].astype(int)
df["quarter"] = df["quarter"].astype(int)

df = df.merge(df_b, on=["year", "quarter"], how="left")
df = df.merge(df_itc, on=["year", "quarter"], how="left")
df = df.merge(df_e, on=["year", "quarter"], how="left")

# cols = ["chat", "ihat", "lhat", "ghat", "trhat", "bhat", "tau_ic_hat", "e_tau_hat"]
# print(df[cols].describe().T[["mean", "std", "min", "max"]])


# ====================
# export full data
df = df.rename(columns={'r': 'r_obs', 'pi': 'pi_obs', 'xhat': 'trhat'})

df.to_excel(rf"{output}\US_Data_Matlab_year_quarter.xlsx", index=False)

# ====================
# export two versions (1959Q1 to 2008Q4), select columns
df["QPER"] = pd.PeriodIndex.from_fields(year=df["year"], quarter=df["quarter"], freq="Q")
df_samp = df[(df["QPER"] >= pd.Period("1959Q1", freq="Q")) & (df["QPER"] <= pd.Period("2008Q4", freq="Q"))].copy()

# your full dataset
cols_all = ['chat', 'ihat', 'r_obs', 'pi_obs', 'ghat',
            'bhat', 'tau_k_hat', 'tau_l_hat',
            'lhat', 'trhat', 'tau_c_hat', 'tau_d_hat',
            'tau_itc_hat', 'e_tau_hat']

# Zubairy baseline observables (replication)
cols_zub = ['chat', 'ihat', 'r_obs', 'pi_obs', 'ghat',
            'bhat', 'tau_k_hat', 'tau_l_hat']

df_all = df_samp[cols_all].copy()
df_zub = df_samp[cols_zub].copy()

df_zub.to_excel(rf"{output}\US_Data_Matlab_Zubiary.xlsx", index=False)
df_all.to_excel(rf"{output}\US_Data_Matlab.xlsx", index=False)

from scipy.io import savemat
savemat(rf"{output}\US_Data_Matlab_Zubiary.mat", {c: df_zub[c].to_numpy() for c in df_zub.columns})
savemat(rf"{output}\US_Data_Matlab.mat", {c: df_all[c].to_numpy() for c in df_all.columns})




#%%

from scipy.io import loadmat
mat = loadmat(rf"{output}\US_Data_Matlab_Zubiary.mat")
arr = mat['bhat'].flatten()
print(f"bhat  mean={arr.mean():.4f}  std={arr.std():.4f}")
# Target: std should be around 8-12, not 26

#%%










