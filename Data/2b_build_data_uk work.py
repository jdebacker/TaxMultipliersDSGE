"""
1b_build_data_uk.py
===================
UK counterpart to 1b_build_data.py. ONS-quarterly tax-rate version. Pulls quarterly UK macro and fiscal data
from the ONS Time Series API and Bank of England Interactive Statistical
Database, constructs the same observable transformations as the US pipeline,
and writes UK_Data_Full.xlsx.

What's covered (mirrors US pipeline where the UK data are comparable):
    chat, ihat, ghat, bhat, lhat, trhat   -- real per-capita levels, log-linear-detrended
    pi_obs                                 -- gross GDP-deflator inflation
    r_obs                                  -- Bank Rate / 400 (quarterly decimal)
    tau_l, tau_k, tau_c                    -- UK average tax rates using ONS national accounts
                                              and public-sector finance data, with both
                                              detrended and nodetrend hat variants.

What's NOT covered (UK-specific data limits or different tax regime):
    tau_d, tau_i                           -- no clean UK TAXSIM equivalent.
    tau_itc, e_tau                         -- UK capital allowances are not directly comparable
                                              to US ITC + bonus depreciation.

Coverage: 1955Q1 to current. Some series (hours worked) start later (1971Q1).
The build_*() functions return NaN outside their available range; the merge
preserves the longest available window per series.
"""

#%%
# =============================================================================
# imports and config

import os
import time
import requests
import pandas as pd
import numpy as np

OUT_DIR = r"./output"
RAW_DIR = r"./raw"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(RAW_DIR, exist_ok=True)


# =============================================================================
# X-13ARIMA-SEATS seasonal adjustment helper
#
# ONS Public Sector Finances (PUSF) tax-revenue series and UKEA household
# property income/rental series are published NSA only. We seasonally adjust
# them with X-13 before they enter tax-rate construction so that
# (NSA revenue) / (SA base) doesn't carry numerator seasonality into the
# tax-rate hats.
#
# X-13 binary: download from https://www.census.gov/data/software/x13as.html
# Place the executable on PATH, or set the X13PATH environment variable.
# (statsmodels searches both.) On Windows: x13as.exe; on macOS/Linux: x13as.
#
# If the binary is not found, this raises a clear ImportError pointing to
# the install instructions, rather than silently falling back to a different
# method.

def _stl_sa(nz_dt, name=""):
    """
    Pure-Python seasonal adjustment fallback using STL. No external binary,
    so it can't hit X-13's automatic-model-selection failures. Used only when
    X-13 fails (e.g. degenerate regARIMA residuals on flat/interpolated
    series like corporation tax).
    """
    from statsmodels.tsa.seasonal import STL
    res = STL(nz_dt, period=4, robust=True).fit()
    # SA series = observed minus seasonal component
    return nz_dt - res.seasonal


def x13_sa(series, name=""):
    """
    Seasonally adjust a quarterly pd.Series. Tries X-13ARIMA-SEATS first;
    if X-13 fails for data reasons (not install reasons), falls back to STL.
    Returns a pd.Series indexed like the input. NaN head/tail preserved.
    Series with fewer than 16 non-NaN observations are returned unchanged.
    """
    from statsmodels.tsa.x13 import x13_arima_analysis, X13Error, X13NotFoundError

    s = series.copy().astype(float)
    nz = s.dropna()
    if len(nz) < 16:
        print(f"  [x13_sa] {name}: only {len(nz)} obs, skipping SA")
        return s

    # X-13 requires DatetimeIndex; convert PeriodIndex if needed
    if isinstance(nz.index, pd.PeriodIndex):
        nz_dt = nz.copy()
        nz_dt.index = nz_dt.index.to_timestamp(how="start")
        nz_dt.index.freq = "QS"
    else:
        nz_dt = nz

    sa = None
    # Attempt X-13 with a couple of spec variants.
    for i, kw in enumerate((
        dict(outlier=True,  trading=False),
        dict(outlier=False, trading=False),
    )):
        try:
            sa = x13_arima_analysis(nz_dt, freq="Q", **kw).seasadj
            if i > 0:
                print(f"  [x13_sa] {name}: X-13 (outlier off)")
            break
        except X13NotFoundError:
            # Install/path problem -- this is fatal, not a per-series data issue.
            raise RuntimeError(
                f"X-13 binary not found while adjusting {name}. Get it from "
                "https://www.census.gov/data/software/x13as.html and set X13PATH."
            )
        except (X13Error, Exception):
            continue

    # If X-13 failed for data reasons, fall back to STL.
    if sa is None:
        sa = _stl_sa(nz_dt, name)
        print(f"  [x13_sa] {name}: X-13 failed on this series; used STL fallback")

    # Convert back to PeriodIndex matching input
    if isinstance(series.index, pd.PeriodIndex):
        sa.index = sa.index.to_period("Q")

    out = pd.Series(np.nan, index=series.index, dtype=float)
    out.loc[sa.index.intersection(out.index)] = sa.loc[sa.index.intersection(out.index)]
    return out

START_Q = pd.Period("1955Q1", "Q")
END_Q   = pd.Period(pd.Timestamp.today(), "Q") - 1  # one quarter back; ONS lags a quarter


#%%
# =============================================================================
# ONS Time Series API helper
#
# URL format: https://www.ons.gov.uk/{topic_path}/timeseries/{cdid}/{dataset}/data
# Returns JSON with "quarters" list of {date, value, year, quarter} entries.
# We only need cdid + dataset_id; topic_path is in the JSON response anyway,
# but the request also works if we use a search-based fallback.

def ons_q(cdid, dataset_id="qna", topic_path=None, monthly_how="mean"):
    """
    Pull a quarterly ONS series by CDID. Returns pd.Series indexed by Period[Q].

    Args:
        cdid: 4-letter ONS code, e.g. "YBHA" (nominal GDP), "ABMI" (real GDP)
        dataset_id: dataset publication, usually "qna" (quarterly national accounts),
                    "ct" (consumer trends), "drsi" (gov receipts/expenditure),
                    "mret" (monthly retail), "lms" (labour market)
        topic_path: optional override; if None, tries common paths
    """
    candidate_paths = [topic_path] if topic_path else [
        "economy/nationalaccounts/satelliteaccounts",
        "economy/grossdomesticproductgdp",
        "economy/nationalaccounts/uksectoraccounts",
        "economy/governmentpublicsectorandtaxes/publicspending",
        "economy/governmentpublicsectorandtaxes/publicsectorfinance",
        "economy/inflationandpriceindices",
        "employmentandlabourmarket/peopleinwork/employmentandemployeetypes",
        "economy/nationalaccounts/balanceofpayments",
    ]
    for path in candidate_paths:
        url = f"https://www.ons.gov.uk/{path}/timeseries/{cdid.lower()}/{dataset_id}/data"
        for attempt in range(3):
            try:
                r = requests.get(url, timeout=120,
                                 headers={"Accept": "application/json"})
                if r.status_code == 404:
                    break  # try next path
                r.raise_for_status()
                payload = r.json()

                quarters = payload.get("quarters", [])
                if quarters:
                    rows = []
                    for q in quarters:
                        yr = int(q["year"])
                        qn = int(q["quarter"].lstrip("Q"))
                        val = pd.to_numeric(q.get("value", ""), errors="coerce")
                        rows.append((pd.Period(year=yr, quarter=qn, freq="Q"), val))
                    s = pd.Series({p: v for p, v in rows}).sort_index()
                    s.name = cdid.upper()
                    return s

                months = payload.get("months", [])
                if months:
                    rows = []
                    for m in months:
                        yr = int(m["year"])
                        month_num = pd.to_datetime(m["month"], format="%B").month
                        val = pd.to_numeric(m.get("value", ""), errors="coerce")
                        rows.append((pd.Timestamp(year=yr, month=month_num, day=1), val))
                    d = pd.Series({t: v for t, v in rows}).sort_index()
                    d.index = d.index.to_period("Q")
                    if monthly_how == "sum":
                        s = d.groupby(level=0).sum(min_count=1)
                    elif monthly_how == "mean":
                        s = d.groupby(level=0).mean()
                    else:
                        raise ValueError("monthly_how must be 'mean' or 'sum'")
                    s.name = cdid.upper()
                    return s

                years = payload.get("years", [])
                if years:
                    rows = []
                    for y in years:
                        yr = int(y["year"])
                        val = pd.to_numeric(y.get("value", ""), errors="coerce")
                        for qn in range(1, 5):
                            rows.append((pd.Period(year=yr, quarter=qn, freq="Q"), val))
                    s = pd.Series({p: v for p, v in rows}).sort_index()
                    s.name = cdid.upper()
                    return s

                break
            except (requests.ConnectionError, requests.Timeout):
                if attempt < 2:
                    time.sleep(2.0 ** attempt)
                else:
                    raise
    raise RuntimeError(f"Could not retrieve ONS series {cdid} from any candidate path")


def try_ons_q(cdid, dataset_ids, *, monthly_how="mean", topic_path=None):
    """Try the same CDID across multiple ONS datasets."""
    last_err = None
    for ds in dataset_ids:
        try:
            return ons_q(cdid, ds, topic_path=topic_path, monthly_how=monthly_how)
        except Exception as e:
            last_err = e
    raise RuntimeError(f"Could not retrieve ONS series {cdid} from datasets {dataset_ids}: {last_err}")


def ons_m(cdid, dataset_id, topic_path=None):
    """
    Pull a monthly ONS series and aggregate to quarterly mean.
    """
    candidate_paths = [topic_path] if topic_path else [
        "employmentandlabourmarket/peopleinwork/employmentandemployeetypes",
        "economy/inflationandpriceindices",
        "economy/grossdomesticproductgdp",
    ]
    for path in candidate_paths:
        url = f"https://www.ons.gov.uk/{path}/timeseries/{cdid.lower()}/{dataset_id}/data"
        for attempt in range(3):
            try:
                r = requests.get(url, timeout=120,
                                 headers={"Accept": "application/json"})
                if r.status_code == 404:
                    break
                r.raise_for_status()
                payload = r.json()
                months = payload.get("months", [])
                if not months:
                    break
                rows = []
                for m in months:
                    yr = int(m["year"])
                    mn = m["month"]  # text like "January"
                    month_num = pd.to_datetime(mn, format="%B").month
                    val = pd.to_numeric(m.get("value", ""), errors="coerce")
                    rows.append((pd.Timestamp(year=yr, month=month_num, day=1), val))
                d = pd.Series({t: v for t, v in rows}).sort_index()
                # Aggregate monthly to quarterly mean
                d.index = d.index.to_period("Q")
                return d.groupby(level=0).mean()
            except (requests.ConnectionError, requests.Timeout):
                if attempt < 2:
                    time.sleep(2.0 ** attempt)
                else:
                    raise
    raise RuntimeError(f"Could not retrieve ONS monthly series {cdid}")


def _load_api_keys(path):
    """Read KEY=value pairs if the file exists; return an empty dict otherwise."""
    keys = {}
    if not os.path.exists(path):
        return keys
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            keys[k.strip()] = v.strip().strip('"').strip("'")
    return keys


def fred_q(series_id):
    """Pull a FRED series through the official FRED API and convert to quarterly mean."""
    keys = _load_api_keys(os.path.join(RAW_DIR, "api_keys.txt"))
    if "FRED_KEY" not in keys:
        raise RuntimeError(
            f"FRED_KEY not found in {os.path.join(RAW_DIR, 'api_keys.txt')}. "
            "Use the same raw/api_keys.txt file as the US data script."
        )

    fred_key = keys["FRED_KEY"]
    url = "https://api.stlouisfed.org/fred/series/observations"
    params = {
        "series_id": series_id,
        "api_key": fred_key,
        "file_type": "json",
        "observation_start": "1900-01-01",
    }

    last_exc = None
    for attempt in range(4):
        try:
            r = requests.get(url, params=params, timeout=120)
            r.raise_for_status()
            obs = r.json().get("observations", [])
            if not obs:
                raise RuntimeError(f"FRED API returned no data for {series_id}")

            df = pd.DataFrame(obs)
            df["date"] = pd.to_datetime(df["date"], errors="coerce")
            df["value"] = pd.to_numeric(df["value"], errors="coerce")
            df = df.dropna(subset=["date", "value"])

            s = df.set_index("date")["value"]
            s.index = s.index.to_period("Q")
            s.name = series_id
            return s.groupby(level=0).mean().sort_index()

        except Exception as e:
            last_exc = e
            if attempt < 3:
                time.sleep(2.0 ** attempt)
            else:
                raise RuntimeError(f"FRED API failed for {series_id}: {last_exc}")

def build_bank_rate_series(start_q=START_Q, end_q=END_Q):
    """
    Bank Rate series:
    - FRED BOERUKQ through 2016Q4.
    - Official BoE Bank Rate history from 2017Q1 onward if available.
    - Prints an overlap check when both sources load.
    """
    fred_rate = pd.Series(dtype=float)
    boe_rate = pd.Series(dtype=float)

    try:
        fred_rate = fred_q("BOERUKQ")
        fred_rate = fred_rate.loc[start_q:min(end_q, pd.Period("2016Q4", "Q"))]
        if len(fred_rate) > 0:
            print(f"  FRED BOERUKQ: {fred_rate.index.min()} -> {fred_rate.index.max()}")
    except Exception as e:
        print(f"  WARNING: FRED BOERUKQ unavailable ({e})")

    try:
        boe_rate = boe_bank_rate_history_q(start_year=1975)
        boe_rate = boe_rate.loc[start_q:end_q]
        if len(boe_rate) > 0:
            print(f"  BoE official Bank Rate: {boe_rate.index.min()} -> {boe_rate.index.max()}")
    except Exception as e:
        print(f"  WARNING: BoE official Bank Rate unavailable ({e})")

    if len(fred_rate) > 0 and len(boe_rate) > 0:
        overlap = fred_rate.index.intersection(boe_rate.index)
        if len(overlap) > 0:
            diff = (fred_rate.loc[overlap] - boe_rate.loc[overlap]).abs().dropna()
            if len(diff) > 0:
                print(f"  FRED/BoE overlap check, max abs diff: {diff.max():.4f} pp")

    pieces = []
    if len(fred_rate) > 0:
        pieces.append(fred_rate.loc[:pd.Period("2016Q4", "Q")])
    if len(boe_rate) > 0:
        pieces.append(boe_rate.loc[pd.Period("2017Q1", "Q"):])

    if not pieces:
        raise RuntimeError("No Bank Rate source was available")

    out = pd.concat(pieces).sort_index()
    out = out[~out.index.duplicated(keep="last")]
    return out


#%%
# =============================================================================
# Bank of England Interactive Database helper
#
# URL format with date range:
#   https://www.bankofengland.co.uk/boeapps/database/_iadb-fromshowcolumns.asp?
#       Travel=NIxAZxSUx&FromSeries=1&ToSeries=50&DAT=RNG&FD=1&FM=Jan&FY=1955
#       &TD=31&TM=Dec&TY=2026&FNY=&CSVF=TT&html.x=66&html.y=26&SeriesCodes={code}
#       &UsingCodes=Y&Filter=N&title=Bank+Rate&VPD=Y
# Returns CSV with two columns: DATE (dd Mon yyyy), <SeriesCode>.

def boe_bank_rate_history_q(start_year=1975):
    """Pull official Bank Rate changes from BoE and convert to quarterly mean."""
    from io import StringIO

    url = "https://www.bankofengland.co.uk/boeapps/database/Bank-Rate.asp"
    headers = {
        "User-Agent": "Mozilla/5.0",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    }

    try:
        for attempt in range(3):
            try:
                r = requests.get(url, timeout=120, headers=headers)
                r.raise_for_status()
                tables = pd.read_html(StringIO(r.text))
                break
            except Exception:
                if attempt < 2:
                    time.sleep(2.0 ** attempt)
                else:
                    raise

        rate_table = None
        for t in tables:
            cols = [str(c).strip().lower() for c in t.columns]
            if any("date" in c for c in cols) and any("rate" in c for c in cols):
                rate_table = t.copy()
                break
        if rate_table is None:
            raise RuntimeError("Could not find Bank Rate table on BoE page")

        rate_table.columns = [str(c).strip() for c in rate_table.columns]
        date_col = [c for c in rate_table.columns if "date" in c.lower()][0]
        rate_col = [c for c in rate_table.columns if "rate" in c.lower()][0]

        changes = pd.DataFrame({
            "date": pd.to_datetime(rate_table[date_col].astype(str).str.strip(),
                                   format="%d %b %y", errors="coerce"),
            "rate": pd.to_numeric(rate_table[rate_col], errors="coerce"),
        }).dropna().sort_values("date")

        if changes.empty:
            raise RuntimeError("BoE Bank Rate table parsed but no usable rows found")

        start_date = pd.Timestamp(year=start_year, month=1, day=1)
        end_date = pd.Timestamp.today().normalize()
        daily_idx = pd.date_range(start_date, end_date, freq="D")
        daily = pd.Series(index=daily_idx, dtype=float)

        before_start = changes.loc[changes["date"] <= start_date]
        if not before_start.empty:
            daily.iloc[0] = before_start.iloc[-1]["rate"]

        for d, rr in changes.itertuples(index=False):
            if start_date <= d <= end_date:
                daily.loc[d] = rr

        daily = daily.ffill().dropna()
        daily.index = daily.index.to_period("Q")
        return daily.groupby(level=0).mean()

    except Exception as e:
        print(f"  BoE page blocked/unavailable ({e}); using FRED BOERUKQ fallback")
        s = fred_q("BOERUKQ")
        return s.loc[pd.Period(f"{start_year}Q1", "Q"):]


#%%
# =============================================================================
# helpers shared with the US pipeline (copied to keep file standalone)

def log_linear_detrend(x, scale=100.0):
    arr = np.asarray(x, dtype=float)
    log_arr = np.log(arr)
    t = np.arange(len(log_arr))
    X = np.column_stack([np.ones_like(t, dtype=float), t.astype(float)])
    beta = np.linalg.lstsq(X, log_arr, rcond=None)[0]
    return scale * (log_arr - X @ beta)

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
# Pull ONS series
#
# Mapping decisions (CVM series are used directly where available):
#   YBHA  Gross Domestic Product (current prices)            [GDP nominal]
#   ABMI  GDP chained volume measure                          [real GDP -- not used]
#   YBGB  GDP deflator (current price/chained volume × 100)   [pi]
#   ABJR+HAYO  Household + NPISH final consumption, CVM      -> chat
#   NPQT      Gross fixed capital formation, CVM              -> ihat
#   NMRY      General government final consumption, CVM        -> ghat
#   EBAQ      Quarterly interpolated population                -> per-capita scaler
#
# Hours worked: YBUS is total actual weekly hours worked, all workers,
# seasonally adjusted, available from the Labour Market Statistics dataset.

print("Pulling UK quarterly national accounts from ONS...")

# National accounts aggregates
gdp_nom  = ons_q("YBHA", "qna")   # GDP, current prices, £m, SA
gva_real = ons_q("ABMI", "qna")   # GDP, chained volume measure, £m, SA

# Consumption follows the Bank-of-England-style UK convention:
# households + NPISH, using CVM levels directly (do not deflate again).
hhfce_real = ons_q("ABJR", "qna")   # Household final consumption, CVM, £m, SA
npish_real = ons_q("HAYO", "qna")   # NPISH final consumption, CVM, £m, SA
hhfce_nom  = ons_q("ABJQ", "qna")   # Household final consumption, current prices, £m, SA
npish_nom  = ons_q("HAYE", "qna")   # NPISH final consumption, current prices, £m, SA
cons_real  = hhfce_real + npish_real
cons_nom   = hhfce_nom + npish_nom
cons_deflator = (cons_nom / cons_real) * 100.0

# These ONS QNA series are already chained-volume measures.
gfcf_real  = ons_q("NPQT", "qna")   # Gross fixed capital formation, CVM, £m, SA
ggfce_real = ons_q("NMRY", "qna")   # General government final consumption, CVM, £m, SA
ggfce_nom  = try_ons_q("NMRP", ["qna", "ukea", "pn2"])  # General government final consumption, current prices, £m, SA

# Population: EBAQ is already quarterly interpolated, in thousands.
pop_q = ons_q("EBAQ", "bb") * 1000.0

# Debt and transfers from Public Sector Finances.
# HF6W is PSND excluding public sector banks, £bn CPNSA. Convert £bn -> £m.
# Keep the source coverage as reported by ONS; do not manually delete early quarters.
debt_uk = ons_q("HF6W", "pusf") * 1000.0

# Social benefits: PUSF net social benefits, nominal £m.
soc_ben = try_ons_q("RPGG", ["ukea"])  # GG D.62 social benefits paid (excluding in-kind), CP SA, £m

# Hours worked: total actual weekly hours, all workers, seasonally adjusted.
try:
    hours = ons_q("YBUS", "lms")  # total actual weekly hours worked, millions, SA
except Exception:
    print("  WARNING: hours series YBUS not available; lhat will be skipped")
    hours = pd.Series(dtype=float)

# Tax-rate inputs. Revenue series are public-sector finance flows; use
# quarterly values when available and sum monthly values if ONS only returns
# months. Income bases are current-price national accounts series.
print("\nPulling UK tax inputs...")
EC      = try_ons_q("DTWM", ["qna", "pn2"])       # compensation of employees, CP SA £m
W       = try_ons_q("DTWL", ["qna", "pn2"])       # wages and salaries, CP SA £m
PRI     = try_ons_q("ROYH", ["ukea", "qna"])      # HH+NPISH gross mixed income, CP SA £m
CP      = try_ons_q("CGBZ", ["qna", "pn2"])       # corporations gross operating surplus, CP SA £m

# Household net interest (Jones 2002 NI component). UKEA dataset.
# I69P: Households interest receivable, £m, NSA
# HACY: Households interest payable,    £m, NSA
# NSA at this disaggregation level; the seasonality mismatch with the SA
# series above is small for postwar UK household net interest.
NI_recv = try_ons_q("I69P", ["ukea"], monthly_how="sum")
NI_paid = try_ons_q("HACY", ["ukea"], monthly_how="sum")
NI      = NI_recv - NI_paid

# Household imputed rental (Jones 2002 Rental component). UKEA dataset.
# HABM: Households (S.14) Operating Surplus, gross (B.2g) Resources, £m, NSA
# Per ONS, B.2g for the household sector "is comprised almost entirely of
# imputed rental, a national accounts concept which captures the value to
# owner-occupiers of living in their own home." This is the matching
# denominator term for the council tax (CTAX) revenue in the numerator.
RENTAL = try_ons_q("HABM", ["ukea"], monthly_how="sum")

CSI     = try_ons_q("AIIH", ["pusf"], monthly_how="sum")  # compulsory social contributions / NICs
try:
    income_tax_plus_nics = try_ons_q("KSS8", ["pusf"], monthly_how="sum")
    PIT = income_tax_plus_nics - CSI
except Exception:
    # Fallback if the combined series is unavailable: PAYE + other income tax
    # + self-assessed income tax.
    paye = try_ons_q("MS6W", ["pusf"], monthly_how="sum")
    other_income_tax = try_ons_q("MF6X", ["pusf"], monthly_how="sum")
    self_assessed = try_ons_q("LISB", ["pusf"], monthly_how="sum")
    PIT = paye + other_income_tax + self_assessed

CT      = try_ons_q("ACCD", ["pusf", "bb", "qna"], monthly_how="sum")  # corporation tax
BRATES  = try_ons_q("CUKY", ["pusf"], monthly_how="sum")  # business rates
CTAX    = try_ons_q("NMHM", ["pusf"], monthly_how="sum")  # council tax
PT      = BRATES + CTAX
# Prefer the quarterly national-accounts version of D.2 when available.
# NMYE is general government current receipts: taxes on production (D.2).
# It is the closest UK analog to the US NIPA "taxes on production and imports"
# aggregate used in the MRT/Jones consumption-tax formula.
TPI     = try_ons_q("NMYE", ["ukea", "pusf"], monthly_how="sum")

print("  Tax inputs pulled: DTWM, DTWL, ROYH, CGBZ, HABM, I69P-HACY, AIIH/PIT, ACCD, CUKY, NMHM, NMYE")

# Apply X-13ARIMA-SEATS seasonal adjustment to all NSA inputs that feed into
# tax-rate construction. The denominators (EC, W, PRI, CP) are already SA from
# the qna dataset and do NOT need adjustment. PUSF tax-receipt series and UKEA
# household property-income series are NSA-only, so they are SA'd here.
print("\nApplying X-13 seasonal adjustment to NSA tax-revenue and property-income inputs...")
NI_recv = x13_sa(NI_recv, "NI_recv (I69P)")
NI_paid = x13_sa(NI_paid, "NI_paid (HACY)")
NI      = NI_recv - NI_paid
RENTAL  = x13_sa(RENTAL,  "RENTAL (HABM)")
CSI     = x13_sa(CSI,     "CSI (AIIH)")
PIT     = x13_sa(PIT,     "PIT (KSS8 - CSI)")
CT      = x13_sa(CT,      "CT (ACCD)")
BRATES  = x13_sa(BRATES,  "BRATES (CUKY)")
CTAX    = x13_sa(CTAX,    "CTAX (NMHM)")
PT      = BRATES + CTAX
TPI     = x13_sa(TPI,     "TPI (NMYE)")

print("  GDP nominal (YBHA):", gdp_nom.index.min(), "->", gdp_nom.index.max(), f"({len(gdp_nom)} qtrs)")
print("  Consumption real: ", cons_real.index.min(), "->", cons_real.index.max())
print("  GFCF (NPQT):      ", gfcf_real.index.min(), "->", gfcf_real.index.max())
print("  GGFCE (NMRY):     ", ggfce_real.index.min(), "->", ggfce_real.index.max())
print("  Debt (HF6W):      ", debt_uk.dropna().index.min() if debt_uk.notna().any() else "N/A",
      "->", debt_uk.dropna().index.max() if debt_uk.notna().any() else "N/A")
print("  Hours (YBUS):     ", hours.index.min() if len(hours) else "N/A",
      "->", hours.index.max() if len(hours) else "N/A")


#%%
# =============================================================================
# Construct GDP deflator from nominal/real ratio
#
# pi[t] = deflator[t] / deflator[t-1]
# where deflator = (nominal GDP / real GDP) * 100 (so SS deflator ~ 100 in base year)

deflator = (gdp_nom / gva_real) * 100.0
deflator = deflator.dropna()
print(f"\nGDP deflator constructed: range {deflator.min():.2f} to {deflator.max():.2f}")


#%%
# =============================================================================
# Pull Bank Rate
#
# Use FRED BOERUKQ through 2016Q4, then splice official BoE Bank Rate
# from 2017Q1 onward. The overlap is checked against BoE official history.

print("\nPulling Bank Rate...")
try:
    bank_rate_q = build_bank_rate_series(START_Q, END_Q)
    print(f"  Bank Rate: {bank_rate_q.index.min()} -> {bank_rate_q.index.max()}")
except Exception as e:
    print(f"  ERROR pulling Bank Rate: {e}")
    bank_rate_q = pd.Series(dtype=float)

# Optional override/supplement file
pre75_path = os.path.join(RAW_DIR, "uk_bank_rate_pre1975.csv")
if os.path.exists(pre75_path):
    pre75 = pd.read_csv(pre75_path)
    pre75["QPER"] = pd.PeriodIndex(pre75["QPER"], freq="Q")
    pre75 = pre75.set_index("QPER")["bank_rate"]
    bank_rate_q = pd.concat([pre75, bank_rate_q]).sort_index()
    bank_rate_q = bank_rate_q[~bank_rate_q.index.duplicated(keep="last")]
    print(f"  Local Bank Rate supplement loaded: {len(pre75)} obs")


#%%
# =============================================================================
# Per-capita scaler
#
# EBAQ is already quarterly interpolated and was converted from thousands
# to persons above. No additional interpolation is needed.

print("\nUsing quarterly interpolated ONS population (EBAQ).")


#%%
# =============================================================================
# Build observable hats
#
# For CVM aggregates (chat, ihat, ghat), use the ONS real series directly.
# For nominal stock/flow series (bhat, trhat), deflate by the GDP deflator.
# Then divide by population and take 100 * log-linear-detrended deviation.

idx = pd.period_range(START_Q, END_Q, freq="Q")

def build_real_pc_hat_from_real(real, pop, name):
    """Build a real per-capita log-linear-detrended hat series from a CVM series."""
    common = real.index.intersection(pop.index)
    if len(common) == 0:
        return to_yq_frame(idx, **{name: np.full(len(idx), np.nan)})
    real_pc = (real.loc[common] / pop.loc[common]).dropna()
    if len(real_pc) < 2:
        return to_yq_frame(idx, **{name: np.full(len(idx), np.nan)})
    hat = log_linear_detrend(real_pc.values)
    out = pd.Series(np.nan, index=idx)
    out.loc[real_pc.index] = hat
    return to_yq_frame(idx, **{name: out.to_numpy()})


def build_real_pc_hat_from_nominal(nominal, deflator, pop, name):
    """Build a real per-capita log-linear-detrended hat series from a nominal series."""
    common = nominal.index.intersection(deflator.index).intersection(pop.index)
    if len(common) == 0:
        return to_yq_frame(idx, **{name: np.full(len(idx), np.nan)})
    real = nominal.loc[common] / (deflator.loc[common] / 100.0)
    return build_real_pc_hat_from_real(real, pop.loc[common], name)


df_c  = build_real_pc_hat_from_real(cons_real,  pop_q, "chat")
df_i  = build_real_pc_hat_from_real(gfcf_real,  pop_q, "ihat")
df_g  = build_real_pc_hat_from_real(ggfce_real, pop_q, "ghat")
df_b  = build_real_pc_hat_from_nominal(debt_uk, deflator, pop_q, "bhat")
df_tr = build_real_pc_hat_from_nominal(soc_ben, deflator, pop_q, "trhat")

# Hours: log-linear-detrended in levels (no deflator, no per-capita scaling
# beyond what's implicit in the index). The US lhat uses HOANBS (index, sa).
# UK YBUS is total weekly hours; we treat it the same way.
if len(hours) > 0:
    hpc = hours.dropna()
    hat = log_linear_detrend(hpc.values)
    out = pd.Series(np.nan, index=idx)
    common = idx.intersection(hpc.index)
    out.loc[common] = hat[:len(common)]
    df_l = to_yq_frame(idx, lhat=out.to_numpy())
else:
    df_l = to_yq_frame(idx, lhat=np.full(len(idx), np.nan))

# Inflation: gross deflator ratio
def build_pi(deflator, idx):
    s = deflator.dropna()
    pi = s / s.shift(1)
    out = pd.Series(np.nan, index=idx)
    common = idx.intersection(pi.index)
    out.loc[common] = pi.loc[common].values
    return to_yq_frame(idx, pi=out.to_numpy())

df_pi = build_pi(deflator, idx)

# Bank Rate: percent annualized -> quarterly decimal
def build_r(bank_rate, idx):
    out = pd.Series(np.nan, index=idx)
    common = idx.intersection(bank_rate.index)
    out.loc[common] = bank_rate.loc[common].values / 400.0
    return to_yq_frame(idx, r=out.to_numpy())

df_r = build_r(bank_rate_q, idx)


#%%
# =============================================================================
# Build UK average tax rates and hat variants
#
# UK counterpart to the US Jones/Mendoza-Razin-Tesar tax construction.
# The UK capital-income base is necessarily a proxy because the compact ONS
# time-series API does not expose all Jones components in the same way as BEA
# NIPA; use corporate gross operating surplus plus half of mixed income.

def clean_rate(s):
    s = s.replace([np.inf, -np.inf], np.nan)
    return s.where((s >= 0.0) & (s <= 1.0))


def build_tax_hat_frame(rate, name):
    level = pd.Series(np.nan, index=idx)
    common = idx.intersection(rate.index)
    level.loc[common] = rate.loc[common].values
    level = clean_rate(level)

    obs = level.dropna()
    detrended = pd.Series(np.nan, index=idx)
    if len(obs) >= 2:
        detrended.loc[obs.index] = linear_detrend_level(obs.values)
        mu = obs.mean()
        hat_dt = (detrended - mu) / mu
        hat_nd = (level - mu) / mu
    else:
        mu = np.nan
        hat_dt = pd.Series(np.nan, index=idx)
        hat_nd = pd.Series(np.nan, index=idx)

    return to_yq_frame(
        idx,
        **{
            name: level.to_numpy(),
            f"{name}_detrend": detrended.to_numpy(),
            f"{name}_hat_detrend": hat_dt.to_numpy(),
            f"{name}_hat_nodetrend": hat_nd.to_numpy(),
        }
    )


CI = RENTAL + CP + NI + PRI / 2.0
tau_p_uk = PIT / (W + PRI / 2.0 + CI)

tau_l_uk = (tau_p_uk * (W + PRI / 2.0) + CSI) / (EC + PRI / 2.0)
tau_k_uk = (tau_p_uk * CI + CT + PT) / CI

# Consumption tax: UK analog of the US MRT/Jones formula.
# US benchmark in 1b_build_data.py:
#     tau_c = (taxes on production/imports - property taxes)
#             / (private consumption + government consumption - numerator)
#
# For the UK, NMYE is the broad D.2 taxes-on-production aggregate.
# Business rates (CUKY) are property-tax-like taxes included in D.2, so they
# are removed from the consumption-tax numerator. Council tax is not subtracted
# here because it is not part of NMYE's production-tax aggregate; it is kept in
# PT for the capital-tax formula.
property_tax_like_in_D2 = BRATES
indirect_tax = TPI - property_tax_like_in_D2
tau_c_uk = indirect_tax / (cons_nom + ggfce_nom - indirect_tax)

df_tl = build_tax_hat_frame(tau_l_uk, "tau_l")
df_tk = build_tax_hat_frame(tau_k_uk, "tau_k")
df_tc = build_tax_hat_frame(tau_c_uk, "tau_c")

print("\nUK tax-rate summary:")
for nm, rr in [("tau_l", tau_l_uk), ("tau_k", tau_k_uk), ("tau_c", tau_c_uk)]:
    s = clean_rate(rr).dropna()
    if len(s):
        print(f"  {nm:6s}  {s.index.min()} -> {s.index.max()}  mean={s.mean():.4f}  obs={len(s)}")
    else:
        print(f"  {nm:6s}  no usable observations")


#%%
# =============================================================================
# Merge and export

frames = [df_c, df_i, df_g, df_b, df_pi, df_r, df_l, df_tr, df_tl, df_tk, df_tc]
df_out = frames[0]
for f in frames[1:]:
    df_out = df_out.merge(f, on=["year", "quarter"], how="outer")
df_out = df_out.sort_values(["year", "quarter"]).reset_index(drop=True)

out_path = os.path.join(OUT_DIR, "UK_Data_Full.xlsx")
df_out.to_excel(out_path, index=False)

print(f"\nWrote {out_path}")
print(f"  rows: {len(df_out)}  cols: {df_out.shape[1]}")
print(f"  range: {df_out['year'].min()}Q{df_out['quarter'].iloc[0]} -> "
      f"{df_out['year'].max()}Q{df_out['quarter'].iloc[-1]}")
print(f"  columns: {df_out.columns.tolist()}")

print("\nNon-NaN count by series:")
for col in df_out.columns:
    if col in ("year", "quarter"):
        continue
    n = df_out[col].notna().sum()
    print(f"  {col:8s}  {n:>4d} obs")
