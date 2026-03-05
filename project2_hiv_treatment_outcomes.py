"""
HIV Treatment Outcomes Dashboard — OCS Viral Suppression & ART Adherence
Author: Nicholas Steven
Target Role: Analyst, Health Data — OHTN
Repo: github.com/nicholasstevenr/OHTN-health-data-project

Computes viral suppression rates (stratified), time-to-suppression (KM),
ART adherence proxy (MPR), PHQ-9/suppression logistic regression,
and site-level funnel plot benchmarking.
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")

SUPPRESSION_THRESHOLD  = 200    # copies/mL
CONFIRMED_SUPPRESSION  = 2      # consecutive measurements below threshold
MPR_ADHERENCE_THRESHOLD = 0.80  # <80% = non-adherent
PHQ9_HIGH_BURDEN       = 15     # PHQ-9 ≥ 15 = severe depression


# ── Load ──────────────────────────────────────────────────────────────────────

def load(vl_path: str, phq_path: str, enrolment_path: str,
          pharmacy_path: str):
    vl   = pd.read_csv(vl_path,         parse_dates=["viral_load_date","art_start_date"])
    phq  = pd.read_csv(phq_path,        parse_dates=["phq9_date"])
    enrl = pd.read_csv(enrolment_path,  parse_dates=["enrolment_date","last_contact_date"])
    rx   = pd.read_csv(pharmacy_path,   parse_dates=["dispense_date","next_dispense_date"])
    print(f"VL records: {len(vl):,}  |  PHQ-9: {len(phq):,}  |  Rx: {len(rx):,}")
    return vl, phq, enrl, rx


# ── 1. Viral Suppression Rates ────────────────────────────────────────────────

def viral_suppression(vl: pd.DataFrame, enrl: pd.DataFrame) -> pd.DataFrame:
    """Most-recent VL per participant; flag suppressed if < 200 copies/mL."""
    most_recent = (
        vl.sort_values("viral_load_date")
        .groupby("participant_id")
        .last()
        .reset_index()
    )
    most_recent["suppressed"] = most_recent["viral_load_copies"] < SUPPRESSION_THRESHOLD
    merged = most_recent.merge(
        enrl[["participant_id","site_code","age_group","housing_status"]],
        on="participant_id", how="left"
    )

    overall = merged["suppressed"].mean() * 100
    print(f"\n── Viral Suppression ──")
    print(f"  Overall: {overall:.1f}%")

    strats = []
    for col in ["site_code", "age_group", "housing_status"]:
        if col in merged.columns:
            s = (merged.groupby(col)["suppressed"]
                 .agg(n="count", suppressed_n="sum")
                 .assign(suppression_pct=lambda x: (x["suppressed_n"]/x["n"]*100).round(1))
                 .reset_index()
                 .rename(columns={col: "stratum_value"}))
            s["stratifier"] = col
            strats.append(s)

    return merged, pd.concat(strats, ignore_index=True)


# ── 2. Time-to-Suppression (Kaplan-Meier) ────────────────────────────────────

def time_to_suppression(vl: pd.DataFrame) -> pd.DataFrame:
    """
    For ART-initiators: days from art_start_date to first confirmed suppression.
    Returns per-participant TTS and regimen-level summary.
    """
    vl_sorted = vl.sort_values(["participant_id", "viral_load_date"])
    vl_sorted["suppressed"] = vl_sorted["viral_load_copies"] < SUPPRESSION_THRESHOLD
    vl_sorted["prev_supp"] = vl_sorted.groupby("participant_id")["suppressed"].shift(1)
    # First date with 2 consecutive suppressed measurements
    confirmed = vl_sorted[vl_sorted["suppressed"] & vl_sorted["prev_supp"].fillna(False)]
    first_confirmed = (confirmed.groupby("participant_id")["viral_load_date"]
                       .min().reset_index()
                       .rename(columns={"viral_load_date": "first_suppression_date"}))

    art_init = vl_sorted[vl_sorted["art_start_date"].notna()][
        ["participant_id","art_start_date","art_regimen_class"]].drop_duplicates("participant_id")

    tts = art_init.merge(first_confirmed, on="participant_id", how="left")
    tts["days_to_suppression"] = (
        tts["first_suppression_date"] - tts["art_start_date"]
    ).dt.days
    tts["achieved_suppression"] = tts["days_to_suppression"].notna()

    # Regimen-level summary
    reg_summary = (
        tts.groupby("art_regimen_class")
        .agg(
            n                    = ("participant_id", "count"),
            n_suppressed         = ("achieved_suppression", "sum"),
            median_days_to_supp  = ("days_to_suppression", "median"),
            p25_days             = ("days_to_suppression", lambda x: x.quantile(0.25)),
            p75_days             = ("days_to_suppression", lambda x: x.quantile(0.75)),
        )
        .round(1)
        .reset_index()
    )
    print(f"\n── Time-to-Suppression by Regimen ──")
    print(reg_summary.to_string(index=False))
    return tts, reg_summary


# ── 3. Medication Possession Ratio (Adherence Proxy) ─────────────────────────

def medication_possession_ratio(rx: pd.DataFrame) -> pd.DataFrame:
    """
    MPR = days_supply_dispensed / days_in_period for rolling 6-month windows.
    """
    rx = rx.sort_values(["participant_id", "dispense_date"])
    rx["days_gap"] = (rx["next_dispense_date"] - rx["dispense_date"]).dt.days
    rx["days_supply"] = rx["days_supply"].fillna(30)   # default 30-day supply

    mpr = (
        rx.groupby("participant_id")
        .apply(lambda g: pd.Series({
            "total_days_supply": g["days_supply"].sum(),
            "total_follow_up_days": (g["next_dispense_date"].max() - g["dispense_date"].min()).days + 1,
        }))
        .reset_index()
    )
    mpr["mpr"] = (mpr["total_days_supply"] / mpr["total_follow_up_days"].replace(0, np.nan)).clip(0, 1).round(3)
    mpr["non_adherent"] = mpr["mpr"] < MPR_ADHERENCE_THRESHOLD

    print(f"\n── Medication Possession Ratio ──")
    print(f"  Median MPR:         {mpr['mpr'].median():.3f}")
    print(f"  Non-adherent (<0.8): {mpr['non_adherent'].mean()*100:.1f}%")
    return mpr


# ── 4. PHQ-9 & Suppression Logistic Regression ───────────────────────────────

def phq9_suppression_model(vl_merged: pd.DataFrame,
                             phq: pd.DataFrame) -> dict:
    """
    Logistic regression: viral suppression ~ PHQ-9 quartile + age + ART regimen
    """
    # Most recent PHQ-9 per participant
    recent_phq = (phq.sort_values("phq9_date")
                  .groupby("participant_id")["phq9_total"]
                  .last().reset_index())
    recent_phq["phq9_quartile"] = pd.qcut(recent_phq["phq9_total"], q=4, labels=[1,2,3,4]).astype(float)
    recent_phq["high_depression"] = (recent_phq["phq9_total"] >= PHQ9_HIGH_BURDEN).astype(int)

    model_df = vl_merged.merge(recent_phq[["participant_id","phq9_total","high_depression"]],
                                on="participant_id", how="inner").dropna(
        subset=["suppressed","phq9_total","high_depression"])

    if len(model_df) < 50:
        return {"error": "Insufficient data"}

    features = ["high_depression"]
    if "age_numeric" in model_df.columns:
        features.append("age_numeric")

    X = StandardScaler().fit_transform(model_df[features].values)
    y = model_df["suppressed"].astype(int).values

    lr = LogisticRegression(max_iter=500, random_state=42)
    lr.fit(X, y)

    coef_df = pd.DataFrame({
        "feature":     features,
        "coefficient": lr.coef_[0].round(3),
        "odds_ratio":  np.exp(lr.coef_[0]).round(3),
    })
    print(f"\n── PHQ-9 / Suppression Model ──")
    print(coef_df.to_string(index=False))
    return {"coefficients": coef_df}


# ── 5. Site Funnel Plot ───────────────────────────────────────────────────────

def site_funnel_plot(vl_merged: pd.DataFrame) -> pd.DataFrame:
    """
    Funnel plot: site suppression rate vs. expected based on case mix.
    Flags sites >2 SD from expected.
    """
    if "site_code" not in vl_merged.columns:
        return pd.DataFrame()

    overall_rate = vl_merged["suppressed"].mean()
    site_stats = (
        vl_merged.groupby("site_code")
        .agg(n=("suppressed","count"), observed_supp=("suppressed","sum"))
        .reset_index()
    )
    site_stats["observed_rate"] = site_stats["observed_supp"] / site_stats["n"]
    # Expected = overall_rate; control limits = ±2 SE
    site_stats["expected_n_supp"] = site_stats["n"] * overall_rate
    site_stats["se"] = np.sqrt(overall_rate * (1 - overall_rate) / site_stats["n"])
    site_stats["upper_2sd"] = overall_rate + 2 * site_stats["se"]
    site_stats["lower_2sd"] = overall_rate - 2 * site_stats["se"]
    site_stats["outlier_low"]  = site_stats["observed_rate"] < site_stats["lower_2sd"]
    site_stats["outlier_high"] = site_stats["observed_rate"] > site_stats["upper_2sd"]
    site_stats = site_stats.round(4)

    outliers = site_stats[site_stats["outlier_low"] | site_stats["outlier_high"]]
    print(f"\n── Site Funnel Plot ──")
    print(f"  Overall suppression rate: {overall_rate*100:.1f}%")
    print(f"  Outlier sites: {len(outliers)}")
    return site_stats


# ── Export ────────────────────────────────────────────────────────────────────

def export_all(results: dict, outdir: str = "output") -> None:
    import os; os.makedirs(outdir, exist_ok=True)
    for name, obj in results.items():
        if isinstance(obj, pd.DataFrame) and len(obj):
            path = f"{outdir}/{name}.csv"
            obj.to_csv(path, index=False)
            print(f"  Exported → {path}")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    vl, phq, enrl, rx = load(
        "data/ocs_viral_loads_synthetic.csv",
        "data/ocs_phq9_synthetic.csv",
        "data/ocs_enrolment_synthetic.csv",
        "data/ocs_pharmacy_synthetic.csv",
    )

    vl_merged, supp_strats = viral_suppression(vl, enrl)
    tts_df, tts_summary    = time_to_suppression(vl)
    mpr_df                 = medication_possession_ratio(rx)
    phq_model              = phq9_suppression_model(vl_merged, phq)
    funnel_df              = site_funnel_plot(vl_merged)

    export_all({
        "suppression_stratified": supp_strats,
        "time_to_suppression":    tts_df,
        "tts_by_regimen":         tts_summary,
        "medication_mpr":         mpr_df,
        "site_funnel_plot":       funnel_df,
    })
