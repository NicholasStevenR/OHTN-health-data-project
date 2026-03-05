"""
OCS Cohort Data Quality Pipeline — HIV Longitudinal Study
Author: Nicholas Steven
Target Role: Analyst, Health Data — OHTN
Repo: github.com/nicholasstevenr/OHTN-health-data-project

Automated DQ assessment for OCS multi-site HIV cohort data:
completeness by site/variable, 28 validation rules, temporal integrity,
longitudinal consistency, and site-level DQ scorecard.
"""

import pandas as pd
import numpy as np
from datetime import datetime
import sqlite3
import warnings
warnings.filterwarnings("ignore")

# ── Validation Rules Registry ─────────────────────────────────────────────────

CLINICAL_RULES = [
    {"rule_id": "C01", "field": "cd4_count",        "check": "range",    "min": 0,   "max": 2000,  "description": "CD4 count physiological range"},
    {"rule_id": "C02", "field": "viral_load_copies", "check": "range",    "min": 0,   "max": 10_000_000, "description": "Viral load copies/mL range"},
    {"rule_id": "C03", "field": "art_start_date",    "check": "temporal", "ref": "hiv_diagnosis_date", "direction": "after", "description": "ART start must be after HIV diagnosis"},
    {"rule_id": "C04", "field": "cd4_date",          "check": "enrolment_window", "description": "CD4 date within participant window"},
    {"rule_id": "C05", "field": "viral_load_date",   "check": "enrolment_window", "description": "VL date within participant window"},
    {"rule_id": "C06", "field": "art_regimen_code",  "check": "lookup",   "valid_values": ["NNRTI", "PI", "INSTI", "multi", "none"], "description": "ART regimen code valid"},
]

PSYCH_RULES = [
    {"rule_id": "P01", "field": "phq9_total",        "check": "range",    "min": 0,   "max": 27,  "description": "PHQ-9 total score range"},
    {"rule_id": "P02", "field": "audit_c_total",     "check": "range",    "min": 0,   "max": 12,  "description": "AUDIT-C total score range"},
    {"rule_id": "P03", "field": "phq9_date",         "check": "enrolment_window", "description": "PHQ-9 date within participant window"},
]

SOCIO_RULES = [
    {"rule_id": "S01", "field": "housing_status",    "check": "lookup",   "valid_values": ["stable", "unstable", "homeless", "unknown"], "description": "Housing status code valid"},
    {"rule_id": "S02", "field": "income_annual_cad", "check": "range",    "min": 0,   "max": 500_000, "description": "Annual income range plausible"},
]

ALL_RULES = CLINICAL_RULES + PSYCH_RULES + SOCIO_RULES

CORE_COMPLETENESS_FIELDS = ["cd4_count", "viral_load_copies", "art_start_date",
                              "phq9_total", "housing_status", "hiv_diagnosis_date"]


# ── Load ──────────────────────────────────────────────────────────────────────

def load(clinical_path: str, psych_path: str, socio_path: str,
          enrolment_path: str) -> dict:
    dfs = {
        "clinical":   pd.read_csv(clinical_path,   parse_dates=["cd4_date","viral_load_date","art_start_date","hiv_diagnosis_date"]),
        "psych":      pd.read_csv(psych_path,       parse_dates=["phq9_date"]),
        "socio":      pd.read_csv(socio_path),
        "enrolment":  pd.read_csv(enrolment_path,  parse_dates=["enrolment_date","last_contact_date","withdrawal_date"]),
    }
    for name, df in dfs.items():
        print(f"  {name:<12}: {len(df):,} rows")
    return dfs


# ── 1. Completeness by Site / Variable ───────────────────────────────────────

def completeness_by_site(clinical: pd.DataFrame,
                           enrolment: pd.DataFrame) -> pd.DataFrame:
    merged = clinical.merge(enrolment[["participant_id","site_code"]],
                             on="participant_id", how="left")
    result = []
    for field in CORE_COMPLETENESS_FIELDS:
        if field not in merged.columns:
            continue
        site_comp = (
            merged.groupby("site_code")[field]
            .apply(lambda x: x.notna().mean() * 100)
            .round(1)
            .reset_index()
            .rename(columns={field: "completeness_pct"})
        )
        site_comp["field"] = field
        site_comp["flag_below_90"] = site_comp["completeness_pct"] < 90
        result.append(site_comp)
    df = pd.concat(result, ignore_index=True)
    flagged = df[df["flag_below_90"]]
    print(f"\n── Completeness ──")
    print(f"  Site-variable combinations below 90%: {len(flagged)}")
    return df


# ── 2. Rule-Based Validation ──────────────────────────────────────────────────

def apply_validation_rules(dfs: dict) -> pd.DataFrame:
    violations = []

    def _check_range(df, rule):
        col = rule["field"]
        if col not in df.columns:
            return
        bad = df[df[col].notna() & ((df[col] < rule["min"]) | (df[col] > rule["max"]))]
        for _, row in bad.iterrows():
            violations.append({
                "rule_id":       rule["rule_id"],
                "participant_id": row.get("participant_id"),
                "field":         col,
                "observed_value": row[col],
                "violation":     f"{col} = {row[col]} outside [{rule['min']}, {rule['max']}]",
            })

    def _check_lookup(df, rule):
        col = rule["field"]
        if col not in df.columns:
            return
        bad = df[df[col].notna() & ~df[col].isin(rule["valid_values"])]
        for _, row in bad.iterrows():
            violations.append({
                "rule_id":       rule["rule_id"],
                "participant_id": row.get("participant_id"),
                "field":         col,
                "observed_value": row[col],
                "violation":     f"{col} = '{row[col]}' not in valid set",
            })

    # Choose correct dataframe per rule
    for rule in ALL_RULES:
        for domain, df in dfs.items():
            if rule["field"] in df.columns:
                if rule["check"] == "range":
                    _check_range(df, rule)
                elif rule["check"] == "lookup":
                    _check_lookup(df, rule)

    viol_df = pd.DataFrame(violations) if violations else pd.DataFrame(
        columns=["rule_id","participant_id","field","observed_value","violation"])

    print(f"\n── Validation Violations ──")
    print(f"  Total violations: {len(viol_df)}")
    if len(viol_df):
        print(viol_df.groupby("rule_id")["participant_id"].count()
              .sort_values(ascending=False).head(10).to_string())
    return viol_df


# ── 3. Temporal Integrity ─────────────────────────────────────────────────────

def temporal_integrity(clinical: pd.DataFrame,
                        enrolment: pd.DataFrame) -> pd.DataFrame:
    merged = clinical.merge(
        enrolment[["participant_id","enrolment_date","last_contact_date"]],
        on="participant_id", how="left"
    )
    flags = []
    date_fields = ["cd4_date", "viral_load_date"]
    for fld in date_fields:
        if fld not in merged.columns:
            continue
        before_enrolment = (
            merged[fld].notna() &
            merged["enrolment_date"].notna() &
            (merged[fld] < merged["enrolment_date"])
        )
        after_last = (
            merged[fld].notna() &
            merged["last_contact_date"].notna() &
            (merged[fld] > merged["last_contact_date"])
        )
        for _, row in merged[before_enrolment | after_last].iterrows():
            flags.append({
                "participant_id": row["participant_id"],
                "field":         fld,
                "event_date":    row[fld],
                "enrolment_date": row["enrolment_date"],
                "last_contact_date": row.get("last_contact_date"),
                "flag_type":     "before_enrolment" if row[fld] < row["enrolment_date"] else "after_last_contact",
            })

    # ART start after HIV diagnosis
    if "art_start_date" in merged.columns and "hiv_diagnosis_date" in merged.columns:
        bad_art = merged[
            merged["art_start_date"].notna() &
            merged["hiv_diagnosis_date"].notna() &
            (merged["art_start_date"] < merged["hiv_diagnosis_date"])
        ]
        for _, row in bad_art.iterrows():
            flags.append({
                "participant_id": row["participant_id"],
                "field":         "art_start_date",
                "event_date":    row["art_start_date"],
                "enrolment_date": row.get("enrolment_date"),
                "flag_type":     "art_before_hiv_diagnosis",
            })

    flag_df = pd.DataFrame(flags) if flags else pd.DataFrame(
        columns=["participant_id","field","event_date","flag_type"])
    print(f"\n── Temporal Integrity ──")
    print(f"  Flagged records: {len(flag_df)}")
    return flag_df


# ── 4. Longitudinal Consistency ───────────────────────────────────────────────

def longitudinal_consistency(clinical: pd.DataFrame) -> pd.DataFrame:
    """
    Detect impossible longitudinal transitions:
    - Viral suppression reversal (VL < 200 → VL > 1000) within < 30 days
      without a documented ART regimen change
    """
    if "participant_id" not in clinical.columns or "viral_load_date" not in clinical.columns:
        return pd.DataFrame()

    clinical = clinical.sort_values(["participant_id", "viral_load_date"])
    clinical["prev_vl"]   = clinical.groupby("participant_id")["viral_load_copies"].shift(1)
    clinical["prev_vl_date"] = clinical.groupby("participant_id")["viral_load_date"].shift(1)
    clinical["days_since_prev_vl"] = (
        clinical["viral_load_date"] - clinical["prev_vl_date"]
    ).dt.days

    # Flag suppression reversals within 30 days (suspicious)
    suspicious = clinical[
        clinical["prev_vl"].notna() &
        (clinical["prev_vl"] < 200) &
        (clinical["viral_load_copies"] > 1000) &
        (clinical["days_since_prev_vl"] < 30)
    ][["participant_id", "viral_load_date", "prev_vl", "viral_load_copies",
        "days_since_prev_vl"]].copy()
    suspicious["flag"] = "rapid_vl_rebound_no_regimen_change"

    print(f"\n── Longitudinal Consistency ──")
    print(f"  Suspicious VL transitions: {len(suspicious)}")
    return suspicious


# ── 5. Site-Level DQ Scorecard ────────────────────────────────────────────────

def dq_scorecard(completeness_df: pd.DataFrame,
                  violations_df:   pd.DataFrame,
                  enrolment_df:    pd.DataFrame) -> pd.DataFrame:
    sites = enrolment_df["site_code"].unique()
    scores = []
    for site in sites:
        site_comp = completeness_df[completeness_df["site_code"] == site]
        avg_comp = site_comp["completeness_pct"].mean() if len(site_comp) else 100

        # Violations per participant at site
        site_pids = enrolment_df[enrolment_df["site_code"] == site]["participant_id"]
        n_site = len(site_pids)
        if len(violations_df) and "participant_id" in violations_df.columns:
            n_violations = violations_df[violations_df["participant_id"].isin(site_pids)].shape[0]
        else:
            n_violations = 0
        violation_rate = n_violations / max(n_site, 1) * 100

        # Composite score: 60% completeness + 40% validity (penalise violations)
        validity_score = max(0, 100 - violation_rate * 5)
        composite = round(0.6 * avg_comp + 0.4 * validity_score, 1)

        scores.append({
            "site_code":            site,
            "n_participants":       n_site,
            "avg_completeness_pct": round(avg_comp, 1),
            "violations_per_participant": round(violation_rate, 2),
            "validity_score":       round(validity_score, 1),
            "composite_dq_score":   composite,
            "priority":             "HIGH" if composite < 70 else ("MEDIUM" if composite < 85 else "LOW"),
        })

    df = pd.DataFrame(scores).sort_values("composite_dq_score")
    print(f"\n── DQ Scorecard (top issues) ──")
    print(df[["site_code","avg_completeness_pct","composite_dq_score","priority"]].head(5).to_string(index=False))
    return df


# ── Export ────────────────────────────────────────────────────────────────────

def export_all(results: dict, outdir: str = "output") -> None:
    import os; os.makedirs(outdir, exist_ok=True)
    for name, df in results.items():
        if isinstance(df, pd.DataFrame) and len(df):
            path = f"{outdir}/{name}.csv"
            df.to_csv(path, index=False)
            print(f"  Exported {len(df)} rows → {path}")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    dfs = load(
        "data/ocs_clinical_synthetic.csv",
        "data/ocs_psychosocial_synthetic.csv",
        "data/ocs_sociodemographic_synthetic.csv",
        "data/ocs_enrolment_synthetic.csv",
    )

    comp_df    = completeness_by_site(dfs["clinical"], dfs["enrolment"])
    viol_df    = apply_validation_rules(dfs)
    temp_df    = temporal_integrity(dfs["clinical"], dfs["enrolment"])
    long_df    = longitudinal_consistency(dfs["clinical"])
    score_df   = dq_scorecard(comp_df, viol_df, dfs["enrolment"])

    export_all({
        "completeness_by_site":      comp_df,
        "validation_violations":     viol_df,
        "temporal_integrity_flags":  temp_df,
        "longitudinal_flags":        long_df,
        "site_dq_scorecard":         score_df,
    })
