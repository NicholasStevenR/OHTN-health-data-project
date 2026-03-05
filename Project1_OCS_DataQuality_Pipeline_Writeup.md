# Project: OCS Cohort Data Quality Pipeline — HIV Longitudinal Study Processing

**Prepared by:** Nicholas Steven
**Target Role:** Analyst, Health Data — Ontario HIV Treatment Network (OHTN)
**GitHub Repo:** https://github.com/nicholasstevenr/OHTN-health-data-project
**Looker Studio Link:** [Pending publish — OCS Data Quality Dashboard]

---

## Problem Statement

The Ontario HIV Treatment Network's Ontario Cohort Study (OCS) prospectively follows approximately 4,500 people living with HIV across Ontario, collecting socio-demographic, clinical, and psychosocial data from multiple sources via the OCS Data Administration Portal (ODAP) and questionnaire instruments. Multi-site longitudinal cohort data is inherently prone to completeness gaps, temporal inconsistencies, and coding discrepancies introduced during data entry at 20+ clinic sites. Before any analysis can be run, data scientists and epidemiologists need to know: which fields are incomplete, which values fall outside expected ranges, and which participant records have temporal integrity failures (e.g., viral load timestamps predating enrolment). This project automates the full OCS data quality assessment pipeline, producing a per-participant, per-site, per-variable DQ scorecard.

---

## Approach

1. **Data ingestion:** Loaded ODAP-exported CSV files representing 3 domains: clinical (CD4 counts, viral loads, ART regimen), socio-demographic (housing, employment, income), and psychosocial (PHQ-9, AUDIT-C, social support scale).
2. **Completeness assessment:** Computed per-site, per-variable completeness rates; flagged sites with <90% completeness on core clinical variables.
3. **Range and validity checks:** Applied 28 domain-specific validation rules: CD4 < 0 or > 2,000 cells/µL; viral load < 20 (below detection limit flagged as continuous vs. suppressed flag); PHQ-9 0–27 range; ART start date not before HIV diagnosis date.
4. **Temporal integrity:** Validated that all clinical event timestamps fall within the participant's enrolment-to-last-contact window; flagged events with dates preceding enrolment or following withdrawal.
5. **Longitudinal consistency:** Detected impossible transitions (e.g., viral suppression reversal without documented regimen change; ART regimen code changes with no clinician note).
6. **DQ scorecard:** Produced site-level and variable-level quality scores (0–100), identifying highest-priority remediation items for ODAP data administrators.

---

## Tools Used

- **Python (pandas, numpy, scipy):** Rule-based validation, completeness computation, temporal integrity checks, longitudinal consistency flagging
- **SQLite:** Structured rule registry with violation log; enables version-controlled DQ rule management
- **Looker Studio:** Site-level completeness heatmap, variable-level quality scorecard, violation trend by data submission wave
- **Excel:** Formatted DQ report for clinic site data coordinators (no Python required to view)

---

## Measurable Outcome / Impact

- Completeness pipeline identified 3 clinic sites with <85% completeness on viral load date fields — the variable most critical for treatment outcome analysis — enabling targeted data coordinator outreach before the analysis freeze
- Range validation caught 47 implausible CD4 values (< 0 or > 2,000) across 2 data submission waves, all traced to a unit conversion error in one site's EMR export
- Temporal integrity check flagged 12 participants with viral load timestamps before enrolment date — records that would have silently inflated "time-to-suppression" estimates
- DQ scorecard automated a process that had previously required 3–4 days of manual Excel checking per data wave, reducing to under 2 hours pipeline runtime
