# Project: HIV Treatment Outcomes Dashboard — OCS Viral Suppression & ART Adherence Analytics

**Prepared by:** Nicholas Steven
**Target Role:** Analyst, Health Data — Ontario HIV Treatment Network (OHTN)
**GitHub Repo:** https://github.com/nicholasstevenr/OHTN-health-data-project
**Looker Studio Link:** [Pending publish — OCS Treatment Outcomes Dashboard]

---

## Problem Statement

OHTN researchers and program staff need to monitor HIV treatment outcomes across the OCS cohort — tracking viral suppression rates, time-to-suppression, ART adherence proxies, and mental health co-morbidity burden over time and across Ontario clinic sites. Policy stakeholders need to know: Are suppression rates improving? Are certain subgroups (by age, housing status, mental health burden) experiencing worse outcomes? Are clinic sites performing significantly above or below the provincial average? This project built an outcomes analytics dashboard answering these questions for the OCS longitudinal cohort.

---

## Approach

1. **Viral suppression analysis:** Defined viral suppression as VL < 200 copies/mL at the most recent measurement; computed suppression rates at cohort level and stratified by clinic site, age group (<30, 30–50, >50), diagnosis year, and housing status.
2. **Time-to-suppression:** For participants initiating ART within the OCS follow-up window, computed days from ART initiation to first confirmed suppression (VL < 200 on two consecutive measurements); Kaplan-Meier survival curves by ART regimen class (INSTI vs. NNRTI vs. PI).
3. **ART adherence proxy (prescription refill gap):** Computed medication possession ratio (MPR) from pharmacy refill date intervals as an adherence proxy; flagged participants with MPR < 0.8 over any 6-month window.
4. **Mental health co-morbidity burden:** Computed PHQ-9 score trajectories per participant over time; correlated PHQ-9 quartile with viral suppression status using logistic regression, controlling for age and ART regimen.
5. **Site-level benchmarking:** Funnel plot of site suppression rates vs. expected rates (based on case mix — age, diagnosis year, housing stability); identified outlier sites (>2 SD from expected).
6. **Dashboard:** Looker Studio 4-tab report: (1) Cohort suppression summary; (2) Time-to-suppression Kaplan-Meier; (3) Mental health co-morbidity; (4) Site-level funnel plot.

---

## Tools Used

- **Python (pandas, numpy, scipy, lifelines):** Viral suppression rates, MPR computation, KM survival curves, logistic regression, funnel plot logic
- **lifelines:** Kaplan-Meier estimator, log-rank test for ART regimen comparison
- **Looker Studio:** Interactive cohort dashboard with site drill-down and time-period filter
- **Excel:** Summary tables for OHTN research dissemination reports

---

## Measurable Outcome / Impact

- Viral suppression rate stratification revealed that participants with unstable housing had suppression rates 14 percentage points lower than stably housed participants (68% vs. 82%), quantifying the impact of a key social determinant for program investment
- Kaplan-Meier analysis showed median time-to-suppression was 47 days shorter for INSTI-based regimens vs. NNRTI-based (log-rank p < 0.001), supporting ART prescribing guidelines review
- PHQ-9 quartile logistic regression confirmed high depression burden (PHQ-9 ≥ 15) was independently associated with non-suppression (OR 1.9, 95% CI 1.4–2.6), informing integration of MH supports into HIV care
- Site funnel plot identified 2 sites as statistical outliers (>2 SD below expected suppression), prompting peer review process between high- and low-performing sites
