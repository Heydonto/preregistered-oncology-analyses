# Audit of R04 — 2026-08-07

## Verdict: R04 stands, but only after the audit reversed its headline.

This is the audit that most changed a report. As first computed, R04 said the molecular
panel (C-index 0.797) far outperforms imaging (0.602 from R02). **That comparison was
invalid and would have misdirected the glioma arm.**

| # | Attack | Finding | Consequence |
|---|---|---|---|
| D1 | **Are the cohorts comparable?** | TCGA pan-glioma is 589 glioblastoma, 171 oligodendroglioma, 167 astrocytoma, 113 oligoastrocytoma. **UPENN-GBM (R02) is glioblastoma only.** | Comparison invalid as framed |
| D2 | **Restrict TCGA to glioblastoma, matching UPENN** | IDH only 0.719 → **0.541**; IDH+MGMT 0.747 → **0.561**; IDH+MGMT+TERT not estimable (n=26) | **Headline reversed.** Within GBM the panel (0.541–0.561) is *below* imaging (0.602) |
| D3 | **Is the gap just overfitting?** | 5-fold CV: 0.719 → 0.707; 0.797 → 0.759 | No — optimism is modest; cohort composition is the driver |
| D4 | **Positive control** | IDH: median 89.7 vs 14.0 months, p=5.7e-71, correct direction | Pipeline sound |
| D5 | **Biological sanity of co-occurrence** | IDH×TERT mutually exclusive (OR 0.186); IDH×MGMT co-occurring (OR 17.1) — both match known biology (G-CIMP) | Joins and data sound |
| D6 | **Independent replication** | MSK n=923: TERT HR 1.69 (p=2.8e-6), CDKN2A HR 2.37 (p=4.3e-14) | Markers prognostic in a second cohort |

## Why the original framing was wrong

IDH status is close to a proxy for whether a diffuse glioma is low-grade or glioblastoma,
and that distinction dominates survival. A C-index computed across grades therefore measures
mostly "can we tell a low-grade glioma from a GBM", which is not the clinical question the
imaging work addresses. Restricting to one disease removes the tautology and the apparent
advantage disappears.

**The corrected finding is more interesting than the wrong one:** within glioblastoma the
standard molecular panel is weakly prognostic, which is precisely why imaging is worth
studying there.

## Gate failure, retained

First attempt halted at G2: expected 1,046, observed 1,040. Six patients have
`OS_MONTHS = 0`, unusable in a Cox model. The expectation was too tight, not the data wrong.
Corrected to 1,040 with the reason in source; failed run kept as
`gates_FAILED_first_attempt.json`.

## What the audit could not check

- **No within-cohort head-to-head.** Imaging and molecular results come from different
  populations. UPENN does carry IDH1 and MGMT, so a direct head-to-head on those 574
  patients is possible and is the obvious next analysis — it would settle the comparison
  properly rather than by matching on disease alone.
- TERT within glioblastoma cannot be evaluated at all (n=26).
- TCGA treatment is heterogeneous and pre-dates current standards; survival differences
  partly reflect era.
- Panel C-indices come from Cox on binary markers, not a tuned learner, so they are a floor
  rather than a ceiling for what the markers could achieve.
