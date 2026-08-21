# Audit of R02 — 2026-08-07

Five checks run after the result, against the shipped evidence bundle.

## Verdict: R02 stands. The positive finding survives every attack tried.

| # | Attack | Finding | Verdict |
|---|---|---|---|
| B1 | **The signal is just age.** Tumour appearance correlates with age; age is prognostic. | Spearman(prediction, age) = **−0.125**. Age alone C=0.622, radiomics alone C=0.602 — similar magnitude but weakly correlated predictions. | Not explained by age |
| B2 | **No incremental value.** Even if uncorrelated, radiomics may add nothing over age. | Cox: age only C=0.622 (HR 1.35, p=5.8e-11); radiomics only C=0.602 (HR 0.75, p=1.1e-10); **age + radiomics C=0.647** with both terms significant. | Adds value |
| B3 | **Formal test of increment.** | Likelihood-ratio χ²(1) = **36.35, p = 1.6e-9**. | Radiomics adds beyond age |
| B4 | **Convergence.** ElasticNet emitted ConvergenceWarnings at max_iter=5,000. | C-index **0.5978 identically** at 5,000 / 50,000 / 200,000 iterations. The warnings are cosmetic; the solution is stable. | Immaterial |
| B5 | **Leakage.** A positive result demands a null check. | 200 label permutations: null mean **0.4998**, 95% interval [0.459, 0.542], correctly centred on chance. Observed 0.602, **p = 0.005**. | No leakage |

## What the audit changed

Nothing in the numbers. Two things were added to the report because of it: the
age-adjustment analysis (B1–B3) is now reported in the results table rather than left
implicit, and the convergence check is disclosed rather than silently ignored.

## What the audit could not check

- **External validity.** No cohort in our holdings pairs imaging with survival, so the
  result cannot be validated outside UPENN. Stated as the first limitation in the report.
- **Segmentation quality.** Automatic segmentations were used and not inspected.
- **Selection from excluding survivors.** The 17 patients alive have no recorded survival
  time and so are absent. This removes censoring but introduces mild selection that cannot
  be corrected with the data available.
- Whether a different model class would do better. Only the pre-registered elastic net was
  run for the primary; no comparator sweep was performed for survival.

## Relationship to R01

R02 supplies the positive control R01's audit identified as missing. The pipeline recovers
two known associations (age p=2.2e-18, MGMT p=5.8e-8) and finds survival signal at p=0.005,
using the same features and machinery that found nothing for MGMT classification. The R01
null is therefore attributable to the absence of MGMT information in these features, not to
a defective pipeline.
