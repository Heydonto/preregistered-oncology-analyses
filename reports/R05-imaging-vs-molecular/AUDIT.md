# Audit of R05 — 2026-08-07

## Verdict: the reported conclusion is correct, but for a firmer reason than the run gave.

The run produced LR p = 0.054 and, by the pre-specified rule, concluded imaging does not add
beyond the molecular markers. The audit shows that conclusion is right but that **a single
p-value near the threshold should never have been the basis for it.**

| # | Attack | Finding | Consequence |
|---|---|---|---|
| E1 | **Is the "molecular panel" really two markers?** | IDH mutant n = **5 of 247**. Adding IDH to MGMT: LR χ²=1.45, **p=0.228** | The panel is **MGMT alone**. Report must say so |
| E2 | **Is p=0.054 stable?** Borderline results are often seed artefacts | Across 8 CV seeds LR p spans **0.0045–0.274**, median **0.089**; clears 0.05 in only **2 of 8**. Imaging C-index 0.512–0.549 | **p=0.054 was coincidental.** The honest summary is the seed distribution, not one draw |
| E3 | **Method consistency with R02** | R02 used 10× *repeated* 5-fold and averaged predictions; R05 used a single 5-fold | R05's protocol was **less robust than R02's** — my inconsistency, and the direct cause of the seed sensitivity in E2 |

## What the report says as a result

Not "p=0.054, therefore no". Instead: **across eight cross-validation seeds the incremental
value of imaging beyond MGMT is not established (median p=0.089, significant in 2 of 8)**,
and imaging alone reaches only C-index 0.512–0.549 in this subset. The multi-seed
distribution *is* the result.

## Reconciling with R02

R02 reported imaging C-index 0.602 on 574 patients; here it is ~0.535 on 247. Two
differences: the subset requires a definite MGMT call, and it is 43% of the size. R02's
figure came from 10 averaged repeats and remains the better estimate **for its cohort**;
R05's is the better estimate for the head-to-head subset. Both are reported, neither
supersedes the other, and the discrepancy is stated rather than smoothed.

## What the audit could not check

- Whether the MGMT-defined subset is systematically different from the other 327 patients in
  ways that weaken imaging — plausible and untested.
- Whether a repeated-CV rerun of R05 would shift the median p. E2 approximates this with 8
  seeds; a full 10×5 repeated design would be the clean fix and is recommended if this
  analysis is ever put in a manuscript.
- IDH cannot be evaluated at all here (5 mutants).
