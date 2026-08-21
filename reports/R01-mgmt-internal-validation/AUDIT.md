# Audit of R01 — adversarial check, 2026-08-07

Nine checks, run against the shipped evidence bundle only (no re-fitting). The aim was to
break the result, particularly the below-chance point estimate.

## Verdict: R01 stands. Two clarifications added below; no number changes.

| # | Check | Finding | Verdict |
|---|---|---|---|
| A1 | **Label inversion** — a below-chance AUROC can mean the label is flipped | Flipping labels gives elastic net 0.570, RF 0.548, but gradient boosting **0.485**. A genuine inversion would improve *all three*. It does not. | **Not an inversion** |
| A2 | **Degenerate model** — near-constant predictions make AUROC meaningless | Elastic net predictions span 0.264–0.782, sd 0.100, 256 distinct values. RF and GBM similar. | Models are non-degenerate |
| A3 | **Duplicate subjects** | 256 ids, 256 unique, 0 duplicates; all carry the `_11` suffix | Clean |
| A4 | **Label balance vs gate** | 108 methylated / 148 unmethylated — matches gate G1 exactly | Consistent |
| A5 | **Is below-chance significant?** | z = −1.95 vs 0.50, two-sided p = **0.051** | Borderline — see clarification 1 |
| A6 | **Repeat stability** | Elastic net per-repeat interval 0.404–0.511 **spans 0.5**; GBM 0.465–0.545 spans 0.5; RF 0.435–0.492 all below | Not a uniform below-chance effect |
| A7 | **CI method independence** — DeLong is analytic; bootstrap is not | DeLong 0.4299, CI [0.3595, 0.5002]. Bootstrap (4,000 resamples) 0.4290, CI [0.3589, 0.5022] | Methods agree to ~0.001 |
| A8 | **Metric definition** | AUROC of averaged predictions = 0.430; mean of per-repeat AUROCs = 0.459 (Δ +0.029) | See clarification 2 |
| A9 | **Determinism artefact** | Primary vs rerun predictions differ by 0.00e+00 | Bit-identical |

## Clarification 1 — how to describe the below-chance estimate

Taken alone, p = 0.051 is borderline. Taken together with A1 (no consistent inversion), A6
(elastic-net repeats span 0.5), the three-model comparison without multiplicity correction,
and the permutation test (null centred 0.498, observed p = 1.000), the coherent reading is:

> **consistent with chance, with the point estimate happening to fall below it** — not
> evidence of an inverse association.

The report should not, and does not, claim an inverse relationship. Language checked:
R01 says "at or below chance" and "indistinguishable from — indeed slightly worse than —
random", which is accurate. No change required.

## Clarification 2 — the primary metric is the more conservative of two

The pre-specified primary is the DeLong AUROC on out-of-fold predictions **averaged across
repeats** (0.430). The alternative summary, the mean of per-repeat AUROCs, is 0.459. The
averaged-prediction figure is the lower of the two, so a reader could ask whether the worse
number was selected.

It was not selected post hoc: `config.yaml` fixes
`"ci_method": "DeLong on out-of-fold predictions averaged across repeats"` and its SHA-256
(`0ec0a1df…`) is embedded in `results.json`, written before any label was read. Both figures
are reported in the results table.

**The decision is robust either way:** 0.430 and 0.459 are both far below the 0.65
stopping threshold, as are all three models' CI upper bounds (0.500, 0.522, 0.586).

## What the audit could not check

- Whether the CaPTk features were correctly extracted by the data providers. We consumed
  them as published and verified only the archive's byte size and hash.
- Whether the automatic segmentations are anatomically sound. Not inspected; a poor
  segmentation would plausibly contribute to a null, and this is stated as a limitation.
- Whether a different feature family (deep features, perfusion, diffusion) carries signal.
  Out of scope by design.

## Consequence for the programme

R01 is confirmed. The stopping rule stands and M2/M3 remain unjustified **for the MGMT
endpoint**. One methodological gap this audit exposes: a null is only as credible as the
pipeline that produced it, and R01 contains no positive control demonstrating the machinery
can detect an association that is known to exist. **R02 addresses that directly** by testing
age and MGMT against overall survival, both established prognostic factors in
glioblastoma. If those recover and radiomics does not, the R01 null is materially
strengthened.
