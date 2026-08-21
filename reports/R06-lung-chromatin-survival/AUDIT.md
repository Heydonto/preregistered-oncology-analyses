# Audit of R06 — 2026-08-07

## Verdict: the finding stands for ONE cohort. The replication claim is withdrawn.

| # | Attack | Finding | Consequence |
|---|---|---|---|
| F1 | **Event adequacy in the small cohort** | `luad_mskcc_2020`: **45 events / 604 patients (7%)**. Events among mutants: KMT2D 0, SMARCA4 9, KEAP1 5, STK11 8, TP53 18 | Severely underpowered |
| F2 | **Degenerate fit** | KMT2D HR reported as 0.000 [0–inf]: **zero events among 28 mutants** | Not estimable; must not be printed as a number |
| F3 | **Does the small cohort recover its own controls?** | TP53 HR 1.18 (**p=0.58**), STK11 HR 1.04 (**p=0.91**) — both established, both undetected | **Cohort cannot corroborate anything** |
| F4 | **Is SMARCA4's HR 3.97 credible there?** | Built on 9 events among 33 mutants, CI 1.90–8.29 | Small-numbers estimate, not independent confirmation |
| F5 | **Direction consistency** | KEAP1 HR 0.87 (small) vs 2.07 (large) — opposite directions | Signature of underpowering |
| F6 | **Large cohort controls** | STK11 q=0.0004, TP53 q=0.0074, both correct direction and magnitude | Large cohort trustworthy |

## What changed in the report

The run's `_decision` field says "chromatin-regulator mutation(s) replicate as prognostic",
because the pre-specified rule ("significant in both cohorts") was mechanically satisfied.
**The audit overrides this.** A cohort that cannot detect TP53 or STK11 provides no
corroborative evidence, regardless of what its p-value says for SMARCA4. The report states
the result as single-cohort.

This is the third consecutive report where the audit materially altered the headline (R04
reversed, R05 shown seed-dependent, R06 replication withdrawn). The consistent lesson: a
pre-specified rule protects against post-hoc *selection*, but it does not protect against a
cohort being unfit for the question. **Positive controls, not p-values, are what establish
fitness** — and they should gate each cohort individually, not just the analysis as a whole.

## Recommended change to the harness

Future multi-cohort designs should require each cohort to pass its own positive-control gate
*before* its result is allowed to contribute to a replication claim. R06's G3 gate checked
controls only in the larger cohort; had it checked both, the replication claim would never
have been emitted.

## What the audit could not check

- Whether ctDNA detectability confounds the effect sizes (mutation calls require tumour
  shedding, itself adverse). Not separable with these data.
- Treatment effects — no regimen or line-of-therapy adjustment was possible.
- Whether KEAP1 and SMARCA4 effects are independent of each other or of TP53/STK11; only
  univariate models were fitted.
