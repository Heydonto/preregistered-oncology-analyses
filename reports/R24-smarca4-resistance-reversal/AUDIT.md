# R24 audit — the rule returned PARTIAL, and the replication said no

Audited 20 August 2026. Plan `PREREG_PLAN_R24.yaml`, SHA-256 `bffbaad89936c522…`, hashed after
reading only the TPM column headers and before any expression value was compared between groups.
All eight gates passed. Held for IP.

## Verdict

`PARTIAL_one_direction_only`. The knockdown worked (SMARCA4 log₂FC −2.37). Resistance-UP genes fell
below the null as the thesis predicts; resistance-DOWN genes fell too, which is not reversal. In the
second model every sign inverts.

The useful content of this report is that a pre-registered rule returned an inconvenient answer and
a built-in replication then undercut even that.

## Finding 1 — the rule bound, and it bound against the interesting result

Had the decision rule been "resistance-UP genes move below the null", the verdict would be
**THESIS_SUPPORTED** — UP genes sit at −0.0043 against a null of [+0.028, +0.072], comfortably below,
p = 0.0005. That is a clean, publishable-sounding sentence.

The rule declared in advance required **both** directions: UP below the null *and* DOWN above it.
DOWN genes came in at −0.0795, also below. So the mechanical verdict is PARTIAL.

Writing the two-sided rule cost nothing when I wrote it and cost the headline when it fired. That is
the whole argument for pre-registration in one example, and it is worth more to this programme than
a positive result would have been.

**Why the two-sided rule is the right one:** a treatment that lowers *everything* associated with
resistance is not reprogramming the state, it is suppressing a correlated set. Only the up-down
pattern distinguishes reversal from generic suppression, and only the two-sided rule tests for it.

## Finding 2 — the replication arm is what actually settles it

YU005C was in the plan as H3 with no threshold declared, purely as "does this reproduce".

| | resistance UP | resistance DOWN |
|---|---|---|
| PC9-OR | −0.0043 (below null) | −0.0795 (below null) |
| YU005C | **+0.0284 (above null)** | **+0.0322 (above null)** |

Every sign inverts. All four are separable from their nulls at p ≤ 0.004, so this is not noise —
it is two models disagreeing about direction.

**One model showing partial reversal and a second showing the opposite is not a finding about
SMARCA4.** Had the plan included only PC9-OR, this report would have said "partial reversal,
p = 0.0005" and been indefensible. The replication arm cost one extra comparison and changed the
conclusion.

## Finding 3 — the premise is weaker than the test

Cross-line Spearman correlation of the resistance log₂FC vector is **0.17, 0.17 and 0.26** across
PC9/H1975/HCC827. A resistance programme shared across EGFR-mutant lines would sit near 1.

So before asking whether epigenetic modulation reverses *the* resistance programme, one should ask
whether there is one. On this evidence there mostly is not — three lines becoming
osimertinib-resistant change largely different genes.

This was H1 in the plan with no threshold declared, reported descriptively. It should probably have
been a **gate**: if the programme is not shared, the therapeutic thesis has an ill-defined target and
H2 tests something narrower than it appears to. Recorded as a design lesson.

## Finding 4 — the one clean signal, and why it is not the thesis

Binning by magnitude of resistance change gives a monotonic response across five bins:
+0.062, +0.023, −0.083, −0.106, **−0.232**. The more a gene moved on becoming resistant, the more
knockdown lowers it.

Five bins in order, spanning 14,221 down to 171 genes, is not something noise produces easily, and
it is the most interesting number here. But it acts on resistance-associated genes **regardless of
direction**, so it describes suppression of a programme rather than reversal of it — which is
precisely the distinction the two-sided rule was written to catch.

**Recorded as `dose_trend.json`, not only as a figure panel.** A monotonic trend visible in a bar
chart and absent from the evidence files is a claim that cannot be re-checked.

## Finding 5 — effect sizes nobody should get excited about

−0.004 to −0.079 log₂FC. Separable from a 2,000-permutation null and biologically close to nothing.
The report leads with the verdict rather than the p-values for that reason. A significant mean shift
of 0.08 in log₂ space across a gene set is not a therapeutic signal, and describing it as one would
be the same failure as the synthetic figures this work replaced.

## Finding 6 — what this run inherited from R23

R23's plan specified 200 permutations and made its own hypotheses arithmetically untestable under
Holm correction. This plan therefore carries **G5**, a gate asserting that the permutation floor
permits significance at the declared alpha, checked before any test runs. It passed at 0.0005.

A mistake caught in one report becoming a gate in the next is the intended behaviour of this series.

## Limitations to state to the colleague

1. **Two to three replicates per group**, TPM group-mean comparison, no count model.
2. **shRNA off-target effects.** Two hairpins agreeing is supportive, not conclusive.
3. **H2 is not independent** — sets defined in PC9, applied in PC9-OR. YU005C was the independent
   test and it failed.
4. **Transcriptional reversal is not therapy.** No viability or IC50 endpoint here; the plan's
   refusals forbid claiming one.
5. **Cell lines are not patients.**
6. **Hash-before-outcome only** — plan author and executor were the same person.
7. **H1 should have been a gate, not a description.** Recorded so the next plan of this shape makes
   the shared-programme assumption a halt condition.
