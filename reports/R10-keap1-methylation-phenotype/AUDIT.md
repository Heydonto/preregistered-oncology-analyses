# R10 audit — the KEAP1 methylation phenotype

Audited 7 August 2026. Run dir:
`paper-draft/02-revised-protocol-open-datasets/runs/20260807T135447Z-R10keap1pheno-3f8bece`
Config hash (written before any outcome was read):
`229f8cb6cc4312e865d78e3560c143f26a3a5693af86b3fde4477c1960db0f80`
Audit script: `evidence/audit_m11.py` → `evidence/audit_addendum.json`

## Verdict

All four claims survive adversarial checking, and no correction was needed. This is the first
report in the series where that is true, so it is worth saying why the checks were not toothless:
each was aimed at a specific mechanism by which the claim could have been spurious, and two of
them *strengthened* the result rather than merely failing to break it.

## Gates

| Gate | Expected | Observed | Result |
|---|---|---|---|
| G0_probes | 392,489 (identical to R08's filter) | 392,489 | PASS |
| G1_luad | 471 LUAD tumours | 471 | PASS |
| G2_keap1_luad | 84 mutants | 84 | PASS |
| G3_signature_reproduces_R08 | within 5pp of 0.294 | 0.2944 | PASS |
| G4_positive_control_sex | AUROC > 0.95 | 1.000 | PASS |
| G5_annotation_coverage | > 90% annotated | 390,543/392,489 (99.5%) | PASS |

G3 matters more than it looks: R08 estimated the signature from a 20,000-probe sample and got
29.35%; the full 392,489 probes give 29.44%. R08's sampling error was ~0.1pp, so its ±0.6pp
stated uncertainty was conservative.

## Check A — is the NRF2 enrichment a gene-size artefact?

The run's null drew random gene sets matched on the **number of genes**, not on probes-per-gene.
If NRF2 genes carry more probes than average, and probe count correlates with the differential
fraction, the enrichment would be a size effect.

| | Mean fraction differential | p |
|---|---|---|
| NRF2 pathway (33 genes) | 0.357 | — |
| Unmatched null (run) | 0.243 | 0.006 |
| **Probe-count-matched null (audit)** | 0.239 | **0.003** |

The enrichment **strengthened** under matching. The premise is also weak on its own terms:
probes-per-gene correlates with the differential fraction at Spearman ρ = 0.089 — significant
at this n (p=1.3e-32) but explaining under 1% of variance.

Not corrected. The report quotes the matched p-value (0.003) as primary, which is the more
conservative construction.

## Check B — Q4's overlap needed a chance baseline

"12.4% shared" is meaningless without knowing what two independent signatures of these sizes
would share anyway. The run reported the fraction and an odds ratio but not the expected count.

| | Probes |
|---|---|
| KEAP1 signature | 115,557 |
| Smoking signature (KEAP1-wildtype only) | 37,008 |
| Observed overlap | 14,306 |
| **Expected if independent** | **10,896** |
| Maximum possible (= smaller signature) | 37,008 |

Observed/expected = **1.31×**; observed/maximum = 0.39. So the two signatures are only slightly
more similar than chance. The report now states the expected count alongside the observed one,
which is the number that makes 12.4% interpretable.

## Check C — is the classifier reading smoking rather than KEAP1?

The decisive test, and it does not depend on the smoking signature being complete: delete every
one of the 37,008 smoking-associated probes and retrain from scratch.

| | AUROC |
|---|---|
| All 392,489 probes | 0.910 |
| **355,481 probes, all smoking-shared deleted** | **0.904** |

A 0.006 loss. The KEAP1 signal does not live in the probes it shares with smoking.

## Check D — is 0.910 a fold-split artefact?

R05's entire conclusion turned on a seed artefact (LR p=0.054 was one draw from a 0.0045–0.274
range), so this is now a standing check.

AUROC across 5 seeds: 0.910, 0.868, 0.910, 0.892, 0.898 — median **0.898**, range 0.868–0.910.

Stable. The report quotes 0.910 as the primary (seed 20260807, pre-registered) and shows the
full range in the figure, so the reader is not misled by the best draw.

## Why Q4 worked when R08 said it could not

Worth recording as a method note. R08 concluded the smoking confound was unresolvable because
only 4 of 81 never-smokers carried a KEAP1 mutation — true, and fatal for any design that holds
smoking constant while varying KEAP1.

Inverting the design removes the constraint: measure the smoking signature in KEAP1-wildtype
tumours (where never-smokers are plentiful: 65 of 373), then ask whether the KEAP1 signature is
made of those probes. This needs no never-smoker who carries the mutation.

R08's limitation was a limitation of the design it had reached for, not of the data. That is
worth remembering for the other confounds this series has declared unresolvable.

## Limitations to state to the colleague

1. **The smoking signature rests on 65 never-smokers**, so it is measured less precisely than
   its 37,008 probes suggest. If systematically under-detected, the true overlap exceeds 12.4%.
   Two bounds: the limiting group sizes are comparable (65 never-smokers vs 84 KEAP1 mutants), so
   the signatures are not grossly unequally powered; and Check C does not require the smoking
   signature to be complete.
2. **Association, not mechanism.** Nothing here shows KEAP1 mutation *causes* the methylation
   pattern. Tumour purity, copy-number and proliferation differences between mutant and wildtype
   tumours are not adjusted for, and any could contribute.
3. **The top-gene list is descriptive only.** The ranking favours genes with dense probe
   coverage, and no gene-level claim is made. It is included because readers will ask.
4. **NRF2 enrichment is a set-level result.** It does not establish that any individual NRF2
   target is differentially methylated in a functionally meaningful way, and the analysis says
   nothing about direction (hyper- vs hypo-methylation) — that was not tested.
5. **No prognostic implication.** R08 showed methylation adds nothing beyond mutation, stage and
   age (LR p=0.85). A strong phenotype and no prognostic value coexist without contradiction:
   KEAP1 status is already knowable from sequencing, so predicting it from methylation has no
   standalone clinical use.
6. **Single cohort, single platform.** TCGA-LUAD 450k only. Nothing here is replicated, and the
   series has learned (R09) not to assume a replication cohort exists.
