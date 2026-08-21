# R08 audit — TCGA lung methylation and the chromatin arm

Audited 7 August 2026. Run dir:
`paper-draft/02-revised-protocol-open-datasets/runs/20260807T130716Z-R08tcgameth-3f8bece`
Config hash (written before any outcome was read):
`07be4547c08128d0b7d02a1e0208a55eec2f90f7972619df8fa26c4de448a240`
Audit script: `evidence/audit_m8.py` → `evidence/audit_addendum.json`

## Verdict

Both headline findings stand: KEAP1 (and SMARCA4) carry a genuine within-histology methylation
phenotype, and methylation adds nothing prognostically. One claim was withdrawn (KMT2D), one
crash was diagnosed and gated, and one confounder was declared unresolvable.

## Gates

| Gate | Expected | Observed | Result |
|---|---|---|---|
| G0_files | 919 | 919 | PASS |
| G1_shape | (486427, 919) | (486427, 919) | PASS |
| G2_replicates | 4 | 4 | PASS |
| G3_tissue_split | 74 normals | 74 | PASS |
| G4_probes_kept | > 300,000 | 392,489 | PASS |
| PC1_tumour_vs_normal | > 2,000 / 20,000 | 16,609 | PASS |
| G5_clinical_join | > 95% | 900/915 | PASS |
| PC2_sex | AUROC > 0.95 | 0.9994 (`cg12026625`) | PASS |
| G6_mutation_api | > 500 records | 1,228 | PASS |
| G7_survival_n | > 700 | 811 | PASS |
| G8_stage_parsed | late fraction 0.20–0.60 | 0.480 | PASS |
| G9_variance_* (×4) | covariate non-constant | 2 levels each | PASS |

PC3 and PC4 were logged rather than gated, by design — they test biology, not data integrity,
so they should not be able to halt the run. Both held anyway (PC3 p=0.0056; PC4 AUROC 0.955).

## Finding 1 — KMT2D withdrawn: the signal was entirely histology

The most consequential finding. Pooled across histologies, KMT2D looked like a strong hit.

| Gene | Pooled | Within LUAD | Within LUSC | Mutation rate LUAD vs LUSC |
|---|---|---|---|---|
| KEAP1 | 24.2% | **29.4%** (n=84) | 4.5% (n=33) | 17.8% vs 8.9%, p=1.9e-4 |
| SMARCA4 | 5.8% | **11.3%** (n=41) | — (n=11) | 8.7% vs 3.0%, p=4.9e-4 |
| KMT2D | 20.7% | **0.0%** (n=33) | 0.03% (n=78) | 7.0% vs 21.1%, p=2.5e-9 |

The mechanism is explicit: PC4 shows LUAD and LUSC are separable from methylation at AUROC
0.955, and KMT2D is three times more common in LUSC. So a pooled mutated-versus-wildtype
comparison is substantially a LUSC-versus-LUAD comparison. Stratifying removes it and nothing
remains.

KEAP1 moves the other way — 24.2% pooled → 29.4% within LUAD — which is the signature of a real
effect that pooling was *diluting* rather than creating.

**Precedent:** this is the R04 failure mode (a molecular panel appeared to beat imaging until
grade mix was controlled) recurring in a different dataset. Both were caught by asking whether
the comparison groups differ in something other than the variable of interest. That check should
be standard in any future report in this series.

## Finding 2 — the first run crashed, and the crash was the good outcome

The stage covariate was parsed with `.replace("Stage ", "")`, but cBioPortal writes
`"STAGE IIIA"` in uppercase. The parse silently failed for every patient, `late_stage` became
constant zero, and the Cox design matrix went singular — `LinAlgError` / `ConvergenceError`.

Two things worth recording:

1. **It failed loudly, but only by luck.** A constant covariate happens to make Cox singular.
   Had the same bug hit a covariate used only in a correlation or a mean comparison, it would
   have produced a plausible wrong number instead of a crash.
2. **Relying on that luck is not a control**, so two gates were added: `G8_stage_parsed`
   (late-stage fraction must fall in 0.20–0.60; observed 0.480) and `G9_variance_*` (every
   model covariate must be non-constant). Either would have caught it directly.

The crashed run is retained at `evidence/gates_FIRST_CRASHED_RUN.json`. A secondary fix in the
same pass: `cv_cox` now drops covariates that are constant *within a fold* and logs the drop,
because SMARCA4 has only 52 carriers and a fold can legitimately contain none.

## Finding 3 — smoking cannot be separated from KEAP1 mutation

Declared as an unresolved limitation rather than adjusted away.

| Gene | Mutated in ever-smokers | in never-smokers | OR | p |
|---|---|---|---|---|
| KEAP1 | 109/726 | 4/81 | 3.40 | 0.011 |
| KMT2D | 104/726 | 4/81 | 3.22 | 0.016 |
| SMARCA4 | 48/726 | 2/81 | 2.80 | 0.22 |

KEAP1 mutation is associated with ever-smoking, and smoking demonstrably alters methylation in
this very cohort (PC3, p=0.0056). Stratifying by smoking is impossible: only **4** of 81
never-smokers carry a KEAP1 mutation, so the never-smoker stratum cannot support a
genome-wide comparison.

Therefore the within-LUAD KEAP1 methylation phenotype is reported as **real but not established
as independent of smoking**. Its magnitude (29.4% of probes, min q=1.7e-20) makes a pure
smoking effect implausible, but this report cannot demonstrate that, and it should not be
written as though it had.

## Finding 4 — Q3 is seed-stable (no correction needed)

R05's conclusion turned on a seed artefact (LR p=0.054 was one draw from 0.0045–0.274), so the
same check was pre-registered here. Across 8 cross-validation seeds the combined model held at
C-index 0.586–0.609 (median 0.595), and the likelihood-ratio p stayed at 0.854. The Q3 null is
not a split artefact.

## Finding 5 — Q2's null is informative, unlike R07's

Worth stating explicitly because the series now contains both kinds of null.

| | R07 | R08 Q2 |
|---|---|---|
| Cohort | 69 patients, 29 events-equivalent | 796 patients, 313 events |
| Observed | AUROC 0.555 [0.42, 0.69] | C-index 0.486 [0.450, 0.523] |
| Detectable at 80% power | ≥ 0.685 | ≥ 0.552 |
| Reading | inconclusive | **genuine null** |

R08 excludes any methylation prognostic signature above C-index 0.552. That is a substantive
negative result, and it is the basis for the recommendation not to fund further methylation
acquisition for NSCLC prognosis.

## Reproducibility checks

- Config written and hashed before any outcome was read
- Acquisition: all 919 files MD5-verified against GDC's published values; a single probe vector
  SHA-256-verified identical across all 919; fewer than 919 was a fatal condition, not a
  smaller matrix
- Platform restriction (450k only) applied *before* download, not after seeing results
- Technical replicates resolved by a deterministic rule (lowest `file_id`) fixed in the config
- Probe NA filter computed on tumours only, so normals cannot influence probe inclusion
- Feature selection and hyperparameters fitted strictly inside training folds
- Both runs retained, including the crashed one
- All 25 numbers in the report verified programmatically against `evidence/*.json`

## Limitations to state to the colleague

1. **KEAP1's methylation phenotype is confounded with smoking** and cannot be separated in this
   cohort (Finding 3). This is the main caveat on the positive result.
2. **KEAP1's *prognostic* effect does not appear here** (HR 1.085, p=0.61) against R06's HR 2.07
   in the MSK ctDNA cohort. R08 is not the right test — its model is stage- and
   histology-adjusted, and TCGA is largely resected early-stage disease against MSK's advanced
   disease. R09 tests it on R06's own terms.
3. **Q1 uses 20,000 sampled probes, not all 392,489**, for tractability. The sample is drawn
   once under a fixed seed and used identically for every gene and every stratum, so
   comparisons between them are internally consistent; absolute percentages carry sampling
   error of roughly ±0.6pp.
4. **SMARCA4 within LUSC could not be tested** (11 carriers). Its within-histology result rests
   on LUAD alone.
5. **Array generations were not combined.** 311 HM27 and 53 EPIC v2 files exist for these
   projects and were excluded. Including them would add samples at the cost of a platform
   confound; that trade was declined in advance.
