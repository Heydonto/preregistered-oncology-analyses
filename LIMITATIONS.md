# Limitations, consolidated

Pulled verbatim from each report's audit. Nothing here is a summary or a paraphrase; if a
limitation reads awkwardly it reads that way in the source.


## R07-methylation-tki-response

1. **Underpowered by design, and quantified.** 69 patients against ~465k probes. Detectable
   only at AUROC ≥ 0.685. To reach AUROC 0.60 at 80% power needs 263 patients.
2. **SD patients excluded** (10 of 79), pre-declared. Response is dichotomised PR vs PD, so
   the analysis speaks to the extremes of the response distribution, not to a gradient.
3. **Single cohort, single platform.** No replication cohort exists in our holdings; no batch
   or array-position covariates are available in the series metadata, so technical confounding
   cannot be adjusted for — though it would tend to *create* apparent signal, not remove it,
   which makes the null more credible rather than less.
4. **Smoking control unresolved.** Left as-is rather than substituted for a better one, since
   selecting a control after seeing it fail would defeat its purpose.

## R08-tcga-methylation-chromatin

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

## R09-replication-attempt

1. **No institutional independence was achievable.** Both eligible cohorts were MSK, as R06's
   was. Patients are disjoint and gated as such, but that is a weaker claim and the report says
   so explicitly rather than implying more.
2. **The two eligible cohorts are heavily selected** — brain metastases and PDX-derived tumours.
   Matching R06's median OS does not make them equivalent populations, and their control failure
   may reflect selection rather than setting.
3. **This is not evidence against R06.** No cohort demonstrated the competence to test it. The
   report states the bounded claim.
4. **The screen covers cBioPortal only.** Cohorts in dbGaP, EGA or institutional holdings were
   not screened, because they are not openly accessible. A qualifying cohort may well exist
   behind controlled access — §5 of the report specifies exactly what to look for.
5. **`lung_msk_pdx` overlap could not be verified by a shared identifier scheme** beyond patient
   ID matching, which returned 0. If that study re-identifies patients differently, undetected
   overlap is possible; it would only bias toward *agreement* with R06, and none was observed.

## R10-keap1-methylation-phenotype

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

## R12-federated-simulation

1. **The federation is simulated.** Weights never crossed a real institutional boundary; no IT,
   governance or scanner heterogeneity was involved. The 301 MiB-per-run communication figure is
   a real measurement of the protocol, not of a deployment.
2. **Only subtype supports the gap-recovery claim.** EGFR and OS are weak signals for every
   method tested, and the report is right not to claim recovery there.
3. **EGFR is power-limited** at 74 mutants (Finding 2).
4. **Five silos, one collection** (Finding 5).
5. **Attention-MIL over Phikon-v2 features throughout** — conclusions may not transfer to other
   encoders, and the site signature in particular could be encoder-specific.

## R13-site-adaptation

1. **Subtype k=0 is a strong, clean result** and can be stated without hedging.
2. **EGFR: report the paired delta (+0.042, 7/8 seeds), not the means**, and always with the
   36–43 patient evaluation size.
3. **EGFR is power-limited** across the programme (R12: 74 mutants, MDE 0.584).
4. **One external cohort.** CPTAC only; no claim about site-transfer in general.
5. **The methylation track is cross-platform** (450k→EPIC, 93.53% of loci transferable) and its
   own audit verified R10 reproduction bit-identically — but it is a different data modality from
   the WSI track and the two should not be pooled into a single "adaptation" claim.

## R14-autonomous-generation

1. **A pre-registered hard gate failed and was reinterpreted afterwards** (Finding 1). Must be
   stated as a post-hoc reinterpretation, in those words.
2. **One sealed cohort, 247 patients.** Three sealed holdouts remain unconsumed (ACRIN PET, LUSC
   imaging mass spec, Anti-PD-1-Lung imaging); the generalisation of the process claim rests on
   running them.
3. **H1 was structurally underpowered from the start.** The generator's own sharpest self-criticism:
   modality availability is confounded with cohort membership, so the n=73 all-modality
   intersection was never going to have power. A stronger design would have pre-registered the
   pairwise intersections as primary.
4. **Two AUROC definitions and a permutation-statistic mismatch** on PC3 (Finding 2).
5. **The firewall is self-imposed.** For R14 it was one agent restraining itself. Paper 2's claim is
   materially stronger if the generator is a structurally separate agent that cannot see labels —
   which is how the three remaining holdouts will be run.

## R15-wsi-methylation

1. **Association, not measurement.** Partial ρ 0.135–0.244 explains a few percent of variance.
   Survives adjustment for tumour purity (Finding 5).
2. **A slide cannot replace the assay** — KEAP1 0.664 from morphology against 0.910 from the
   array.
3. **One cohort, one encoder, one platform.** TCGA lung / Phikon-v2 / Illumina 450k, nothing
   replicated.
4. **The Arm-1 negative is bounded** to subtype-supervised embeddings and must not be quoted as a
   general result about H&E.
5. **Naive validation would overstate this by ~0.085 in ρ**, which is the number that matters if
   anyone else claims a similar capability from single-collection cross-validation.

## R16-acrin-arm-a

1. **Single trial, tabular only.** Core-lab PET assessments, not the 135.5 GiB of PET/CT imaging.
   The generalisation claim is about an unseen cohort and data shape, not an unseen imaging pipeline.
2. **The exposure is 59% floored.**
3. **Two of four controls did not fully corroborate.**
4. **G7 passed at exactly its threshold.**
5. **Not an independent replication.** The published ACRIN 6668 analysis used the same data; the
   claim is about the *process* recovering a known result blind, not about establishing the biology.
6. **The central negative result is about our own gates.** Twelve gates now pass; eight passed while
   the reported primary was a control's hazard ratio. Anyone reusing this harness should assume its
   gate set is still incomplete in ways not yet discovered.

## R17-percpg-methylation-supervised

1. **This does not overturn R15's headline.** Aggregate partial ρ 0.135–0.244 beyond subtype,
   surviving purity adjustment, stands. What is bounded is the genome-wide per-CpG claim.
2. **One architecture, one alpha.** A different readout or per-CpG tuning is not excluded, but it
   would need pre-declaration to count.
3. **Site-grouped folds only.** No random-fold variant; the leakage question was settled in R15.
4. **The negative is about resolution, not existence.** Morphology tracks a global methylation
   axis weakly. The claim it does not support is CpG-level prediction.
5. **Two sentences in R15's report and one finding in R15's audit are now wrong** and have been
   corrected in place, with the original text quoted so the change is visible.

## R18-encoder-robustness

1. **Two encoders is not a survey.** This separates "property of Phikon-v2" from "property of the
   archive". It says nothing about the pathology foundation-model field. UNI, H-optimus-0 and
   Virchow2 are behind licence acceptance that was not ours to give.
2. **Patch size is not matched** — 14 against 16, so 256 tokens per tile against 196. Class,
   width, depth, objective, readout and normalisation are matched; patch size cannot be without
   retraining.
3. **Absolute performance is not a comparison of encoders.** A natural-image model doing worse at
   lung subtype is expected and uninteresting. The claim is about the gap *within* each encoder.
4. **One cohort** — TCGA lung, the same 760 patients, 67 sites, same slides, same tiles.
5. **The external site probe was not re-run under dinov2.** Features are staged
   (`/vol/dinov2-external`, 129 npz) but R19's EAGLE-vs-HTMCP probe remains Phikon-v2 only. That is
   the obvious next arm and it is not done.
6. **Paper 1's methylation-leakage claim is now encoder-specific** and has been corrected in place.

## R19-external-site-probe

1. **Stain came from filename suffixes** (`HE` / `CHR` / `P16`), the only stain field in the staged
   features. Source mislabelling would propagate.
2. **75 slides, two institutions.** Bounds existence, not effect size.
3. **Slide-level mean features**, not the attention-MIL representation.
4. **Row 3's 0.014 discrepancy is permanently unattributable.**
5. **Phikon-v2 only** — answered by R20: the signature holds under a second encoder (0.8768,
   p = 0.0050) but not perfectly, so "total separation" is encoder-specific.
6. **Two Paper 1 rows were unevidenced until today**, and the "23/23 verified" statement covering
   them was wrong. Corrected in Paper 1 in place.

## R20-external-probe-two-encoders

1. **Two institutions, 76 H&E slides.** Bounds existence, not effect size.
2. **Slide-level mean features**, not the attention-MIL representation.
3. **Stain from filename suffixes** (`HE` / `CHR` / `P16`) — source mislabelling propagates.
4. **1.000 and 0.877 are not strictly commensurable.** dinov2-large is weaker at every histology
   task here, so a lower separation is expected in direction. What the arms establish is that both
   clear a permutation null, not that the gap between them is a calibrated quantity.
5. **Two encoders is not a survey** — same gating limitation as R18.
6. **Paper 1's "perfectly" and R19's slide-set description are corrected in place.**

## R21-five-encoder-survey

1. **Corpus is confounded with vintage and capacity.** The three pathology models are newer, larger
   and better-trained than Phikon-v2; dinov2-large is the only non-pathology model. Corpus is the
   most parsimonious reading of the split, not the only one.
2. **Normalisation is not constant** — H-optimus-0 prescribes its own. Correct practice, looser
   design than R18's.
3. **CLS-only readout** may understate Virchow2, whose authors recommend CLS + mean-patch.
4. **One archive, one cohort, one set of tiles.** 760 patients, 67 sites.
5. **The fivefold spread among histology encoders** (+0.0173 to +0.0849) is unexplained.
6. **Two of five are CC-BY-NC-ND.** Lead with the Apache-2.0 arm in anything commercial.
7. **R18's secondary verdict and the narrowing note added to R15 are both corrected in place**,
   with the original wording quoted.

## R22-external-probe-five-encoders

1. **Two institutions, 76 H&E slides.** Existence, not effect size.
2. **Ceiling effects make three of five incomparable.** Presence test, not a ranking.
3. **One natural-image encoder**, so the corpus ordering is an observation only.
4. **Feature width is uncontrolled** (1024/1536/1280); wider representations give a linear probe
   more room on 76 samples. Width does not order the results, but it is not held constant.
5. **Slide-level mean features**, not the attention-MIL representation.
6. **Stain from filename suffixes**; source mislabelling propagates.
7. **R20's 0.8768 and every downstream quotation of it are corrected** to the 25-split mean.

## R23-methylation-drug-resistance

1. **Cell lines are not tumours.** Methylation drifts in culture; ln IC50 is a lab phenotype.
2. **The primary question is unresolved.** A corrected pre-registration should build both indices
   on general-sensitivity residuals and declare that in advance.
3. **Per-drug associations are not drug-specific** and are not claimed to be.
4. **Mixed histology** — 58 small cell and 21 mesothelioma among 187; no stratified analysis was
   pre-declared.
5. **Sex-chromosome probes were not excluded**; the 450k manifest is not staged locally, so
   exclusion E5 went unapplied.
6. **Hash-before-label only.** Plan author and executor were the same person — the weaker form of
   the protocol, not R16's structural isolation.

## R24-smarca4-resistance-reversal

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
