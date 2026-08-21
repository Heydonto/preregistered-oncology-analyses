# R14 audit — de-novo generated analysis on a sealed cohort

Audited 14 August 2026. Config hash recorded in `evidence/results.json`; 11 gates in
`evidence/gates.json`; the generator's own verbatim self-assessment at
`nsclc-rwpr-study/armA/GENERATOR_SELF_ASSESSMENT.md`.

## Verdict

**The process claim holds and is the report's real contribution.** Hash-before-label ordering,
the author-code firewall, exclusion accounting, null centring and budget were all verified and
all passed. The system recorded a failed positive control rather than hiding it, and returned
POWER-LIMITED-NULL rather than manufacturing a finding — which is the behaviour the study was
built to test.

**One finding needs to reach Paper 2 prominently**, because Paper 2's claim is precisely about
honest pre-registered conduct: a pre-registered hard gate failed, and it was subsequently
reinterpreted. The reinterpretation is defensible but it is a human judgement made after seeing
the result, and it must be presented as such and not folded into the pre-registration.

## Finding 1 — PC2 failed, it was a hard gate, and its failure has no pre-registered remedy

`evidence/gates.json`: 10 of 11 gates pass. **`G3_pc2` fails**, and
`_decision.gates_required_all_pass = False`.

PC2 was: *seg-derived total baseline tumour volume → OS (univariate Cox)*, expecting C > 0.55 and
p < 0.05. Observed: **C-index 0.5035, HR per SD 0.965, Wald p = 0.698** on n=187 with 111 events.
Flat.

The config states PC2's basis explicitly: *"baseline tumor burden is robustly prognostic; proves I
can read the imaging."* So on its own pre-registered terms, **PC2's failure means imaging-reading
competence was not demonstrated.**

Now compare how the two imaging controls were pre-registered:

| Control | Pre-registered failure rule |
|---|---|
| PC3 (PATHF → TPS) | *"FAILURE ⇒ PATHF endpoints reported 'pipeline not competent', not evidence of absence"* |
| G4_pc3_competence | *"PC3 outcome recorded; failure downgrades PATHF claims only"* |
| **PC2 (volume → OS)** | **none — `G3_pc2: "PC2 passes"`, a bare hard gate** |

The generator anticipated PC3 failing and wrote the remedy in advance. It did not do so for PC2.
So when PC2 failed there was no pre-registered instruction for what to do.

**Why this matters concretely:** `M_rad` — a radiomics model built on that same imaging — is one of
only **two** models surviving Holm correction (raw p = 0.0020, Holm p = 0.0100). The report leans
on M_rad while the control that was supposed to establish imaging competence failed.

**The defence, and its status.** The self-assessment offers one: *"PC2's null (target-lesion volume
≠ total burden)"* — i.e. PC2 actually tested whether RECIST *target-lesion* segmentations
approximate *total* tumour burden, a different and genuinely questionable premise, rather than
whether the imaging can be read at all. That is a reasonable reading, and PC1/PC3/PC4 all passed,
so the pipeline is demonstrably competent elsewhere.

But it is a **post-hoc reinterpretation of a pre-registered control**. Paper 2 must say so in
those words. The paper's contribution is a protocol for honest conduct; presenting an after-the-
fact rescue of a failed gate as though it were part of the plan would undercut the very thing
being claimed. Stated plainly it costs nothing and is itself an instructive result: pre-registration
surfaced a control whose basis was mis-specified, which is what pre-registration is for.

Note the verdict was unaffected mechanically: H1 required p < 0.05 *and* G3 passing, and H1's
p was 0.6375, so H1 was never going to be confirmed.

## Finding 2 — PC3 is the novel result, and two different AUROCs are reported for it

Phikon-v2 embeddings of PD-L1 IHC slides predict TPS ≥ 50% with **no IHC training**:

| Statistic | Value |
|---|---|
| AUROC of mean OOF (**headline**) | 0.8697 |
| Mean of per-repeat AUROCs | 0.8564 |
| 95% CI (mean OOF) | 0.8111–0.9221 |
| n / positives | 163 / 67 |

The two definitions differ by 1.3 pp and the report quotes the higher. That is a defensible choice
but it must be stated, and both should appear — the same issue R07 hit, where averaging predictions
across repeats before scoring gives a slightly better number than averaging the scores.

The generator separately disclosed a related mismatch: **the permutation statistic is the seed-0 CV
AUROC while the headline is the 5-repeat mean.** Disclosed, not hidden, but any p-value quoted
beside a 5-repeat headline should note that the null was built on a single split.

This is still the most interesting single observation in the report and the one worth carrying
into Paper 2 as a secondary finding.

## Finding 3 — every modality addition made things worse, consistently

`secondary_paired`, all on matched subsets:

| Comparison | n | Δ AUROC | 95% CI |
|---|---|---|---|
| M_rad vs base | 186 | −0.016 | −0.095 to +0.065 |
| M_pathg vs base | 105 | −0.019 | −0.125 to +0.085 |
| M_pathf vs base | 163 | −0.056 | −0.155 to +0.041 |

**Not one is individually significant** — every CI contains zero. But the sign is negative 3 out
of 3, and H1 (M_full vs M_base, n=73) is also negative at −0.036.

The honest statement is the one the report should make: on cohorts of this size, adding imaging or
pathology to a clinical+PD-L1+TMB baseline does not help, and the consistency of the direction
across four independent comparisons is more informative than any single interval. It is **not**
evidence that the modalities carry no information — it is evidence that at n = 73–186 the added
parameters cost more than they return.

## Finding 4 — H3 excludes zero but must not be read as biology

`H3`: C-index 0.485 (full) vs 0.694 (base), Δ = −0.209, CI −0.353 to −0.064 — genuinely excludes
zero. On **n = 73 with 34 events.**

The generator's own reading is correct and should be adopted verbatim: *"real CI exclusion of 0,
but on n=73 it's a small-sample overfitting statement, not biology."* Paper 2 should quote it that
way.

## Finding 5 — the marginal claims were correctly withheld (no action)

Holm correction across the E1 model family leaves only two survivors:

| Model | Raw p | Holm p | Claimable |
|---|---|---|---|
| E1_M_base | 0.0010 | **0.0060** | yes |
| E1_M_rad | 0.0020 | **0.0100** | yes (see Finding 1) |
| E1_M_pathf | 0.0180 | 0.0719 | no |
| E1_M_clin | 0.0200 | 0.0719 | no |
| E1_M_pathg | 0.0869 | 0.1738 | no |
| E1_M_full | 0.2458 | 0.2458 | no |

M_pathf at Holm 0.072 is not claimable and the report does not claim it. The PATHG author-split
0.765 rests on a single evaluation at n=52. Both were flagged by the generator as marginal and
both were kept out. Correct.

## Finding 6 — process gates verified, and this is the paper's contribution

| Gate | Verified |
|---|---|
| G1_prereg_order | config SHA registered before any sealed-content read |
| G5_no_silent_exclusions | per-endpoint accounting reconciles; each exclusion enumerated with a reason |
| G6_nan_degeneracy | no NaN predictions, no constant OOF vector, 100% OOF coverage |
| G7_null_centred | every archived null within 0.02 of 0.5 (observed 0.494–0.495) |
| G8_pretreatment_only | RAD/TB from baseline scans only; PATHF/PATHG from the diagnostic slide |
| G10_no_leakage | passed |
| G11_firewall | author code sealed for the whole arm |

G7 is the one worth highlighting: nulls centred at 0.494–0.495 across every model family is what
makes the permutation p-values trustworthy rather than decorative.

## Limitations to state to the colleague

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
