# R15 audit — H&E and the methylome

Audited 15 August 2026. Run `20260814T211255Z-R15wsimeth-3f8bece`, config
`3d8b9658f679e069bd0cc4fded832dcd636c859cdc5bbf7cc848df5c2c5a5ef3`. Peer-review audit copy. This audit was written as an internal adversarial check before the work
was submitted for review.

## Verdict

Both headline claims stand: morphology carries epigenomic information beyond subtype, and it does
not substitute for the assay. **Two gates were failed and one was re-specified during the work**,
and that history is the most important thing in this file, because a reader of the final numbers
alone would not know the first answer was the opposite one.

## Finding 1 — the sequence of three attempts, recorded in full

This is not tidy and should not be tidied.

| Attempt | Representation | Subtype control | Outcome |
|---|---|---|---|
| 1 | mean-pooled patient embedding | 0.737 | **HALTED** on gate `PC_subtype > 0.90` |
| 2 | attention-MIL, tile-level | 0.799 | **HALTED** on the same gate |
| 3 | attention-MIL, gate re-specified | 0.799 grouped / **0.9703** random | proceeded |

Attempt 1's halt was correct: averaging a whole slide dilutes tumour with stroma, and the
representation was genuinely too weak to interpret anything downstream.

Attempt 2's halt was **the gate's fault, not the pipeline's.** The >0.90 threshold had been taken
from R12's centralized subtype (0.970) and R13's CPTAC transfer (0.979) — both measured with
train and test sharing contributing sites. Demanding that a site-*disjoint* number match a
site-*sharing* benchmark is a category error, and no grouped result was ever going to clear it.

The re-specified gate — reproduce R12's 0.970 under R12's own random-fold protocol — returned
0.9703, within 0.0003. That is a legitimate correction because it changed a *threshold*, not a
criterion, and because the diagnostic that exposed the error is itself reported. Both halted runs
are retained.

**The risk to be honest about:** re-specifying a gate after it fails is exactly the move that
should attract suspicion. What distinguishes this from goalpost-moving is that the replacement
gate is *stricter in kind* (it demands agreement with an external published number rather than
clearing an arbitrary bar), and that the failed threshold's provenance was traced to a specific
identifiable mistake. A reader who disagrees can read the original config in `evidence/`.

## Finding 2 — the first negative was an artefact of our own supervision choice **(PARTLY WITHDRAWN 2026-08-17)**

> **Withdrawn in part by R17.** This finding was titled and written as a general result. It holds
> for the six aggregate targets and **fails genome-wide.** R17 re-ran the identical per-CpG scan on
> the methylation-supervised embedding: 83,369 CpGs above null against 90,137 for the
> subtype-supervised one — not better, marginally worse, median ρ 0.0073 against 0.0057. The
> genome-wide negative is *robust* to supervision.
>
> R17 also supplies the mechanism this finding lacked: a patient's **observed** global mean
> methylation alone puts 74.7% of CpGs above null. Five of the six targets below are means over
> large annotation classes and the sixth is the global mean, so a representation can learn that
> single axis — earning partial ρ 0.135–0.244, which still stands — while carrying nothing at CpG
> resolution.
>
> The text below is left as written, because what it got wrong is the point. See
> `results/R17-percpg-methylation-supervised/`.

The single most consequential correction in this report.

Arm 1 supervised the embedding on subtype, then asked it about methylation, and every readout
agreed it carried nothing:

| Readout | Arm 1 result |
|---|---|
| genome-wide median ρ | 0.006 |
| CpGs above permutation null | 110,212 / 399,579 (27.6%) |
| **subtype label alone** | **220,925 (55.3%)** |
| bulk methylation | ρ = 0.011, p = 0.76 |
| genomic-context enrichment | flat, 0.97–1.01 |

Four independent readings, all consistent, all wrong. The embedding had been optimised to
separate two histological classes; nothing during training asked it to preserve methylation
structure. Supervising directly on methylation gave ρ 0.218–0.317 with partial ρ 0.135–0.244
beyond subtype.

Two things about this deserve recording:

- **A 512-dimensional embedding lost to one binary variable.** That should have been read at the
  time as evidence the representation was wrong for the question, rather than as evidence about
  biology.
- **The flat genomic context was misread** — and then misread a second time. In R10 the KEAP1
  methylation phenotype was strongly structured (island 0.54, open sea 1.44). Arm 1's flatness
  looked like diffuse noise in the *biology*; I reread it as diffuse noise in the *representation*.
  **R17 shows neither reading works:** the methylation-supervised representation is still nearly
  flat (spread 0.212 against 0.047) and tilts the *opposite* way to R10, islands enriched at 1.116
  where R10 had them depleted at 0.54. Flatness of genomic context was not diagnostic of anything
  here, in either direction.

**Generalisable rule:** a negative result from a probe supervised on target A is not evidence
about target B. Before believing any null in this programme, check what the representation was
trained to do.

**Superseding rule (2026-08-17, from R17).** The rule above is sound and I applied it without
testing it. Re-supervising fixed the aggregate targets and changed nothing genome-wide, so
"the supervision was wrong" was a *hypothesis* presented as a conclusion. Two further lessons,
both sharper than the original:

- **Stating a limitation does not bound a claim.** This report's limitations section correctly
  said "per-CpG prediction under methylation supervision was not attempted." That disclosure sat
  three pages behind a summary sentence asserting the answer had been reversed. A reader takes the
  summary. The fix is to narrow the sentence, not to add a caveat.
- **An audit is not automatically the safe document.** This file exists to catch over-reach, and
  it is where the over-reach occurred — the report's own limitations were more careful than its
  audit's headline finding. Being the sceptical document does not confer scepticism.

## Finding 3 — the leakage arm was run because the IP record required it

Not run for the paper; run because a patent should not rest on a number whose naive-analysis
inflation is unknown.

>**Narrowed 2026-08-17 by R18, then re-corrected 2026-08-18 by R21.** This note first said the
> inflation below was "specific to Phikon-v2". That was wrong, and wrong in the same way the finding
> it was correcting had been. R21 measured four histology encoders (Phikon-v2, UNI, Virchow2,
> H-optimus-0): **all four inflate all six targets**, +0.0173 to +0.0849. Only the natural-image
> encoder does not (+0.0004, 3/6). The inflation below is a property of **histology pretraining**,
> not of one vendor, and the number itself is encoder-dependent in magnitude.
> Re-measured on `facebook/dinov2-large` with identical tiles, folds, architecture and
> hyperparameters, the mean becomes **+0.0004** (0.000 on fraction of headroom, against 0.119 here),
> with three of six targets negative. The subtype figure in the same table *does* reproduce and grows
> (+0.2034 vs +0.1710), so this is a genuine dissociation and not a weaker-encoder artefact. Finding
> 3's argument — that the leakage arm was run because a patent should not rest on a number whose
> inflation is unknown — still holds; the inflation is now known to be encoder-dependent, which is a
> sharper version of the same point. See `results/R18-encoder-robustness/`.

| | Mean inflation from site sharing |
|---|---|
| Six methylation targets | **+0.0849** (range +0.0589 to +0.1184) |
| Subtype AUROC | +0.1710 |
| KEAP1 AUROC | +0.0018 |

So a conventional cross-validated analysis of this capability would have reported ρ ≈ 0.34–0.42
where the honest values are 0.22–0.32. Methylation *is* site-confounded — batch effects would
predict that — but less than subtype.

The KEAP1 result being essentially uninflated is a useful internal consistency check: mutation
status has no institutional prevalence structure to exploit, so there is no shortcut available,
and none appears.

## Finding 4 — targets were chosen to be unarguable, and that was the right call

All six are means over **fixed, annotation-defined** probe sets. No component was fitted, no
probe selected on outcome, so there is nothing to accuse of selection bias. The alternative
considered — principal components of the methylation matrix — would have required fitting on the
same patients and invited exactly that objection.

The target means reproduce known biology unprompted (promoters 0.172 and islands 0.226
hypomethylated; open sea 0.668 and gene bodies 0.628 hypermethylated), which is a free check that
the matrix is correctly assembled and joined.

## Finding 5 — tumour purity: I said "no correction available" and was wrong

The first version of this audit named purity as the most likely benign explanation and claimed no
correction was possible, because we lacked per-slide tumour content. That was premature: TCGA
publishes consensus ABSOLUTE purity estimates, and they cover **741 of 760 patients** here.

**The confound is real.** Purity correlates with observed methylation at ρ = −0.306
(p = 1.5e-17) and with the morphology-derived predictions at ρ = −0.218 (p = 2.2e-9). So the
concern was well-founded, not hypothetical.

**Adjusting for it barely moves anything.** Partial ρ controlling for subtype alone versus
subtype + purity:

| Target | subtype only | + purity | change |
|---|---|---|---|
| KEAP1 signature | 0.252 | 0.221 | −0.031 |
| global | 0.218 | 0.210 | −0.008 |
| open sea | 0.228 | 0.218 | −0.010 |
| TSS200 | 0.210 | 0.190 | −0.021 |
| gene body | 0.220 | 0.221 | +0.001 |
| CpG island | 0.134 | 0.136 | +0.002 |

All six remain significant (p ≤ 2.1e-4; four below 1e-8). The largest loss is 0.031.

**Lesson for this series:** "no correction available" was an assertion I had not checked. The
data existed, publicly, in a well-known file. Before writing that phrase again the right move is
a five-minute search — especially in a document intended as an IP record, where an unaddressed
confound is a weakness someone else will find.

**Still unadjusted:** stromal composition, necrosis, proliferation index. Purity was the one with
the strongest prior and the best available measurement, but it is not the whole set.

**Two targets where subtype beats the slide.** Promoter (0.345 vs 0.303) and island (0.288 vs
0.218). The partial correlations remain positive, which supports "different information rather
than more of it" — but it is a weaker claim than the headline numbers suggest in isolation, and
the report states it.

**Mean-pooling across a patient's slides.** No weighting by tumour content or slide quality.

## Limitations to state to the colleague

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
