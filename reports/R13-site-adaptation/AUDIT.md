# R13 audit — site adaptation on a new hospital's data

Audited 14 August 2026. Two tracks with separate evidence bundles: `evidence/armB-wsi/` (CPTAC,
432 patients) and `evidence/armB-meth/` (anti-PD-1 EPIC methylation).

## Verdict

**The subtype claim is fully supported and unusually clean.** The EGFR claim is directionally
supported but its point estimates are fragile, and the report must carry the evaluation-set size
wherever it quotes an EGFR number. The methylation track was already audited on 2026-08-10 and
that audit is stronger than most in this series.

## Finding 1 — the methylation track's existing audit is exemplary (no action)

`evidence/armB-meth/audit_armb_meth.json` verifies R10 reproduction **bit-identically**:

| | Value |
|---|---|
| AUROC recomputed from OOF | 0.9098375784422296 |
| R10 published | 0.9098375784422296 |
| Absolute difference | 0.0 |
| OOF predictions bit-identical to R10 | **True** (sha256 `976b2527…cf5b60` both sides) |
| Config hash matches file bytes | True |

It also independently recounts the array overlap: 485,577 450k loci, 867,926 EPIC, 454,181
shared — 93.53% of 450k transferable. That is the right way to verify a cross-platform claim, and
it was done before I looked. Nothing to add.

## Finding 2 — the subtype "k = 0" claim is solid

Eight seeds, evaluation set fixed at 173 patients throughout:

| Method | k | Mean AUROC | SD | Range |
|---|---|---|---|---|
| zero-shot | 0 | **0.981** | 0.008 | 0.970–0.991 |
| fine-tune | 10 | 0.982 | 0.008 | 0.972–0.994 |
| fine-tune | 100 | 0.985 | 0.008 | 0.975–0.996 |
| linear probe | 100 | 0.985 | 0.008 | 0.973–0.994 |

The seed SD is 0.008 and the entire k=0→100 gain is +0.004 — within noise. So "a new hospital
needs zero labelled cases for subtype" is exactly what the data say. The fixed 173-patient
evaluation set across all seeds is what makes this comparison clean.

## Finding 3 — the EGFR curve: direction real, point estimates fragile

This is the finding that needs to reach the paper.

| Method | k | Mean | SD | Min–Max | Range | n_eval |
|---|---|---|---|---|---|---|
| zero-shot | 0 | 0.712 | 0.088 | 0.562–0.872 | 0.309 | 36–43 |
| fine-tune | 10 | 0.736 | 0.092 | 0.586–0.885 | 0.300 | 36–43 |
| fine-tune | 25 | 0.764 | 0.090 | 0.613–0.913 | 0.300 | 36–43 |
| fine-tune | 100 | 0.778 | 0.085 | 0.633–0.913 | 0.280 | 36–43 |

**The unpaired seed spread (0.30) is seven times the k-effect (0.042).** Taken at face value that
would make the learning curve uninterpretable — the R05 failure mode.

It is rescued by pairing, and this is the correct analysis: the seed variance is *shared* across k
because each seed fixes the evaluation set, so a within-seed comparison removes it.

- paired deltas k=10→k=100: +0.026, +0.097, +0.028, −0.010, +0.047, +0.041, +0.039, +0.067
- **7 of 8 seeds improved**, mean paired delta **+0.042**

So more labels do help EGFR, consistently. What is *not* supported is any particular absolute
value: with 36–43 evaluation patients and an SD of 0.088, "0.764" is a point on a wide
distribution.

**Two requirements for the paper.** Quote the paired delta and the seed count, not the bare means.
And state the evaluation-set size (36–43 patients) every time an EGFR AUROC appears — a reader
seeing "0.778" beside subtype's "0.985" will otherwise assume comparable precision, when one rests
on 173 patients and the other on ~40.

## Finding 4 — the EGFR evaluation set shrinks with k, and varies by seed

`n_eval_patients` is fixed at 173 for subtype but moves 36–43 for EGFR. The k labelled patients
are drawn from the evaluation pool, so each seed evaluates a different, smaller set. Pairing
handles the variance (Finding 3), but it means the EGFR arm evaluates on roughly a tenth of the
cohort that the subtype arm uses.

Combined with R12's independent EGFR power finding — 74 mutants, minimum detectable AUROC 0.584 —
the consistent picture across both reports is that **EGFR is the weak endpoint in this programme
and every EGFR statement should be power-qualified.**

## Finding 5 — what "a new hospital" means here

CPTAC is a genuinely independent institution set, which is the strength of this report and the
thing R12 could not offer. But it is one external cohort, so "what does the model do on day one
at a new hospital" is answered for *a* new hospital, not for hospitals in general. The paper
should say CPTAC, not "a new site", wherever the claim is made.

## Limitations to state to the colleague

1. **Subtype k=0 is a strong, clean result** and can be stated without hedging.
2. **EGFR: report the paired delta (+0.042, 7/8 seeds), not the means**, and always with the
   36–43 patient evaluation size.
3. **EGFR is power-limited** across the programme (R12: 74 mutants, MDE 0.584).
4. **One external cohort.** CPTAC only; no claim about site-transfer in general.
5. **The methylation track is cross-platform** (450k→EPIC, 93.53% of loci transferable) and its
   own audit verified R10 reproduction bit-identically — but it is a different data modality from
   the WSI track and the two should not be pooled into a single "adaptation" claim.
