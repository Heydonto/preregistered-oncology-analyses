# R22 audit — a gate halted the run, and it was right to

Audited 19 August 2026. Run `20260819T003001Z-R22extfive-def4440`, config
`025add55f381b264bd7d8d3e449e18b47c32ba0de0bb443f9cc311784c7c3bad`. Peer-review audit copy.

## Verdict

`EXTERNAL_SIGNATURE_PRESENT_MAGNITUDE_VARIES`, on R20's pre-declared thresholds. All five gates
passed in the final run. **The first attempt halted**, and that halt is the substance of this report.

## Finding 1 — the reproduction gate caught a statistic, not a bug

The first attempt used a fresh seed. Gate G3 compared its dinov2-large SET_HE value against R20's
0.8768, observed **0.9240**, and stopped the run.

Nothing was broken. The seed sets the `StratifiedKFold` split, and on 76 slides with ~15 per fold a
handful of reassignments move balanced accuracy by several points. Measured over 25 splits:

| Encoder | R20's single split | mean ± sd | range |
|---|---|---|---|
| Phikon-v2 | 1.0000 | 1.0000 ± 0.0000 | 1.0000–1.0000 |
| **dinov2-large** | **0.8768** | **0.9176 ± 0.0247** | **0.8583–0.9528** |
| UNI | 1.0000 | 0.9937 ± 0.0110 | 0.9630–1.0000 |
| H-optimus-0 | 1.0000 | 0.9961 ± 0.0078 | 0.9713–1.0000 |
| Virchow2 | 0.9528 | 0.9601 ± 0.0175 | 0.9157–0.9815 |

R20's 0.8768 is 1.6 standard deviations below its own mean. It was the pre-declared seed, so it was
not chosen to flatter anything — it was simply an unlucky draw reported as a fixed quantity.

**The consequence is a real overstatement.** R20's headline contrast, and the sentence I added to
Paper 1's abstract, said the gap was 1.000 versus 0.877. The honest gap is **1.000 versus 0.918**.
The conclusion — signature present everywhere, total only under some encoders — survives intact.

**Generalisable rule:** a cross-validated metric on fewer than a few hundred samples is a random
variable, and quoting one draw of it is a category error. Report the distribution over splits. This
is distinct from the three earlier failures in this series, which were all "generalised past the arm
I ran"; this one is "treated a noisy estimate as a point."

**What made it visible:** only the gate. Had I not pinned R20's value as a hard expectation, the
fresh-seed 0.9240 would have been written up as the result and the fragility never noticed. The gate
was there to catch a coding difference and caught a statistical one instead.

## Finding 2 — three of five are saturated, which limits what the probe can say

Phikon-v2, UNI and H-optimus-0 all reach 1.0000, and Phikon-v2 does so on all 25 splits with
standard deviation exactly zero. Differences among saturated measurements are not interpretable.

So the external probe is a **presence test**, not a ranking, and the report says so. Any narrative
that orders these five encoders by external separability is reading ceiling noise.

## Finding 3 — the corpus ordering is suggestive and I am not claiming it

Histology means are 0.9601, 0.9937, 0.9961, 1.0000; the natural-image encoder is 0.9176. Same
direction as R21's methylation split.

But R21's split was categorical — 6/6 targets versus 3/6, +0.08 versus +0.0004. Here everything is
above 0.91 and the whole spread is 0.08 of balanced accuracy across ceiling-bound measurements, on
**one** non-pathology point. That is not enough to attribute the gap to corpus, and after being
wrong three times about exactly this kind of inference the report states it as an observation and
stops.

## Finding 4 — the stain shortcut, reproduced deliberately

Adding the archive's 53 IHC slides raises dinov2-large from 0.8768 to 0.9444. A different stain
makes "tell the institutions apart" easier. R12's audit predicted this in writing; R19 partly fell
into it; here it is run on purpose and labelled an upper bound. Running a known trap and showing the
gap is more useful than describing it.

## Limitations to state to the colleague

1. **Two institutions, 76 H&E slides.** Existence, not effect size.
2. **Ceiling effects make three of five incomparable.** Presence test, not a ranking.
3. **One natural-image encoder**, so the corpus ordering is an observation only.
4. **Feature width is uncontrolled** (1024/1536/1280); wider representations give a linear probe
   more room on 76 samples. Width does not order the results, but it is not held constant.
5. **Slide-level mean features**, not the attention-MIL representation.
6. **Stain from filename suffixes**; source mislabelling propagates.
7. **R20's 0.8768 and every downstream quotation of it are corrected** to the 25-split mean.
