# R21 audit — the third time I generalised from too few arms

Audited 18 August 2026. Run `20260818T124806Z-R21fiveenc-03974b7`, config
`b23ba9307a6081c4f2d18460252862254084d899d9c344163dd4f0fc13ebcc54`. Held for IP.

## Verdict

All gates passed. `SUBTYPE_LEAKAGE_UNIVERSAL` and
`METHYLATION_LEAKAGE_IS_A_PROPERTY_OF_HISTOLOGY_PRETRAINING`.

The second one corrects R18, and the correction has the same shape as two earlier ones. That
pattern is the most useful thing in this file.

## Finding 1 — the same error, three times, by me

| Report | What I concluded | What it took to break it |
|---|---|---|
| R15 audit | The genome-wide null was "an artefact of our own supervision choice" | R17: re-running the scan under methylation supervision (83,369 vs 90,137) |
| R18 | Methylation leakage was "a property of Phikon-v2 features" | R21: three more histology encoders, all inflating 6/6 |
| R19 | The external signature was "real, **total**, and unchanged" | R20: dinov2-large reaching 0.877, not 1.000 |

Each time: a real measurement, a defensible-sounding generalisation one step beyond it, and the
generalisation wrong. Each time the limitation that would have caught it was *already written down
in the same document*. R18 says "two encoders is not a survey" and then attributes the effect to a
vendor.

**Generalisable rule, and this is the third attempt at phrasing it.** A limitation section is not a
hedge on the summary — it is a specification of which sentences are not yet licensed. If a
limitation names a comparison you have not run, no sentence anywhere in the document may assume its
result. The fix is not more caveats; it is deleting the sentence until the arm exists.

## Finding 2 — what four histology encoders bought

Under the pre-declared rule (≥5 of 6 targets inflated), reused verbatim from R18's hashed config:

| Encoder | Corpus | Mean meth. inflation | Targets |
|---|---|---|---|
| Phikon-v2 | histology | +0.0849 | 6/6 |
| UNI | histology | +0.0606 | 6/6 |
| Virchow2 | histology | +0.0773 | 6/6 |
| H-optimus-0 | histology | +0.0173 | 6/6 |
| **dinov2-large** | **natural images** | **+0.0004** | **3/6** |

Four for four on all six targets against one at 3/6 and effectively zero. The split is clean by
corpus and it was not obtainable with two encoders, whichever two.

**H-optimus-0 is the honest wrinkle.** At +0.0173 it is a quarter of Phikon-v2's, so "histology
encoders leak methylation" covers a fivefold spread. The verdict rule keys on 6/6 targets rather
than magnitude, which is why it fires — and the report states the range rather than a single number.

## Finding 3 — reusing a pre-declared rule versus writing a new one

The verdict thresholds here were not chosen for this run. They are R18's, fixed and hashed before
R18's own arms were read, applied unchanged to new encoders. That is legitimate: same rule, new
data. Had I set a threshold after seeing five encoders, the verdict would be worth nothing.

Worth being explicit because the temptation was real. H-optimus-0's +0.0173 would fail a
magnitude-based bar of, say, 0.05 — and I could have written that bar today and reported a mixed
result. The rule that existed keys on target count, so it fires, and the spread goes in the
limitations where it belongs.

## Finding 4 — the pattern I refused to claim

Relative inflation falls as site-disjoint capability rises: Spearman ρ = −0.80 across the five
encoders. Panel C makes it look like a law.

It is not reported as a finding. It was noticed *after* reading all five arms, and at n = 5 no rank
correlation can reach p < 0.05 by this test — the observed p is 0.104. Reporting it as a result
would be exactly the garden-of-forking-paths move Paper 2 is about, committed inside the report
that documents the discipline.

It is in `results.json` under `post_hoc_observation` with `STATUS` set to
"POST-HOC OBSERVATION, NOT A TESTED FINDING", the figure panel says so in red, and the text says
what it licenses: a pre-registration on encoders not used here.

**This is the most valuable thing in the run and I am not allowed to conclude it.** That is
uncomfortable and it is the correct outcome.

## Finding 5 — subtype leakage is universal but its magnitude is not transferable

All five clear the 0.05 relative bar, so the phenomenon is real everywhere. But relative inflation
runs 0.444 to 0.852 and site-disjoint AUROC runs 0.7356 to 0.9497. Virchow2 reaches 0.9497 under
site-disjoint folds where Phikon-v2 manages 0.7992.

So Paper 1's "+0.171 from fold assignment alone" is a statement about one encoder, not a constant of
the archive. The direction generalises across five encoders; the magnitude does not generalise past
the one it was measured on. Paper 1 is updated accordingly.

## Finding 6 — KEAP1 held, unprompted, five times

R15 argued on mechanistic grounds that KEAP1 mutation has no institutional prevalence structure to
exploit, so no shortcut should be available, and measured +0.0018. Across five encoders chosen for
entirely unrelated reasons the range is +0.0018 to +0.0436 — small everywhere.

A prediction made for a reason, surviving on models selected for other reasons, is worth more than
the same number repeated. It is also the only claim in this series that has never needed correcting.

## Limitations to state to the colleague

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
