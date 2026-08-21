# R17 audit — the genome-wide null, re-examined

Audited 17 August 2026. Run `20260817T002058Z-R17percpgmeth-6626d09`, config
`4779152747f14818e8f1279733a0fe2dd9209dfcbef8eb8662890f6f5b016a42`. Peer-review audit copy. This audit was written as an internal adversarial check before the work
was submitted for review.

## Verdict

The headline stands and it is a correction to our own prior report. **R15's genome-wide negative
is not an artefact of subtype supervision.** All 13 gates passed, including four that had to pass
before the new arm was permitted to count.

The uncomfortable part is not the result. It is that R15's audit — the document whose whole
purpose was to catch over-reach — was the place the over-reach happened.

## Finding 1 — this run exists because the previous audit generalised

R15's Finding 2 was titled "the first negative was an artefact of our own supervision choice" and
called it "the single most consequential correction in this report". It was built on a real
observation: six aggregate methylation targets went from nothing to partial ρ 0.135–0.244 when
supervision changed. From that it derived a general rule — *"a negative result from a probe
supervised on target A is not evidence about target B"* — and applied it to the genome-wide scan.

The rule is sound. Its application here was not tested, and now that it has been, it is wrong at
CpG resolution:

| Arm (512-d, matched folds/ridge) | above null | % | median ρ | max ρ |
|---|---|---|---|---|
| methylation-supervised | 83,369 | 20.9 | 0.0073 | 0.350 |
| subtype-supervised | 90,137 | 22.6 | 0.0057 | 0.485 |

Not better. Marginally worse, on every one of the three summaries.

**What R15 got right, and it matters:** its limitations section states plainly that "per-CpG
prediction under methylation supervision was not attempted." The gap was known and disclosed. The
failure was not concealment — it was that a disclosed gap sat next to a narrative sentence
("reversed the answer", "recovered the signal") which read as though the gap had been closed.

**Generalisable rule, and I would rather have this one than the last one:** stating a limitation
does not bound a claim. If a sentence in the summary is broader than the experiment, the
limitation at the back will not fix it — a reader takes the summary. The remedy is to narrow the
sentence, not to add a caveat.

## Finding 2 — the reimplementation was gated before it was trusted

A correction to an earlier report is worthless if the new code is merely different. Four gates ran
first:

| Gate | Check | Observed |
|---|---|---|
| G4 | reproduces R15's per-CpG vectors, r > 0.999 | **r = 1.000000** both arms |
| G5 | vectorised Spearman vs `scipy.stats.spearmanr` | max dev **2.1e-08** |
| G8 | R15's published counts on R15's own null threshold | **110,212 / 220,925**, exact |
| G1 | `fold_of` equals GroupKFold(5) by site, identical across arms | verified |

G8 is the one that matters most: the two published numbers came back to the unit. So the new
arm's null cannot be attributed to a different pipeline.

G1 deserves a note. Rather than re-deriving folds, this run reads `fold_of` as written by the MIL,
so the ridge folds *are* the representation's folds and nest exactly. Re-deriving them would have
been defensible only after checking they agree — so the check is the gate, not an assumption.

## Finding 3 — the control I nearly failed to include

The methylation-supervised embedding was trained on six aggregate targets, one of which is the
mean over every probe in the scan. My first instinct was to treat that as a fatal circularity and
caveat it in prose.

It is the opposite of fatal, and prose was the wrong instrument. The circularity **favours** the
arm under test — it should have made the methylation-supervised embedding win — and it still lost.
That makes the negative stronger, not weaker.

But writing that down is not evidence, so arm `A_globalmean` was added: predict every CpG from the
patient's *observed* global mean methylation, i.e. ground truth for the dominant axis. Result:

| Predictor | above null | % |
|---|---|---|
| observed global mean methylation (1-d) | **298,320** | **74.7** |
| subtype label (1-d) | 206,412 | 51.7 |
| either 512-d embedding | ~83–90k | ~21–23 |

**Most of what "per-CpG predictability" measures is one scalar axis.** That single number
reconciles R15's two results without contradiction: five of R15's six aggregate targets are means
over large annotation classes and the sixth is the global mean, so a representation can learn that
axis — earning partial ρ 0.135–0.244 — while carrying nothing at CpG resolution. Verdict
`GLOBAL_AXIS_ONLY`.

Without this arm the report would have said "methylation supervision does not help" and left the
reader with no explanation. With it, there is a mechanism.

## Finding 4 — the 512-d vs 1-d gap is confounded, and the report says so

Both scalars beat both embeddings, by a factor of two to three. It is tempting to present that as
a finding about morphology. It is not clean: ridge with 512 correlated dimensions on ~608 training
patients has far higher variance per CpG than one well-chosen scalar, so an unknown share of the
gap is regularisation.

R15 hit the same pattern — "a 512-dimensional embedding lost to one binary variable" — and read it
as evidence the representation was wrong for the question. That reading is not supported, because
the *right* representation loses to the same binary variable by a similar margin.

This is why the report's conclusion rests only on `A_meth` vs `A_sub`, where dimensionality,
architecture, folds, ridge and alpha are all held fixed and supervision is the single difference.

## Finding 5 — genomic context, and not overselling 0.212

The methylation-supervised arm has a context spread of 0.212 against the subtype-supervised arm's
0.047. It would be easy to call that partial recovery of R10's structure.

It is not. R10 had islands **depleted** at 0.54 and open sea **enriched** at 1.44. This arm has
islands **enriched** at 1.116 and open sea **depleted** at 0.918 — the same axis, opposite sign.
A larger spread pointing the wrong way is not a weak version of the right answer, and the report
states it that way.

## Finding 6 — a small defect inherited from R15, fixed and reported

R15's permutation null drew `rng.permutation` twice, so it scored `Y[p1]` against a model fitted
on `Y[p2]`. Both are decoupled from the embedding, so the null is valid — but it is not the
single-permutation null the code reads as, and nobody noticed.

Fixed here (one permutation) and the old variant recomputed alongside it: p95 0.05763 single
versus 0.05953 two-permutation, medians −0.0071 and +0.0021. The effect is immaterial and the
ordering of the arms is unchanged. Recorded because the next person to read that line should not
have to rediscover it.

## Limitations to state to the colleague

1. **This does not overturn R15's headline.** Aggregate partial ρ 0.135–0.244 beyond subtype,
   surviving purity adjustment, stands. What is bounded is the genome-wide per-CpG claim.
2. **One architecture, one alpha.** A different readout or per-CpG tuning is not excluded, but it
   would need pre-declaration to count.
3. **Site-grouped folds only.** No random-fold variant; the leakage question was settled in R15.
4. **The negative is about resolution, not existence.** Morphology tracks a global methylation
   axis weakly. The claim it does not support is CpG-level prediction.
5. **Two sentences in R15's report and one finding in R15's audit are now wrong** and have been
   corrected in place, with the original text quoted so the change is visible.
