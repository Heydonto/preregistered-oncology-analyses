# R19 audit — a manuscript number with no artefact behind it

Audited 17 August 2026. Run `20260817T010512Z-R19extsiteprobe-a6b6353`, config
`5985b7c75d064a204d4d512e8cd3fb7b6389768c596f07d0b4d719788ea92039`. Peer-review audit copy.

## Verdict

The science survives: the external site signature is real and unchanged by restricting to H&E. All
six gates passed. *(Amended 2026-08-17: this line originally said "real, **total**, and unchanged".
R20 shows the totality is a Phikon-v2 property — dinov2-large reaches 0.8768 on the same slides — so
"total" is withdrawn. See Finding 7.)*

The process did not. Two of Paper 1's rows were in the manuscript without evidence, and this
programme had described Paper 1 as "23/23 claims verified".

## Finding 1 — "verified" meant "recomputed from a number I remembered"

Paper 1's probe table has three rows. One is backed by `results/R12-federated-simulation/evidence/
results.json`. The other two — 1.0000 and 0.9957 — were produced in a working session, printed,
and typed into the manuscript. No JSON, no gates, no config hash, no run directory.

Every other quantitative claim in this programme is checked against a file in `evidence/`. These
two were checked against my own recollection, and the verification pass recorded them as verified
because the numbers in the manuscript matched the numbers I believed. That is not verification, it
is agreement with myself.

**What made it invisible:** the claim-checking pass compared manuscript text against reported
values. It never asked *whether a file existed*. A missing artefact and a matching artefact both
read as "consistent".

**Generalisable rule:** a verification pass must resolve every number to a path, and the absence
of a path is a failure, not a neutral outcome. Matching is the second check; existence is the
first. This is now the rule for the series.

## Finding 2 — the file was named for the filter it did not apply

The probe input is `ext_feat/external_he_mean.npz`. Three of its 29 HTMCP slides are not H&E:

| Slide | Stain |
|---|---|
| `HTMCP-02-06-01032-CHR` | chromogenic ISH |
| `HTMCP-02-15-01090-CHR` | chromogenic ISH |
| `HTMCP-02-07-01065-06-P16` | p16 immunohistochemistry |

R12's audit had already identified this risk in writing — HTMCP-LC is majority IHC, only 30 of its
84 slides are H&E, and pooling them would measure a *stain* signature while calling it an
*institutional* one. The mitigation was named, the file was named after the mitigation, and the
mitigation was not applied.

**The naming is the dangerous part.** An unfiltered file called `external_mean.npz` invites the
question. One called `external_he_mean.npz` answers it, wrongly, and nobody looks again. A filename
is not a gate, and here it actively suppressed the check it appeared to document.

Gate G1 in this run is written *inversely* for that reason: it asserts the contamination is
present and names the three slides. If someone silently repairs the input file, that gate fails and
forces a re-read of this report rather than quietly passing.

## Finding 3 — the concern was legitimate and the result is clean

This is the part that could have gone either way, and it went well.

| Arm | n | Bal. acc. | Null mean | Null p95 | p |
|---|---|---|---|---|---|
| published 78 slides | 78 | 1.0000 | 0.5074 | 0.6221 | 0.0050 |
| **H&E only** | 75 | **1.0000** | 0.5058 | 0.6174 | 0.0050 |

**Difference attributable to the three non-H&E slides: 0.0000.** Removing them removes nothing,
because the remaining 75 slides already separate perfectly. The stain shortcut existed and was not
needed.

Reported because a null result from a control is still a result: had I only written prose saying
"three slides cannot matter much", the claim would rest on my judgement. It now rests on an arm.

## Finding 4 — one row does not reproduce, and cannot be made to

Row 3 recomputes to **0.9816** against the published **0.9957**. The gap is 0.014 and its cause is
unrecoverable: the original probe's model, regularisation, cross-validation scheme and number of
subsampling draws were never written down.

I cannot say the published number was wrong. I can only say it is not reproducible from the
surviving inputs, which for a document held as an IP record is the same problem wearing a politer
face. Paper 1 is corrected to the reproducible value with its method stated.

## Finding 5 — a p-value I declined to quote

The 3-way arms return p = 0.0476, which is exactly 1/21 — the floor for 20 subsampling draws. It
is a property of my loop count, not of the data, and quoting it would imply weak evidence for a
0.98-against-0.33 separation.

Paper 1 currently prints p = 0.0050 for this row, from a permutation scheme that is not recorded
and cannot be reconstructed. Neither number should be shown. The report gives the separation, the
null means (0.3345 and 0.3660) and the across-draw standard deviations (0.013 and 0.019), and says
why there is no p.

**The temptation worth naming:** 0.0476 is under 0.05, so it would have passed unremarked. A
threshold-crossing number produced by an arbitrary implementation choice is the easiest kind of
false precision to ship.

## Finding 6 — what a perfect score does and does not license

1.000 on 75 slides with 1024 features is separation, not measurement. The permutation null is what
makes it meaningful — 0.5058, p95 0.6174 — and the report leads with that rather than the 1.000.
The substantive claim is that the separation is *total*, not that it exists; two cohorts differing
in fixation, stain lot, scanner and protocol should be distinguishable, and a representation with
no invariance objective will oblige.

## Finding 7 — a third defect in the same file, found later (R20, 2026-08-17)

`ext_feat/external_he_mean.npz` held **29 of the archive's 80 HTMCP slides**. 51 were absent, the
input having been assembled from whatever was staged locally, and neither this report nor Paper 1
recorded it.

Phikon-v2 returns **1.0000 on the full 76-slide H&E set and on all 129** as well, so the omission did
not bias the headline. The defect is provenance, not correctness — but it completes a pattern in this
one file: no persisted artefact (Finding 1 of R19's summary), a name asserting an H&E filter it did
not apply (Finding 2 above), and now an undocumented subset. None was catchable by reading the
analysis code, which faithfully processed whatever it was given. All three needed going back to the
archive and counting.

**Rule now applied:** any derived aggregate must carry a manifest of contributing IDs inside the
artefact. See `spine/stage_dinov2_volume.py::mean_embeddings`, which writes `slide_ids`, `cohort` and
`is_he` alongside the matrix. Full account in `results/R20-external-probe-two-encoders/`.

## Limitations to state to the colleague

1. **Stain came from filename suffixes** (`HE` / `CHR` / `P16`), the only stain field in the staged
   features. Source mislabelling would propagate.
2. **75 slides, two institutions.** Bounds existence, not effect size.
3. **Slide-level mean features**, not the attention-MIL representation.
4. **Row 3's 0.014 discrepancy is permanently unattributable.**
5. **Phikon-v2 only** — answered by R20: the signature holds under a second encoder (0.8768,
   p = 0.0050) but not perfectly, so "total separation" is encoder-specific.
6. **Two Paper 1 rows were unevidenced until today**, and the "23/23 verified" statement covering
   them was wrong. Corrected in Paper 1 in place.
