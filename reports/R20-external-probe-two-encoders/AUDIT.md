# R20 audit — a third defect in the same file, and a word Paper 1 has to give up

Audited 17 August 2026. Run `20260817T061230Z-R20extprobe2enc-977a613`, config
`ab01b9847dedcfcf3db9cffc98a66a29a04c2ece09f2cf27937b4a42aab0f0e2`. Held for IP.

## Verdict

All five gates passed. Verdict `EXTERNAL_SIGNATURE_WEAKER_UNDER_SECOND_ENCODER`, pre-declared.

The external site signature is real under both encoders — 0.8768 against a null of 0.4990 at
p = 0.0050 for a model that has never seen a pathology slide. But *perfect* separation is a
Phikon-v2 property, not a fact about the archive.

## Finding 1 — the third defect in one intermediate file

`ext_feat/external_he_mean.npz` has now yielded three problems, found on three separate occasions:

| # | Defect | Found by |
|---|---|---|
| 1 | The two probe rows it fed had **no persisted artefact at all** | R19 |
| 2 | Named `_he_` while containing **3 non-H&E slides** | R19 |
| 3 | Held **29 of the archive's 80 HTMCP slides** — 51 absent, undocumented | R20 |

All three are one species: *an intermediate file whose name asserted a property nobody re-checked.*
Not one was caught by reading the analysis code, because the analysis code was correct — it faithfully
processed whatever the file contained. They were caught by going back to the archive and counting.

**Generalisable rule:** a derived input needs a provenance record stating what it was built from and
what was excluded, and that record has to be checked against the source, not against the filename.
For this programme: any `*_mean.npz` or similar aggregate must ship a manifest of contributing IDs.
The new builder (`spine/stage_dinov2_volume.py::mean_embeddings`) writes `slide_ids`, `cohort` and
`is_he` into the artefact for exactly this reason.

Gate G2 asserts the omission is present and that the count is exactly 51, so a silent repair of the
input fails rather than passing quietly.

## Finding 2 — and the defect turned out not to matter, which I could not have known in advance

Phikon-v2 returns **1.0000 on all three sets**: R19's 78, the full 76-slide H&E set, and all 129.
The undocumented subset did not bias R19's headline at all.

That is worth stating plainly and it is also worth not over-reading. The subset was *not* chosen to
flatter the result — it was chosen by which files happened to be on local disk — and it happened to
be harmless. Had HTMCP's omitted 51 slides included the visually atypical ones, the story would be
different. Reporting "no harm done" is honest; treating it as evidence the practice is safe is not.

## Finding 3 — the word Paper 1 has to give up

Paper 1's abstract says two cohorts separate *perfectly*, and its discussion leans on that:
"Two genuinely independent cohorts separate perfectly."

On the same 76 slides with the same probe: Phikon-v2 1.0000, dinov2-large 0.8768. The *existence* of
the signature is encoder-independent and strongly significant under both. The *totality* is not.
Paper 1 is corrected to attach an encoder to the claim.

This is the same shape as R18's methylation finding, and the two together now form a pattern worth
naming: **a histology-specialised encoder registers institutional character more completely than a
generic one.** Subtype leakage was *larger* under dinov2 (+0.2034 vs +0.1710), methylation leakage
vanished (+0.0004 vs +0.0849), and external separation dropped (0.877 vs 1.000). Those three do not
line up in one direction, so the pattern is real but not simple, and no single sentence covers it.

## Finding 4 — SET_ALL behaved exactly as R12's audit predicted

dinov2-large scores **higher** on all 129 slides (0.9444) than on H&E alone (0.8768). Adding 53 IHC
slides makes the task *easier*. That is the stain shortcut, visible in the numbers rather than
argued for.

R12's audit had flagged in writing that HTMCP-LC is majority IHC and that pooling would measure a
stain signature while calling it institutional. R19 then partly reproduced the problem anyway. Here
SET_ALL is included deliberately and labelled as an upper bound, not a measurement — the useful form
of a known trap is to run it and show the gap.

## Finding 5 — what a pre-declared primary arm bought

The config named `dinov2-large / SET_HE` as the primary arm before any number was read. With six arms
available and values from 0.877 to 1.000, choosing afterwards would have let me report either
"the signature is encoder-independent" (SET_ALL, 0.944) or "it is much weaker" (SET_HE, 0.877)
truthfully.

The pre-declaration also fixed the verdict thresholds at 0.95 and 0.75, so 0.8768 landing in the
middle band produced `WEAKER_UNDER_SECOND_ENCODER` mechanically rather than by my judgement about
whether 0.877 counts as "strong".

## Finding 6 — every number here was one CV split (R22, 2026-08-19)

R22 re-measured with 25 splits instead of one:

| Encoder | this report | 25-split mean ± sd |
|---|---|---|
| Phikon-v2 | 1.0000 | 1.0000 ± 0.0000 |
| dinov2-large | 0.8768 | **0.9176 ± 0.0247** |

The 0.8768 was the pre-declared seed, so not cherry-picked — but it is 1.6 sd below its own mean, and
this report quoted it as a fixed quantity. The honest gap is 1.000 vs 0.918. Verdict and conclusion
unaffected; the headline contrast was overstated.

A cross-validated metric on 76 samples is a random variable. Quoting one draw of it is a category
error, and it is a different failure from the three "generalised past the arm I ran" errors elsewhere
in this series. See `results/R22-external-probe-five-encoders/`.

## Limitations to state to the colleague

1. **Two institutions, 76 H&E slides.** Bounds existence, not effect size.
2. **Slide-level mean features**, not the attention-MIL representation.
3. **Stain from filename suffixes** (`HE` / `CHR` / `P16`) — source mislabelling propagates.
4. **1.000 and 0.877 are not strictly commensurable.** dinov2-large is weaker at every histology
   task here, so a lower separation is expected in direction. What the arms establish is that both
   clear a permutation null, not that the gap between them is a calibrated quantity.
5. **Two encoders is not a survey** — same gating limitation as R18.
6. **Paper 1's "perfectly" and R19's slide-set description are corrected in place.**
