# R18 audit — the answer strengthens one claim and narrows another

Audited 17 August 2026. Run `20260817T060045Z-R18encodercmp-5002607`, config
`143b3b9d0a92fabbef0a926121b7a75a1eeacc003bdbb8b0b690e884757b1944`. Peer-review audit copy.

## Verdict

All six gates passed. Two findings, pointing in opposite directions:

- **`SITE_SIGNATURE_IS_ARCHIVE_NOT_ENCODER`** — subtype leakage reproduces under an encoder that
  has never seen a histology slide, and is **larger** there (+0.2034 against +0.1710). Paper 1's
  headline claim survives its own first limitation.
- **`METHYLATION_LEAKAGE_DID_NOT_REPRODUCE`** — under dinov2-large the +0.0849 mean methylation
  inflation becomes +0.0004, with three of six targets moving in the anti-leakage direction.
  *(**Corrected 2026-08-18 by R21.** This bullet originally said the effect was "**specific to
  Phikon-v2**". Wrong. R21 added UNI, H-optimus-0 and Virchow2: all four histology encoders inflate
  all six targets, +0.0173 to +0.0849, and only the natural-image encoder does not. The
  dissociation is by training corpus, not by vendor. See Finding 7.)*

Both verdict strings, and the quantity deciding them, were fixed in the config and hashed before
either dinov2 arm was read.

## Finding 1 — the reproduction gate was the whole ballgame, and it was exact

A cross-encoder comparison is worthless if the harness is not the one that produced the published
numbers. The unified harness was written fresh, so the Phikon-v2 arms were re-run through it and
gated against R15 first:

| Control | Recomputed | Published |
|---|---|---|
| subtype AUROC, site-grouped | 0.7992 | 0.799 |
| subtype AUROC, random folds | **0.9703** | 0.9703 |
| KEAP1 AUROC, site-grouped | 0.6642 | 0.664 |
| mean methylation inflation | **0.0849** | 0.0849 |

All twelve methylation correlations — six grouped, six random — came back **bit-identical** to the
stored values. That is what deterministic seeding should give, and it matters because it means the
dinov2 differences are not absorbing harness noise. Had any of G2–G5 missed, the run would have
halted and the comparison would have been void.

## Finding 2 — the metric that decides was chosen before the results, and it had to be

dinov2-large is not trained on histology, so its absolute performance is lower (subtype 0.7356
site-grouped against 0.7992). That creates an obvious trap: a smaller raw inflation could mean
*less leakage* or merely *less signal to inflate*. Whichever I picked after seeing the numbers, I
could have told a clean story.

So relative inflation — (random − grouped) / (1 − grouped), the fraction of remaining headroom that
fold assignment alone recovers — was declared in the config as the deciding quantity before any
dinov2 arm ran. It is the comparison that survives one encoder being weaker.

As it happens the raw figure was *more* favourable to the conclusion than the relative one
(dinov2's raw inflation exceeds Phikon-v2's; its relative inflation is slightly lower, 0.7694
against 0.8518). Pre-declaring cost nothing here. It would have cost something if the numbers had
landed the other way, which is the point.

## Finding 3 — the dissociation is not a headroom artefact, and I checked

The tempting dismissal of the methylation result is that dinov2's methylation signal is weaker
(mean grouped ρ 0.1658 against 0.2856), so of course its inflation is smaller.

That does not survive the relative figure. On fraction-of-headroom, methylation inflation is
**0.119 for Phikon-v2 and 0.000 for dinov2-large**. A weaker representation with less room to
inflate would show a *reduced* gap. It would not show a gap that vanishes while the same
architecture's subtype gap is the larger of the two. Three of six targets are negative
(KEAP1 signature −0.0157, global −0.0145, gene body −0.0144), which is what noise around zero looks
like, not attenuated leakage.

## Finding 4 — what one arm alone would have told me

Had I run only the subtype arm — the cheaper, headline-relevant one — I would have written "the
site signature reproduces under a second encoder, so it is a property of the archive" and stopped.
That sentence is true of subtype and false of methylation, and Paper 1 currently states a single
site-leakage magnitude covering both.

**Generalisable rule:** when a paper reports one effect measured on two endpoints, a robustness
check has to cover both endpoints or it licenses only half the claim. The temptation is to test the
headline and generalise, which is the same shape of error as R15's audit generalising six aggregate
targets to the methylome.

## Finding 5 — a mechanism I can offer but have not tested

Phikon-v2 is trained on histology and encodes fine-grained stain and preparation character.
Methylation arrays for a given TCGA patient were typically run at the institution that cut the
slide, so batch structure in the assay and batch structure in the image share an index. A
histology-specialised encoder can pick that up; a natural-image encoder captures coarser
morphology — enough for subtype under site-disjoint folds, apparently not enough to carry array
batch.

**This is a hypothesis and the report labels it as one.** Nothing in R18 tests it. Testing it would
mean relating per-site methylation batch identifiers to feature-space site separability, which is a
different experiment.

## Finding 6 — KEAP1 behaved as R15 predicted, unprompted

R15 argued that KEAP1 mutation status has no institutional prevalence structure to exploit, so
there should be no shortcut available, and reported inflation of +0.0018. Under dinov2-large it is
+0.0202 — still small, and small under an encoder chosen for a different reason. A prediction made
for a mechanistic reason holding up on a second encoder is weak evidence, but it is evidence, and it
was not something the run was set up to check.

## Finding 7 — the attribution was wrong, and my own limitation said why (R21, 2026-08-18)

This report concluded methylation leakage was "a property of Phikon-v2 features" from a comparison
of exactly two encoders, one of which was not a pathology model. Limitation 1 below reads "Two
encoders is not a survey." Both sentences are in this document.

R21 ran three pathology foundation models:

| Encoder | Corpus | Mean meth. inflation | Targets |
|---|---|---|---|
| Phikon-v2 | histology | +0.0849 | 6/6 |
| UNI | histology | +0.0606 | 6/6 |
| Virchow2 | histology | +0.0773 | 6/6 |
| H-optimus-0 | histology | +0.0173 | 6/6 |
| dinov2-large | natural images | +0.0004 | 3/6 |

The effect belongs to **histology pretraining**, not to a vendor — a stronger claim than the one I
made, reached by running the arm I had already identified as missing.

Finding 4 of this audit is also superseded in part: it read the 512-d-vs-1-d gap as partly a
variance effect, which stands, but its framing of methylation leakage as encoder-idiosyncratic does
not.

## Limitations to state to the colleague

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
