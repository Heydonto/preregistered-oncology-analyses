# What a reviewer can reproduce

No bulk data is redistributed here. Most of it we are not permitted to redistribute, and the rest
is large enough that pointing at the source is better than copying it. This file states, per
report, whether the underlying data is publicly obtainable.

## Fully reproducible from public data

| Report | Data | Route |
|---|---|---|
| R08, R10, R15, R17 | TCGA-LUAD/LUSC 450k methylation, clinical | GDC open tier |
| R12, R13, R15, R18, R20, R21, R22 | TCGA diagnostic whole-slide images | GDC open tier |
| R16 | ACRIN 6668 clinical | TCIA, CC BY 3.0, **no application required** |
| R23 | GSE68379 450k + GDSC2 dose–response | GEO + cancerrxgene.org |
| R24 | GSE202859 RNA TPMs | GEO |
| R01–R05 | UPENN-GBM, TCIA collections | TCIA, CC BY 4.0 / CC0 |

## Reproducible only with the model weights

R18 and R21–R22 use five tile encoders. Three are freely downloadable
(`owkin/phikon-v2`, `facebook/dinov2-large`, `bioptimus/H-optimus-0` — the last Apache-2.0).
**`MahmoodLab/UNI` and `paige-ai/Virchow2` are CC-BY-NC-ND-4.0**, gated behind terms acceptance,
and we cannot redistribute their weights or the features derived from them. A reviewer must accept
those terms to reproduce those two arms. The other three arms stand alone.

## Not publicly obtainable

| Report | Blocker |
|---|---|
| R14 | Vanguri multimodal ICI cohort — controlled access |
| R06, R09 | MSK cohorts via cBioPortal; terms vary by study, one carries a licence caveat |

## Data integrity warning

**Verify member counts, not just checksums.** Our mirror of `GSE68379_RAW.tar` is byte-complete
against both our storage and GEO's declared size (9,020,579,840), ends with a valid tar terminator,
and raises no error on read — while containing **1,302 of 2,056 members**, 649 of 1,028 cell lines.
Two independent extraction modes gave the identical truncation.

It was caught only because the job asserted an expected array count and halted. R23 therefore
fetches its 396 IDATs per sample from GEO rather than from the archive.

`data/MANIFESTS/` carries expected member counts alongside checksums for this reason. We recommend
anyone reproducing this work do the same, and we assume other bulk archives in our inventory may
have the same defect until counted.

## Environment

`environment/versions.json` and `requirements.txt` record the interpreter and package versions.
Per-run environments are captured verbatim in `reports/*/evidence/env.txt` — these are `pip freeze`
dumps from a shared workstation and therefore include packages unrelated to these analyses. They
are left unedited because they are provenance records.
