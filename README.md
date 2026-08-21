# Pre-registered oncology analyses

Audit copy prepared for peer review. Twenty-three analyses of public oncology data, each with a
pre-declared configuration, mechanical gates, a written record of what the gates returned, and an
adversarial audit written against the analysis after the fact.

Two manuscripts draw on this work:

- **Paper 1** — institutional confounding in whole-slide models, measured across five encoders
- **Paper 2** — pre-registration of machine-generated analyses, and where it fails

## Start here

```bash
python3 verify_claims.py
```

Every load-bearing claim in both manuscripts is pinned to **one file and one field**. The script
resolves all 84, reports `NO_FILE` / `NO_FIELD` / `MISMATCH` as distinct failures, and then re-runs
every claim with a **perturbed expectation** — if any perturbation is accepted, the suite fails. It
currently reports 51/51 for Paper 1, 33/33 for Paper 2, at 100% detection power.

That negative control exists because an earlier version of this check had **zero** power: it
searched all evidence for each number, and with 1.5M distinct numeric strings across 1,356 files
every possible four-decimal AUROC appears somewhere. A test that cannot fail is not a test. The
history is in the module docstring.

## Layout

```
reports/R01…R24/     REPORT.pdf, AUDIT.md, evidence/, make_figure.py
papers/              both manuscripts, PDF and LaTeX source
data/INVENTORY.pdf   every dataset: size, licence, access route
data/MANIFESTS/      accessions, checksums, expected member counts
environment/         pinned versions
tools/               the grounding sweep used across all reports
CLAIMS.md            the claim → file → field manifest, human-readable
CORRECTIONS.md       every correction we made to our own work, with original wording
LIMITATIONS.md       consolidated, per report
ACCESS.md            which results a reviewer can reproduce, and which need controlled data
```

## What is unusual about this repository

**Every report has an audit written against it.** The audits are adversarial by design and several
overturn the report they accompany. `CORRECTIONS.md` lists all of them with the original wording
preserved, including four cases where we generalised one step past the arm we had actually run.

**Gates halt.** Configurations are written and hashed before outcomes are read; a failed gate stops
the run and the failed run is kept. Several reports exist only because a gate fired — one halted
twice on a silently truncated input archive that was byte-complete and checksum-valid.

**Negative and confounded results are reported at the same length as positive ones.** R09 failed to
replicate in every cohort tried. R23's pre-registered primary test returned a significant answer
that its own post-hoc control reversed. R24 returned PARTIAL and then failed to replicate.

## Honest notes about these files

- `reports/*/evidence/env.txt` is a verbatim `pip freeze` from a shared workstation and therefore
  lists packages unrelated to these analyses. It is left unedited: it is a provenance record, and
  tidying one would falsify it.
- Two side-quest analyses present in our internal tree are **excluded** here. Both concern a third
  party's unpublished manuscript and figure code, which is not ours to publish.
- No bulk data is redistributed. See `ACCESS.md` and `data/MANIFESTS/`.

## Correspondence

Reza Nehzati, Ph.D. — VMC MAR COM Inc. DBA HeyDonto, Knoxville, TN, United States
