# R07 audit — methylation and EGFR-TKI response (GSE147377)

Audited 7 August 2026. Run dir:
`paper-draft/02-revised-protocol-open-datasets/runs/20260807T082054Z-R07tki-3f8bece`
Config hash (written before any outcome was read):
`623dff57160f76886a0684ac3e71c42a3fd4274f5a7382bb7893db8f9562fb0c`

## Verdict

The result stands as **inconclusive**. The pipeline is sound, the sample mapping is proven
correct, the negative control behaves correctly, and the cohort is too small to answer the
question. Three corrections were made; none changes that verdict, and one of them changes the
headline wording.

## Gates

| Gate | Expected | Observed | Result |
|---|---|---|---|
| G0_samples | 79 | 79 | PASS |
| G1_response_labels | PR 40 / PD 29 / SD 10 | PR 40 / PD 29 / SD 10 | PASS |
| G2_beta_columns | 79 AVG_Beta columns | 79 | PASS |
| G3_ahrr_present | cg05575921 present | present | PASS |
| G4_mapping_controls | at least one control holds | B holds at AUROC 1.000 | PASS |
| G6_negative_control | permutation null contains 0.50 | mean 0.491, CI 0.354–0.666 | PASS |

## Finding 1 — a gate fired incorrectly on the first run

The pre-registered rule is a **disjunction**: *"if neither control holds, the mapping is wrong
and no response result is reported."* The first implementation gated on control A alone and
halted when the AHRR/smoking control returned p=0.093.

That halt was a coding error, not a real failure. Fixed by evaluating both controls and
halting only if both fail — i.e. the code was changed to match the pre-registered rule; the
rule itself was **not** edited. Its wording is fixed by the config hash of the first (halted)
run, which is retained at `evidence/gates_FIRST_HALTED_RUN.json` for exactly this reason.

This is the failure mode that pre-registration is supposed to catch, working in the intended
direction: the rule was written correctly in advance, so the implementation could be checked
against it rather than the other way round.

## Finding 2 — the reported permutation p was anti-conservative

The observed statistic was the AUROC of probabilities **averaged over 10 repeats**; each null
replicate was a **single** 5-fold split. Averaging reduces variance and raises AUROC, so the
two statistics were not comparable and the mismatch biased p *toward* significance.

Recast like-for-like (mean per-repeat AUROC 0.532 against the same null):

- reported: p = 0.188
- corrected: **p = 0.297**

The correction moves away from significance, so the conclusion is unaffected. The report quotes
0.30 and states the reason. Recorded in `evidence/audit_addendum.json`.

## Finding 3 — the headline claim had to be weakened

The run's `_decision` field reads *"no evidence that methylation predicts TKI response at this
sample size."* Literally true, but readable as a negative finding, which the data do not
support:

- AUROC 0.555, 95% CI **0.42–0.69** — consistent with no effect *and* with a moderate one
- minimum detectable AUROC at 80% power: **0.685**

Anything weaker than 0.685 was invisible by construction. The report therefore states the
bounded claim — *if such a signal exists and is of moderate strength or less, this cohort
could not have seen it* — rather than a negative result.

What makes this readable as genuinely uninformative rather than broken is that the negative
control passed: the permutation null is centred at 0.491, so there is no leakage inflating or
deflating the estimate.

**Precedent note:** this is the third report in the series (with R04 and R06) where the audit
had to overrule the run artefact's `_decision` wording. `_decision` fields in run artefacts are
superseded by the reports and audits wherever they disagree.

## Finding 4 — the sample mapping is proven, not assumed (no correction needed)

Worth recording because it was the largest structural risk in this report. The beta matrix
columns are `EXP_1…EXP_79`, not GSM accessions, so the phenotype join assumed column order. A
silent mis-mapping would have invalidated every number while still producing plausible output.

Control B settles it: probe `cg04744025` separates 46 female from 33 male samples at
AUROC = 1.000. Chance probability for one probe is 1/C(79,33) ≈ 5×10⁻²³; across the 20,000
probes screened, the expected number of chance-perfect separations is ~10⁻¹⁸. The mapping is
correct.

Control A (AHRR) is consistent but underpowered: direction correct (0.641 in smokers vs 0.702
in never-smokers, the established direction), p=0.093, with only 24 smokers among 79 patients
and no current/former or intensity distinction. A power limitation of the control, not evidence
against the mapping.

## Reproducibility checks

- Probe selection and hyperparameter search occur strictly inside training folds
  (`SelectKBest` and `GridSearchCV` are pipeline steps, refit per fold) — no selection on the
  full data
- Config written and hashed before any outcome was read
- Both runs retained, including the halted one
- `MANIFEST.sha256` covers every artefact
- All 13 numbers in the report verified programmatically against `evidence/*.json`
- Figure regenerable from archived evidence via `make_figure.py`

## Limitations to state to the colleague

1. **Underpowered by design, and quantified.** 69 patients against ~465k probes. Detectable
   only at AUROC ≥ 0.685. To reach AUROC 0.60 at 80% power needs 263 patients.
2. **SD patients excluded** (10 of 79), pre-declared. Response is dichotomised PR vs PD, so
   the analysis speaks to the extremes of the response distribution, not to a gradient.
3. **Single cohort, single platform.** No replication cohort exists in our holdings; no batch
   or array-position covariates are available in the series metadata, so technical confounding
   cannot be adjusted for — though it would tend to *create* apparent signal, not remove it,
   which makes the null more credible rather than less.
4. **Smoking control unresolved.** Left as-is rather than substituted for a better one, since
   selecting a control after seeing it fail would defeat its purpose.
