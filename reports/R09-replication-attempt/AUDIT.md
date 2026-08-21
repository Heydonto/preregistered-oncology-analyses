# R09 audit — the replication attempt

Audited 7 August 2026. Four runs, all retained:

| Run | Purpose | Outcome |
|---|---|---|
| `20260807T131713Z-R09tcgarepl` | TCGA, first gate structure | HALT — controls failed (full cohort only) |
| `20260807T131826Z-R09tcgarepl` | TCGA, gate scoped to both strata | HALT — controls failed in both strata |
| `20260807T132209Z-R10matchedrepl` | first matched candidate | HALT — setting gate rejected it (median OS 93.8mo) |
| `20260807T132613Z-R10matchedrepl` | two screened cohorts, pooled | HALT — controls failed, STK11 inverted |

Config hashes are in each run's `config.yaml`; the trace for the two that the report draws on is
copied to `evidence/tcga_*` and `evidence/matched_*`.

## Verdict

The report's conclusion is the correct reading of four halted runs: **R06's finding cannot be
adjudicated with public data.** Every halt was a pre-registered gate doing its job, and none was
waived. One gate's *scope* was corrected mid-report; that change is documented below and did not
alter any outcome.

## Finding 1 — a gate's scope was corrected, and it changed nothing

The TCGA control rule as first written halted on the full cohort alone. But the same config
pre-declared a late-stage subgroup *precisely because* the full cohort might be the wrong disease
setting — so halting on the full cohort blocked the analysis designed to address that failure.

The gate was rescoped to "controls must hold in the full cohort **or** in the pre-declared
late-stage subgroup", i.e. halt only if no stratum can demonstrate known biology.

Why this is not moving the goalposts:

- the late-stage subgroup was in the config **before** the run, fixed by the first run's config
  hash — it is not a post-hoc analysis invented after seeing failure
- the **replication criterion** (q<0.05 and same direction as R06) was not touched
- the rescoped gate **failed anyway**: TP53 q=0.219 and STK11 q=0.110 in the late-stage subgroup

Both runs are retained (`evidence/gates_tcga_FIRST_HALT.json`). This is the same class of
implementation-versus-intent error as R07's control gate, which suggests the lesson generalises:
when a config pre-declares a fallback stratum, the gate must be written over all strata from the
start.

## Finding 2 — the first matched candidate was rejected before any outcome was read

`luad_mskcc_2023_met_organotropism` was the obvious choice: 2,653 samples, nominally metastatic
LUAD. Its median OS is **93.8 months**, nearly twice TCGA's and more than three times R06's
28.8. It is not an advanced-disease cohort in any survival-relevant sense — plausibly because
patients must survive long enough to develop and have metastases characterised, which selects
for long survivors.

The pre-registered setting gate caught this **before any hazard ratio was computed**, so no
outcome could have influenced the decision to exclude it. Retained at
`evidence/gates_matched_CANDIDATE_REJECTED.json`.

## Finding 3 — the control failure is consistent across settings, which is itself the finding

| | Events | TP53 HR (R06: 1.25) | STK11 HR (R06: 1.63) |
|---|---|---|---|
| TCGA full | 385 | 1.08 (q=0.67) | 1.52 (q=0.098) |
| TCGA late-stage | 222 | **0.79** (q=0.22) | 1.70 (q=0.110) |
| Matched pool | 231 | 1.14 (q=0.58) | **0.85** (q=0.58) |

Two inversions of established biology (TP53 protective in TCGA late-stage; STK11 protective in
the matched pool) across 616 events total. The interpretation the report adopts — that these
cohorts cannot test the question — is better supported than "the effect is absent", because an
absent effect would not also erase TP53 and STK11.

This also licenses a caution about R06 itself, which the report states: its result may be
specific to ctDNA-detected metastatic disease at one institution rather than general to NSCLC.
R06 is correctly labelled single-cohort in the series index.

## Finding 4 — the KMT2D near-miss, and why the control gate earned its place

The pooled analysis returned KMT2D at HR 0.52, 95% CI 0.34–0.80, **q=0.014**. Without a control
gate this is publishable-looking: a novel protective association in advanced NSCLC.

It is not interpretable, for four independent reasons:

1. it arises in cohorts that cannot reproduce TP53 or STK11
2. it contradicts R06, which found KMT2D non-significant
3. it is in the protective direction, with no mechanistic rationale offered in advance
4. R08 independently showed KMT2D's apparent *methylation* phenotype was entirely histology
   confounding (20.7% pooled → 0.0% within LUAD)

Both of the strongest-looking KMT2D results generated in this series have dissolved under
scrutiny. That consistency is worth remembering if KMT2D comes up again.

## Finding 5 — the screen is auditable and computes no outcomes

`screen_replication_cohorts.py` evaluates eligibility only — OS endpoint present, mutation
profile present, ≥100 events, median OS nearer 28.8 than 50.3 months, patients disjoint from
R06's cohort. It computes **no hazard ratio**, so eligibility cannot have been steered by
results. Per-study verdicts for all 41 studies are in `evidence/screen_results.json`.

Result: 17 have no OS endpoint, 14 have too few events, 8 are the wrong setting, 2 eligible.
Both eligible cohorts were then analysed (Finding 3), so the screen was not used to stop early.

## Reproducibility checks

- Model identical to R06's (per-gene unadjusted Cox, BH over the same five genes); adding
  covariates would have invalidated the comparison, so the adjusted model is reported as
  context only and never as the replication test
- Only `_sequenced` patients used — an unprofiled patient is not wild-type
- Patient overlap with R06's cohort removed before outcomes were read, and gated (57 removed in
  the rejected candidate, 24 in `bm_nsclc_mskcc_2023`, 0 in `lung_msk_pdx`); a further gate
  confirmed 0 R06 patients survived into the pooled frame
- Pooled model stratified by cohort, so the brain-metastasis and PDX selections are never
  contrasted with each other
- Power for the controls computed **before** running, which is why the control rule is a
  disjunction rather than requiring both — otherwise it would have been a power test dressed up
  as a validity test
- All four runs retained, none deleted
- All 22 numbers in the report verified programmatically against `evidence/consolidated.json`

## Limitations to state to the colleague

1. **No institutional independence was achievable.** Both eligible cohorts were MSK, as R06's
   was. Patients are disjoint and gated as such, but that is a weaker claim and the report says
   so explicitly rather than implying more.
2. **The two eligible cohorts are heavily selected** — brain metastases and PDX-derived tumours.
   Matching R06's median OS does not make them equivalent populations, and their control failure
   may reflect selection rather than setting.
3. **This is not evidence against R06.** No cohort demonstrated the competence to test it. The
   report states the bounded claim.
4. **The screen covers cBioPortal only.** Cohorts in dbGaP, EGA or institutional holdings were
   not screened, because they are not openly accessible. A qualifying cohort may well exist
   behind controlled access — §5 of the report specifies exactly what to look for.
5. **`lung_msk_pdx` overlap could not be verified by a shared identifier scheme** beyond patient
   ID matching, which returned 0. If that study re-identifies patients differently, undetected
   overlap is possible; it would only bias toward *agreement* with R06, and none was observed.
