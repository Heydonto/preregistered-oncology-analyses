# R16 audit — a bug that all twelve gates passed through

Audited 19 August 2026 (retrospectively; the report was written 16 August and this audit was
missing). Run evidence in `evidence/`, plan SHA-256
`752b7462fce5689993dd6c05bca956ff50799794528ec3f4abf9048cb74e6770`, exposure table frozen at
`854378633ea6992c…`. Peer-review audit copy.

## Verdict

`CONFIRMED` stands. HR 1.2023 per doubling of post-treatment peak SUV, permutation p = 0.0005, 166
analysed patients, 125 events, all 12 gates PASS.

But the most important thing in this run is that **an error survived every gate**, and the gate that
now catches it was written afterwards. That is the finding, and it is the one Paper 2 leans on.

## Finding 1 — the verdict was briefly computed from a control's hazard ratio

An earlier execution reported the primary as **HR 1.382 with 95% CI 1.050–1.124** — a point estimate
outside its own interval, which is arithmetically impossible.

Cause: inside the positive-control blocks I reused `hr` and `ci` as scratch variable names. PC2's
values overwrote the primary's before the verdict logic read them. So for one execution the
mechanical verdict was evaluated on **performance status**, a control, rather than on post-treatment
SUV, the exposure the whole plan was about.

**All eight gates in force at the time passed while this was true.** They checked recombination, join
integrity, LE uniqueness, cutoff-row inflation, exposure adequacy, the Stage-1 freeze, exclusion
accounting and permutation calibration. Not one of them looked at whether the number being reported
was the number that had been computed.

What caught it was **an internally impossible value**, noticed by eye. Not a test.

Gate `G9_primary_coherence` was added afterwards and asserts the point estimate lies inside its own
confidence interval. It is now one of the 12 that pass.

**Generalisable rule:** pre-registration constrains *analytic* discretion and does nothing for
*implementation* correctness. They are different failure modes needing different defences — hashes
and frozen plans for the first, internal-consistency assertions for the second. A gate set that
never checks its own outputs for arithmetic possibility is incomplete no matter how many integrity
checks it contains.

## Finding 2 — two coordinator errors that were mine, not the plan's

Recorded because in both cases the plan was ahead of my executor.

**LE uniqueness.** I wrote it as a hard stop. The plan had already anticipated that form LE might not
be one row per patient and specified a duplicate-resolution rule. LE turned out to have 244 rows over
243 patients; the plan's rule handled it and named the single affected patient. My executor would
have halted on a case the plan had solved.

**Field-name case.** The landmark date came back empty because the CSVs capitalise day-offset fields
(`TAe7d`) differently from the dictionary the plan quotes (`tae7d`). Same field, different case, zero
landmark dates, and the failure was silent rather than loud.

Neither required reinterpreting the plan. Both are the coordinator failing to implement a correct
specification, which is a different category from the plan being wrong — and it is the category the
pre-registration machinery does not touch.

## Finding 3 — what the blind generator got right that I would have missed

The plan was written from schema and dictionary alone, with no cell value from any outcome field. It
identified twelve traps in advance; four mattered:

- **Immortal time.** Exposure is measured 8–20 weeks after chemoradiotherapy, so survival from
  registration would guarantee survival to exposure. The plan set the time origin at the
  post-treatment PET and made min(T) ≥ 0 a gate.
- **Vital-status code 9** means "dead, date of death unknown". The obvious `event = (f1e2 == 2)`
  silently undercounts deaths. The plan treated 9 as an event and made a censoring sensitivity a
  *condition* of confirming — S5 duly returns HR 1.202, unchanged.
- **Cutoff forking.** Form SS carries 1,185 rows keyed by an SUV inclusion cutoff — five analyses per
  patient, a ready-made garden of forking paths. The plan locked SS to the most inclusive cutoff and
  forbade cross-cutoff comparison.
- **Undocumented uniqueness.** The recon package never stated LE was one row per patient, yet the
  primary exposure depended on it. The plan required verification and pre-specified the remedy.

I would plausibly have caught immortal time. I do not believe I would have caught code 9 before it
biased the result.

## Finding 4 — the caveat the plan forced into the open

**106 of 180 exposure values are at or below 0.5** — complete metabolic response — and were floored
by the plan's rule. So 59% of the cohort sits at one tied value and the association is carried by the
minority with residual uptake.

The plan required this count to be reported. Without that requirement it is exactly the sort of thing
that ends up in a reviewer's letter rather than in the paper. S3, the binary split at the plan's fixed
3.5 threshold, is arguably the more interpretable form of the same finding (HR 1.789).

## Finding 5 — controls that did not fully corroborate, and were not quietly dropped

| Control | Result | Outcome |
|---|---|---|
| PC0 | post-RT SUV below pre-RT in 94.4% | PASS |
| PC1 | progression by day 365 → shorter OS, HR 2.463 | PASS (p = 3.8e-5) |
| PC2 | performance status ≥1 vs 0, HR 1.382 | **PASS-WEAK** (p = 0.072) |
| PC3 | weight loss ≥5%, 26 exposed vs threshold 30 | **UNINFORMATIVE-BY-N** |

PC3 was declared uninformative by a rule written before the data were seen, not explained away
afterwards — the plan anticipated that trial eligibility screens performance status and weight loss,
so range restriction was expected and given its own outcome branch in advance.

Two of four not fully corroborating is weaker than a clean sweep, and the report says so. Under the
plan's aggregate rule it is sufficient for CONFIRMED because PC1 passed strongly and PC2 pointed the
right way.

## Finding 6 — G7 passed at exactly its threshold

Exposure coverage was required to be ≥ 180 and came in at **exactly 180**. One patient fewer and the
plan's declared fallback chain would have fired, changing the exposure source. Nothing was done
wrong; it is worth knowing that the run sat on the boundary of a pre-declared branch point.

## Limitations to state to the colleague

1. **Single trial, tabular only.** Core-lab PET assessments, not the 135.5 GiB of PET/CT imaging.
   The generalisation claim is about an unseen cohort and data shape, not an unseen imaging pipeline.
2. **The exposure is 59% floored.**
3. **Two of four controls did not fully corroborate.**
4. **G7 passed at exactly its threshold.**
5. **Not an independent replication.** The published ACRIN 6668 analysis used the same data; the
   claim is about the *process* recovering a known result blind, not about establishing the biology.
6. **The central negative result is about our own gates.** Twelve gates now pass; eight passed while
   the reported primary was a control's hazard ratio. Anyone reusing this harness should assume its
   gate set is still incomplete in ways not yet discovered.
