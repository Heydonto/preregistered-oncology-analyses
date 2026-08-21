# Audit of R03 — 2026-08-07

Four checks plus the six in-run gates. A null result is only meaningful if the pipeline is
demonstrably capable of finding a real association, and if the detectable effect size is
stated. Both were tested.

## Verdict: R03 stands, with one required caveat on power.

| # | Check | Finding | Verdict |
|---|---|---|---|
| C0 | **Positive control** — can the pipeline find a known-true association? | AFP gene expression vs independently measured serum AFP: **rho = 0.667, p = 5.5e-32, n = 238**. Two separate assays of the same biology agreeing strongly. | Pipeline and join verified |
| C1 | **Censored clinical inputs** | 43 of 238 serum AFP values are capped (`>1210` etc.). Stripped to the threshold value; rank-based tests preserve ordering, so the control is unaffected. | Handled, disclosed |
| C2 | **Power for the null** — could a real effect have been missed? | MVI absent vs present: observed rank-biserial **0.097**; minimum detectable at 80% power **0.212**. | **Small effects NOT excluded** |
| C3 | **Direction of point estimates** | Median FASN rises across invasion grades (262.7 → 308.7 → 326.1 TPM) in the predicted direction, but means do not (423 → 434 → 383) and the test is null. | Reported, not claimed |
| C4 | **Expression sanity** | FASN 12.9–2763.5 TPM, median 304.6, zero zeros — plausible for a highly expressed hepatic gene, non-degenerate. | Clean |

## The caveat that must travel with this result

C2 is the honest limit. With n=237 split 107/130, this cohort can exclude a **moderate or
larger** association between FASN and vascular invasion but **cannot exclude a small one**.
The report states this explicitly and does not claim proof of no effect. Anyone citing R03
as showing "FASN is unrelated to vascular invasion" would be overreading it; the correct
reading is "no moderate-or-larger association, in two independent cohorts."

## What the audit could not check

- Whether the CLCA expression matrix was correctly normalised upstream. Consumed as
  published; the AFP control gives indirect reassurance.
- Whether MVI grading was consistent across contributing centres.
- Tumour purity as a confounder. `TUMOR_PURITY` is available on all 238 and was **not**
  adjusted for — a purity-driven dilution effect could in principle mask an association.
  This is the most defensible extension if the question is pursued further.

## Relationship to the other reports

R03 completes a consistent pattern across all three reports: **the machinery finds real
associations when they exist** (survival in R02 at p=0.005, age at p=2.2e-18, MGMT survival
at p=5.8e-8, AFP concordance here at p=5.5e-32) **and finds nothing when the claimed
association is absent** (MGMT from radiomics in R01; FASN clinical associations here). The
nulls are therefore credible rather than artefactual.

---

## Addendum, 2026-08-07 — purity confounding closed

The audit above listed tumour purity as the most defensible unexamined confounder. It has
now been tested and the gap is closed.

| Check | Result |
|---|---|
| FASN vs tumour purity | Spearman **+0.202**, p=1.7e-3, n=238 (purity range 0.11–1.00, median 0.60) — a real but weak relationship |
| Does purity differ by invasion group? | Kruskal-Wallis **p=0.280** — no. Purity cannot be masking an MVI association it does not track |
| FASN residualised on purity, vs MVI (M0/M1/M2) | **p=0.583** (unadjusted 0.407) |
| FASN residualised on purity, MVI absent vs present | **p=0.319** (unadjusted 0.201) |

**Conclusion: purity adjustment does not change the null — it slightly strengthens it.** The
R03 finding is not a dilution artefact. No change to the report's numbers is required; this
addendum is the record that the confounder was examined rather than merely acknowledged.

Method: FASN log2-transformed, ordinary least squares on purity, residuals retested. n=237
(one sample lacks an MVI grade).
