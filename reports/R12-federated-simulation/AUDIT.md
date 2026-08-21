# R12 audit — federated learning across simulated hospital silos

Audited 14 August 2026. Config hash `b5200a4a9ece7481…7ac3571`. Evidence read directly from
`evidence/results.json`, `evidence/gates.json` and the 45 out-of-fold prediction files.

## Verdict

The headline claim stands: **on subtype, federation closes the deployment gap** (local models on
other sites' patients 0.837 → FedAvg 0.976, against centralized 0.970). The site-signature
mechanism is properly controlled. Two things need guarding, one of which the report already
handles correctly in prose but leaves exposed in the evidence bundle.

## What the report gets right, and is worth stating

The report does **not** over-generalise federation to all three endpoints. §3 says plainly that
no fraction-of-gap-recovered is claimed for EGFR and OS, and gives the weaker honest statement
instead — FedAvg is never worse than the best alternative on any endpoint (EGFR 0.731 vs
centralized 0.638; OS C-index 0.600 vs 0.581). That is the defensible claim and it is the one
made.

## Finding 1 — the `gap_recovered_*` ratios are unstable and should be marked non-quotable

`evidence/results.json` exposes two derived ratios per endpoint. For subtype they behave:

| Endpoint | `gap_recovered_fedavg` | `gap_recovered_fedavg_homeunion` |
|---|---|---|
| subtype | 1.046 | 0.796 |
| EGFR | **8.706** | **−0.862** |
| OS | 1.689 | 0.110 |

Both are of the form `(fedavg − baseline) / (centralized − baseline)`. When the denominator is
near zero — which is exactly the EGFR and OS situation, where centralized barely beats local —
the ratio explodes or flips sign. EGFR reads 8.71 under one definition and −0.86 under the
other. Neither is interpretable.

The report never quotes these, and says so. But a reader who pulls `results.json` directly — the
whole point of shipping the bundle — can lift "8.7" as though federation recovered 870% of the
gap, or "−0.86" as though federation actively harmed EGFR. **Recommendation:** the two ratios
should carry an explicit `interpretable: false` flag for EGFR and OS in any future bundle, and
the two competing definitions should not be shipped side by side without the caveat that the
home-union variant is the conservative one. Nothing in the report changes.

## Finding 2 — EGFR is underpowered, and the bundle says so

`egfr_power` records 74 mutants against 852 wild-type, `trigger_lt80 = True`, and a minimum
detectable AUROC of **0.584** at the permutation null's 95th percentile. So the EGFR arm could
only ever have detected a fairly strong effect. This is consistent with R13's modest EGFR
numbers (0.712 zero-shot → 0.764 at k=100) and is the correct reading of both.

Any statement about EGFR from R12 or R13 must carry the power bound. It is present in the
evidence; it should be repeated wherever the number is.

## Finding 3 — the site-signature probe is sound (no correction)

This is the mechanism claim on which the whole site-shift story rests, so it was checked closely:

- balanced accuracy **0.706** identifying a slide's silo, against a **0.2** chance rate
- permutation null over 100 shuffles centred at **0.203** — i.e. exactly chance
- null 95th percentile 0.222, observed 0.706, permutation **p = 0.0099**

The null being centred at chance is what makes this credible rather than an artefact of class
imbalance or feature scaling. No correction needed. This is the strongest single result in the
report and it is correctly controlled.

## Finding 4 — a container crash and rerun, properly disclosed

`provenance.cost_estimate_usd.note` records that the grid ran 02:36–05:56 CEST and the
permutation null was rerun 08:05–09:10 after a 7200-second container-timeout crash, cross-
referenced to `gates.json` and the run summary. Cost $18 against a $40 budget.

Disclosed, cross-referenced, and the rerun is the permutation null rather than the primary
result — so it cannot have been a retry-until-significant. No issue. Worth recording only
because an undisclosed rerun would have been serious.

## Finding 5 — what "hospital" means here, and the limit it imposes

The five silos are **partitions of TCGA anchored on real tissue-source sites**, not independent
acquisitions. They therefore share TCGA's central processing, and the site signature they carry
is whatever survives that. A reviewer will lead with this, and correctly.

Two things bound the concern. The signature is real and measurable regardless of its origin
(Finding 3), and R13 provides a genuinely independent institution set (CPTAC, 432 patients)
where the zero-shot result holds. But the word "hospital" is doing work the design does not
fully support, and the paper drafted from this must say "site-anchored partitions of a single
collection" wherever it currently says silo or hospital.

**Follow-up now in hand:** the 2026-08-14 public sweep identified the only two genuinely
independent public NSCLC WSI cohorts — CDDP_EAGLE-1 (49 cases, Italian population-based) and
CGCI-HTMCP-LC (39 cases, HIV+ cohort). Manifests are built. They are far too small to be a
powered external benchmark and will not be presented as one, but they can test whether the site
signature persists across genuinely external acquisition, which is the specific thing this
audit flags as unproven.

## Limitations to state to the colleague

1. **The federation is simulated.** Weights never crossed a real institutional boundary; no IT,
   governance or scanner heterogeneity was involved. The 301 MiB-per-run communication figure is
   a real measurement of the protocol, not of a deployment.
2. **Only subtype supports the gap-recovery claim.** EGFR and OS are weak signals for every
   method tested, and the report is right not to claim recovery there.
3. **EGFR is power-limited** at 74 mutants (Finding 2).
4. **Five silos, one collection** (Finding 5).
5. **Attention-MIL over Phikon-v2 features throughout** — conclusions may not transfer to other
   encoders, and the site signature in particular could be encoder-specific.
