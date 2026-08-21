# R23 audit — a pre-registered verdict killed by a control the plan forgot

Audited 20 August 2026. Plan `PREREG_PLAN.yaml`, SHA-256 `7e0322ba804aef2f…`, hashed before any
IC50 was joined to any methylation value. All eight gates passed. Held for IP.

## Verdict

Three real findings and one honest non-answer.

Methylation predicts ln IC50 for 10 of 11 standard-of-care agents and 14 of 32 chromatin drugs.
The pre-registered primary test returned **THESIS_CONTRADICTED** (ρ = +0.651, p = 0.0005). A
post-hoc control then took that to **−0.160** — sign reversed, below threshold. **The question is
unresolved, and saying so is the result.**

## Finding 1 — the plan made its own hypotheses untestable

The first execution reported **zero** predictable drugs while Vorinostat sat at ρ = 0.73.

Not a data property. The plan specified 200 permutations, so the p-floor is 1/201 = 0.00498. Holm
across 11 drugs cannot go below 0.0547; across 32 it cannot go below 0.159. **Significance was
arithmetically unreachable before a single number was computed.**

I wrote that plan. It passed my own review, it was hashed, and the flaw is visible from the two
numbers it contains — permutation count and family size — without any data at all.

**Generalisable rule:** a pre-registration must be checked for *internal feasibility* before it is
hashed. At minimum: can the declared test, with the declared correction and the declared
permutation count, reject at the declared alpha? That is arithmetic, it takes one line, and it is
now gate G6, which fails loudly if the floor is ever too coarse for the family.

Raising to 2,000 permutations is a departure from the hashed plan. It is logged, and it can only
make rejection *easier* — a change in the direction that costs me something, not one that helps.

## Finding 2 — the primary test measured a confound

Both indices are mean z-scored ln IC50. Cell lines differ enormously in general drug sensitivity:
some are resistant to nearly everything. Two such indices will correlate whether or not chromatin
drugs are special.

Partialling out an index built from **227 GDSC2 drugs in neither set**:

| | ρ | p |
|---|---|---|
| pre-registered, raw | **+0.651** | 0.0005 |
| post-hoc, general sensitivity removed | **−0.160** | 0.033 |

The sign reverses. +0.651 was almost entirely the general axis.

**What I am not entitled to do is swap one for the other.** The pre-registered verdict is
THESIS_CONTRADICTED and that is what the hashed rule returns. The post-hoc figure is not a
pre-registered test and cannot replace it. What both together license is: *the design was
confounded and answers nothing.* Reporting only +0.651 would assert the thesis is wrong; reporting
only −0.160 would assert it is weakly right; both would be misrepresentations.

**And the control was mine, added after seeing +0.651 looked too clean.** Had I not been
suspicious, the report would have said the manuscript's central premise is contradicted at
p = 0.0005 in 187 cell lines. That is a publishable-looking, wrong conclusion, produced by a
correctly executed pre-registration.

**The lesson is uncomfortable and it belongs in Paper 2.** Pre-registration constrains analytic
flexibility. It does not confer construct validity. A hashed plan that measures the wrong thing
measures the wrong thing with excellent discipline.

## Finding 3 — the per-drug table carries the same confound, uncontrolled

Ten of eleven standard-of-care agents predictable, ρ up to 0.607. Fourteen of 32 chromatin drugs,
ρ up to 0.731. Real, Holm-corrected over 2,000 permutations, nulls centred within 0.021.

But a model predicting "this line is broadly resistant" will score well on almost any drug, and
nothing here separates that from drug-specific prediction. The plan did not require it. The report
states the numbers as "methylation carries information about this line's IC50", not as "predicts
response to this agent" — and that is the strongest claim the design supports.

Consistent with this: the drugs with the highest ρ (Vorinostat, Navitoclax) are broadly cytotoxic,
while the lowest (Gemcitabine, 0.177, not significant) is the most schedule- and
metabolism-dependent.

## Finding 4 — our copy of the source archive is corrupt and says nothing about it

`GSE68379_RAW.tar` is byte-complete against both GCS and GEO's declared size
(9,020,579,840), ends with a valid tar terminator, and raises no exception — while enumerating
**1,302 of 2,056 members**: 649 of 1,028 cell lines, 127 of 198 lung. Two extraction modes gave the
identical truncation; a dedicated diagnostic confirmed the archive, not the reader.

Caught **only** because the job asserted an expected array count and refused to proceed. Without
that assertion this report would rest on 127 lines and would never have mentioned it. The 396 lung
IDATs were fetched per sample from GEO instead: 396 of 396, zero failures.

This is the fourth defect in this programme of the form *an input that does not contain what its
label says*, after the three in `external_he_mean.npz`. **The data inventory describes this archive
as 1,028 cell lines — true of GEO's copy, false of ours** — and needs correcting. Every other bulk
archive in the inventory should be assumed suspect until member-counted.

## Three bugs of mine, for completeness

Decompressing by file extension rather than magic number, which killed the first run on one corrupt
member. Blaming streaming tar mode when the archive was at fault. Pinning pandas 2.x against a
methylprep that still calls `DataFrame.append`.

## Limitations to state to the colleague

1. **Cell lines are not tumours.** Methylation drifts in culture; ln IC50 is a lab phenotype.
2. **The primary question is unresolved.** A corrected pre-registration should build both indices
   on general-sensitivity residuals and declare that in advance.
3. **Per-drug associations are not drug-specific** and are not claimed to be.
4. **Mixed histology** — 58 small cell and 21 mesothelioma among 187; no stratified analysis was
   pre-declared.
5. **Sex-chromosome probes were not excluded**; the 450k manifest is not staged locally, so
   exclusion E5 went unapplied.
6. **Hash-before-label only.** Plan author and executor were the same person — the weaker form of
   the protocol, not R16's structural isolation.
