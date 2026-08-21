# Corrections

Every correction we made to our own work, with the original wording preserved. Ordered by date.
Each is also recorded in the audit of the report it affects.

This file exists because a reviewer's first job is to find what a group has hidden. Handing over
the list is cheaper for everyone, and the pattern in it is more informative than any single entry.

---

## 1. A censoring error (R11, 2026-08)

**Original:** used a raw median survival of 11.8 months for a TCGA cohort that is ~50% censored,
making pan-glioma appear to survive worse than GBM.

**Corrected to:** the Kaplan–Meier median, 20.7 months. Gate G13 added to assert that a reported
median is a KM estimate whenever censoring exceeds a threshold.

---

## 2. "The negative was an artefact of our own supervision choice" (R15 audit → R17)

**Original, in R15's audit:** the genome-wide null was attributed to having supervised the
representation on the wrong target, and a general rule was drawn from it.

**What R17 found:** re-running the identical scan under methylation supervision gives 83,369 CpGs
above null against 90,137 for the subtype-supervised one — not better, marginally worse. The
reversal R15 observed is real but specific to its six *aggregate* targets. R17 supplies the
mechanism: a patient's observed global mean methylation alone clears 74.7% of CpGs, so a
representation can learn one axis and carry nothing at CpG resolution.

**Note:** R15's own limitations section stated "per-CpG prediction under methylation supervision
was not attempted." The gap was disclosed and the summary over-claimed anyway. **Stating a
limitation does not bound a claim.**

---

## 3. Two manuscript numbers with no evidence file (R19)

**Original:** Paper 1 quoted balanced accuracies of 1.0000 and 0.9957 for external site probes.
Neither had a persisted artefact — both were computed in a session and typed into the manuscript.
The verification pass had marked them verified *because they matched what we believed*.

**Corrected to:** both re-run and evidenced. 1.0000 reproduces. 0.9957 does **not** — it
recomputes to 0.9816, and because the original probe's hyperparameters were never recorded the
0.014 gap is permanently unattributable. Paper 1 now quotes the reproducible value with its method.

**Also found:** the input file named `external_he_mean.npz` contained three non-H&E slides. The
H&E-only result is unchanged at 1.0000, so the concern was legitimate and the science survived.

---

## 4. "Real, total, and unchanged" (R19 audit → R20/R22)

**Original:** the external site signature was described as *total*.

**Corrected to:** total under Phikon-v2 only. On the same 76 slides a natural-image encoder reaches
0.9176 ± 0.0247 over 25 splits. The word "perfectly" now carries an encoder in Paper 1.

---

## 5. "A property of Phikon-v2 features" (R18 → R21)

**Original:** from two encoders, methylation site-leakage was attributed to one vendor's model.

**Corrected to:** it is a property of **histology pretraining**. All four histology encoders inflate
all six targets (+0.0173 to +0.0849); the one natural-image encoder does not (+0.0004, 3/6).

**Note:** R18's own limitations section reads "two encoders is not a survey," and the inference was
drawn anyway. Same shape as correction 2.

---

## 6. A single cross-validation split quoted as a fixed quantity (R20 → R22)

**Original:** an external separation reported as 0.8768.

**Corrected to:** 0.9176 ± 0.0247 over 25 splits. The 0.8768 was the pre-declared seed — so not
cherry-picked — but sat 1.6 standard deviations below its own mean, and the headline contrast was
overstated. A cross-validated metric on 76 samples is a random variable.

**Found by:** a reproduction gate written to catch a coding difference, which caught a statistical
one instead and halted the run.

---

## 7. A pre-registration that made its own hypotheses untestable (R23)

**Original:** the plan specified 200 permutations. With a p-floor of 1/201, Holm correction across
11 drugs cannot go below 0.0547 and across 32 cannot go below 0.159. **Significance was
arithmetically unreachable before any data was touched**, and the first run reported zero
predictable drugs while one sat at ρ = 0.73.

**Corrected to:** 2,000 permutations via a kernel reformulation. Gate G5/G6 now asserts the
permutation floor permits significance at the declared alpha, before any test runs.

---

## 8. A pre-registered verdict that measured a confound (R23)

**Original:** the primary test returned ρ = +0.651, p = 0.0005 — THESIS_CONTRADICTED under the
hashed rule, every gate passing.

**What a post-hoc control showed:** both indices load on a general drug-sensitivity axis.
Partialling out an index from 227 drugs in neither set takes +0.651 to **−0.160**. The sign
reverses.

**Resolution:** neither number is reported as the answer. The post-hoc figure is not a
pre-registered test and cannot replace one; what the pair licenses is that the design was
confounded and answers nothing. **Pre-registration constrains analytic flexibility; it does not
confer construct validity.**

---

## 9. A figure claim about the wrong figure (internal, 2026-08-20)

**Original:** stated that one figure of a third-party script was "partly" covered by our
drug-response analysis.

**Corrected to:** not covered at all. That figure characterises a compound for which no data exists
in any holding — it requires laboratory work, not a dataset.

---

## The pattern

Corrections 2, 4 and 5 are the same error three times: **a real measurement, a defensible-sounding
generalisation one step beyond it, and the generalisation wrong.** In each case the limitation that
would have caught it was already written in the same document.

Corrections 6, 7 and 8 are a different family: the plan or the statistic was mis-specified rather
than over-read.

The rule we now apply: *a limitations section is not a hedge on the summary, it specifies which
sentences are not yet licensed. If a limitation names a comparison you have not run, no sentence may
assume its result — delete the sentence rather than add a caveat.*
