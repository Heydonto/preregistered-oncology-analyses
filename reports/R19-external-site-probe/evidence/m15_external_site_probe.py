#!/usr/bin/env python3
"""R19 - the external site probe, re-run because Paper 1 quotes it without evidence.

WHY THIS EXISTS. Paper 1 §"Site identity is recoverable" reports three probe rows:

    Five site-anchored partitions within TCGA      0.7062   null 0.2029   p 0.0099
    EAGLE vs HTMCP (two external cohorts)          1.0000   null 0.5052   p 0.0050
    TCGA vs EAGLE vs HTMCP (size-balanced)         0.9957   null 0.3389   p 0.0050

Only the first has a persisted evidence file (R12 `results.json`,
`silo_signature_probe.balanced_accuracy = 0.706242...`). The two external rows -- including the
perfect-separation number that the manuscript leans on hardest -- were computed in a session and
transcribed into the manuscript without an artefact being written. They are unverifiable as they
stand, and "23/23 claims verified" was asserted over them.

AND THE INPUT IS NOT WHAT ITS NAME SAYS. The probe consumed
`quantara-staging/ext_feat/external_he_mean.npz`. Despite `he` in the filename, 3 of its 29 HTMCP
slides are not H&E: `HTMCP-02-06-01032-CHR`, `HTMCP-02-15-01090-CHR` (chromogenic ISH) and
`HTMCP-02-07-01065-06-P16` (p16 immunohistochemistry). A stain difference is the most trivial
shortcut available to a probe asked to tell two cohorts apart, and with balanced accuracy at
exactly 1.000 three mislabelled slides are not obviously harmless. R12's audit flagged this exact
risk -- HTMCP-LC is majority IHC and only its H&E subset is a valid unit -- and the filter was
then not enforced in the file that asserts it was.

So this run does two things: it produces evidence, and it separates the stain confound from the
site signal by running each probe twice, once over the published slide set and once over the
H&E-only subset.

  ARM_ALL      49 EAGLE + 29 HTMCP  = 78 slides   (the published set, 3 non-H&E)
  ARM_HE       49 EAGLE + 26 HTMCP  = 75 slides   (H&E only)
  ARM3_ALL     TCGA vs EAGLE vs HTMCP, size-balanced, published set
  ARM3_HE      the same, H&E only

The original probe's hyperparameters were never recorded, which is part of the finding. This run
declares its own -- standardise, logistic regression C=1.0, 5-fold stratified out-of-fold
balanced accuracy, 200 label permutations, p = (1 + #{null >= obs}) / (B + 1) -- and gates the
reproduction rather than assuming it.

Run:  python3 m15_external_site_probe.py
"""
import collections
import json
import os

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from audit_kit import Run

EXT = "/Users/rezanehzati/quantara-staging/ext_feat/external_he_mean.npz"
TCGA = "/Users/rezanehzati/Projects/quantara/nsclc-rwpr-study/armC/data/runs/mean_embeddings.npz"

SEED = 20260817
B_PERM = 200
N_SUBSAMPLE = 20          # repeats for the size-balanced 3-way, to damp sampling variance
C_REG = 1.0
N_SPLITS = 5

PUBLISHED = {"eagle_vs_htmcp": 1.0000, "eagle_vs_htmcp_null": 0.5052,
             "three_way": 0.9957, "three_way_null": 0.3389,
             "within_tcga": 0.7062}
REPRO_TOL = 0.02
# pre-declared decision thresholds for the H&E-only arm
HE_CONFIRM = 0.95
HE_WEAK = 0.75


def probe(X, y, seed):
    """Out-of-fold balanced accuracy of a standardised logistic probe."""
    y = np.asarray(y)
    n_min = min(collections.Counter(y).values())
    k = min(N_SPLITS, n_min)
    if k < 2:
        return float("nan")
    oof = np.empty(len(y), dtype=y.dtype)
    for tr, te in StratifiedKFold(n_splits=k, shuffle=True, random_state=seed).split(X, y):
        m = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=5000, C=C_REG)).fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    return float(balanced_accuracy_score(y, oof))


def with_null(X, y, seed, label):
    obs = probe(X, y, seed)
    rng = np.random.default_rng(seed)
    null = np.array([probe(X, rng.permutation(y), seed) for _ in range(B_PERM)])
    p = float((1 + int((null >= obs).sum())) / (B_PERM + 1))
    return {"label": label, "balanced_accuracy": obs, "n": int(len(y)),
            "classes": {str(k): int(v) for k, v in sorted(collections.Counter(y).items())},
            "chance": round(1.0 / len(set(y)), 4),
            "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)),
            "null_p95": float(np.quantile(null, 0.95)), "permutation_p": p,
            "n_permutations": B_PERM}


def balanced_three_way(Xt, Xe, Xh, seed):
    """Equal n per cohort, averaged over N_SUBSAMPLE draws of the two larger cohorts."""
    n = min(len(Xt), len(Xe), len(Xh))
    rng = np.random.default_rng(seed)
    accs, nulls = [], []
    for r in range(N_SUBSAMPLE):
        X = np.vstack([Xt[rng.choice(len(Xt), n, replace=False)],
                       Xe[rng.choice(len(Xe), n, replace=False)],
                       Xh[rng.choice(len(Xh), n, replace=False)]])
        y = np.array(["tcga"] * n + ["eagle"] * n + ["htmcp"] * n, object)
        accs.append(probe(X, y, seed + r))
        nulls.append(probe(X, np.random.default_rng(seed + 1000 + r).permutation(y), seed + r))
    accs, nulls = np.array(accs), np.array(nulls)
    p = float((1 + int((nulls >= accs.mean()).sum())) / (len(nulls) + 1))
    return {"balanced_accuracy": float(accs.mean()), "sd_over_subsamples": float(accs.std(ddof=1)),
            "per_cohort_n": int(n), "chance": round(1 / 3, 4),
            "null_mean": float(nulls.mean()), "null_sd": float(nulls.std(ddof=1)),
            "permutation_p": p, "n_subsamples": N_SUBSAMPLE}


def main():
    run = Run("R19extsiteprobe")
    cfg_hash = run.start(
        {"question": "does the external site signature survive restriction to H&E, and do "
                     "Paper 1's two unevidenced probe rows reproduce?",
         "why": "Paper 1 quotes 1.0000 and 0.9957 with no persisted artefact, over an input file "
                "named external_he_mean.npz that contains 3 non-H&E slides (2 CHR, 1 P16).",
         "arms": {"ARM_ALL": "EAGLE vs HTMCP, published 78-slide set",
                  "ARM_HE": "EAGLE vs HTMCP, H&E only (75)",
                  "ARM3_ALL": "3-way size-balanced, published set",
                  "ARM3_HE": "3-way size-balanced, H&E only"},
         "probe": {"model": "StandardScaler + LogisticRegression", "C": C_REG,
                   "cv": f"StratifiedKFold({N_SPLITS}, shuffled)",
                   "metric": "out-of-fold balanced accuracy",
                   "note": "the original probe's hyperparameters were never recorded; these are "
                           "declared here and the reproduction is gated, not assumed"},
         "permutations": B_PERM, "p_formula": "(1 + #{null >= obs}) / (B + 1)",
         "seed": SEED, "published": PUBLISHED, "repro_tol": REPRO_TOL,
         "verdict_rules": [
             f"ARM_HE >= {HE_CONFIRM} -> EXTERNAL_SITE_SIGNATURE_CONFIRMED_HE_ONLY",
             f"ARM_HE >= {HE_WEAK} -> PRESENT_BUT_WEAKER_THAN_PUBLISHED",
             f"ARM_HE < {HE_WEAK} -> PUBLISHED_NUMBER_WAS_STAIN_DRIVEN"],
         "status": "peer-review audit copy"},
        [EXT, TCGA])

    d = np.load(EXT, allow_pickle=True)
    sid = np.array([str(x) for x in d["slide_ids"]], object)
    coh = np.array([str(x) for x in d["cohort"]], object)
    X = d["mean_emb"].astype(np.float64)
    run.gate("G0_input", "78 slides x 1024 finite features",
             {"shape": list(X.shape), "finite": bool(np.isfinite(X).all()),
              "cohorts": {k: int(v) for k, v in sorted(collections.Counter(coh).items())}},
             X.shape == (78, 1024) and bool(np.isfinite(X).all()))

    # H&E determination from the slide-id stain suffix, the only stain field available
    is_he = np.array([(not s.startswith("HTMCP")) or s.rsplit("-", 1)[-1] == "HE" for s in sid])
    nonhe = sorted(sid[~is_he].tolist())
    run.gate("G1_stain_contamination_confirmed",
             "the published set contains non-H&E HTMCP slides (this is the defect, so the gate "
             "asserts it is present and names them)",
             {"n_non_he": len(nonhe), "slides": nonhe}, len(nonhe) == 3,
             "if this ever returns 0 the input file was silently changed and the whole premise "
             "of this run needs re-checking")

    dt = np.load(TCGA, allow_pickle=True)
    Xt = dt["mean_emb"].astype(np.float64)
    run.gate("G2_tcga_input", "TCGA slide-level mean embeddings, 1024-d",
             {"shape": list(Xt.shape)}, Xt.shape[1] == 1024 and Xt.shape[0] > 1000)

    # ---------------- two-cohort probe ----------------
    res = {}
    res["ARM_ALL"] = with_null(X, coh, SEED, "EAGLE vs HTMCP, published 78")
    res["ARM_HE"] = with_null(X[is_he], coh[is_he], SEED, "EAGLE vs HTMCP, H&E only")
    for k in ("ARM_ALL", "ARM_HE"):
        run.log("two_cohort", arm=k, bal_acc=round(res[k]["balanced_accuracy"], 4),
                n=res[k]["n"], null_mean=round(res[k]["null_mean"], 4),
                p=res[k]["permutation_p"])

    run.gate("G3_null_centred_two_cohort",
             "both two-cohort nulls sit at chance 0.5 (within 0.05)",
             {k: round(res[k]["null_mean"], 4) for k in ("ARM_ALL", "ARM_HE")},
             all(abs(res[k]["null_mean"] - 0.5) < 0.05 for k in ("ARM_ALL", "ARM_HE")))

    run.gate("G4_reproduces_paper1_two_cohort",
             f"published 1.0000 recovered on the published slide set, tol {REPRO_TOL}",
             round(res["ARM_ALL"]["balanced_accuracy"], 4),
             abs(res["ARM_ALL"]["balanced_accuracy"] - PUBLISHED["eagle_vs_htmcp"]) <= REPRO_TOL,
             "the original hyperparameters are unrecorded, so a miss here means the manuscript "
             "number is not reproducible from the surviving inputs -- which would itself be the "
             "result, and is why this gate is checked before the H&E arm is interpreted")

    # ---------------- three-cohort probe ----------------
    Xe, Xh = X[coh == "eagle"], X[coh == "htmcp"]
    Xh_he = X[(coh == "htmcp") & is_he]
    res["ARM3_ALL"] = balanced_three_way(Xt, Xe, Xh, SEED)
    res["ARM3_HE"] = balanced_three_way(Xt, Xe, Xh_he, SEED)
    for k in ("ARM3_ALL", "ARM3_HE"):
        run.log("three_way", arm=k, bal_acc=round(res[k]["balanced_accuracy"], 4),
                per_cohort_n=res[k]["per_cohort_n"], null_mean=round(res[k]["null_mean"], 4),
                p=res[k]["permutation_p"])
    run.gate("G5_null_centred_three_way", "both 3-way nulls sit at chance 0.3333 (within 0.05)",
             {k: round(res[k]["null_mean"], 4) for k in ("ARM3_ALL", "ARM3_HE")},
             all(abs(res[k]["null_mean"] - 1 / 3) < 0.05 for k in ("ARM3_ALL", "ARM3_HE")))

    # ---------------- mechanical verdict ----------------
    he = res["ARM_HE"]["balanced_accuracy"]
    if he >= HE_CONFIRM:
        verdict = "EXTERNAL_SITE_SIGNATURE_CONFIRMED_HE_ONLY"
    elif he >= HE_WEAK:
        verdict = "PRESENT_BUT_WEAKER_THAN_PUBLISHED"
    else:
        verdict = "PUBLISHED_NUMBER_WAS_STAIN_DRIVEN"
    delta = res["ARM_ALL"]["balanced_accuracy"] - he
    run.log("VERDICT", primary=verdict, he_only=round(he, 4),
            published_set=round(res["ARM_ALL"]["balanced_accuracy"], 4),
            stain_contribution=round(delta, 4))

    out = {
        "status": "peer-review audit copy",
        "question": "does the external site signature survive H&E restriction, and does Paper 1's "
                    "unevidenced 1.0000 reproduce?",
        "VERDICT": {"primary": verdict,
                    "he_only_balanced_accuracy": he,
                    "published_set_balanced_accuracy": res["ARM_ALL"]["balanced_accuracy"],
                    "difference_attributable_to_non_he_slides": delta,
                    "rules_pre_declared_in": "config.yaml, hash " + cfg_hash},
        "arms": res,
        "stain_contamination": {"n_non_he_in_published_set": len(nonhe), "slides": nonhe,
                                "source_file": EXT,
                                "note": "filename asserts H&E; contents do not"},
        "paper1_rows": {
            "row1_within_tcga": {"published": PUBLISHED["within_tcga"],
                                 "evidence": "R12 evidence/results.json "
                                             "silo_signature_probe.balanced_accuracy",
                                 "status": "already evidenced, not recomputed here"},
            "row2_eagle_vs_htmcp": {"published": PUBLISHED["eagle_vs_htmcp"],
                                    "recomputed_published_set":
                                        res["ARM_ALL"]["balanced_accuracy"],
                                    "recomputed_he_only": he},
            "row3_three_way": {"published": PUBLISHED["three_way"],
                               "recomputed_published_set":
                                   res["ARM3_ALL"]["balanced_accuracy"],
                               "recomputed_he_only": res["ARM3_HE"]["balanced_accuracy"]}},
        "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir},
    }
    np.savetxt(os.path.join(run.dir, "probe_set.csv"),
               np.column_stack([sid, coh, is_he.astype(int)]), fmt="%s", delimiter=",",
               header="slide_id,cohort,is_he", comments="")
    run.write("results.json", out)
    run.finalize()
    print("\nVERDICT:", verdict)
    print(f"  published set {res['ARM_ALL']['balanced_accuracy']:.4f} -> "
          f"H&E only {he:.4f}  (delta {delta:+.4f})")
    print("RUN DIR:", run.dir)


if __name__ == "__main__":
    main()
