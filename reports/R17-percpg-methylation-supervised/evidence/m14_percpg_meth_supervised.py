#!/usr/bin/env python3
"""R17 - genome-wide per-CpG scan under METHYLATION supervision.

R15/m13 ran the genome-wide scan on a SUBTYPE-supervised embedding and found nothing
(median rho 0.006; 110,212/399,579 CpGs above null, against 220,925 for the subtype label
alone). Audit Finding 2 concluded that negative was an artefact of the supervision choice,
not a fact about H&E: the representation had been optimised to separate two histological
classes and nothing in training asked it to preserve methylation structure.

That conclusion was drawn from the six aggregate targets. It was never tested genome-wide,
because the per-CpG scan was only ever run on the wrong representation. This run closes that
gap by re-running the identical scan on the methylation-supervised embedding.

Comparability is the point, so everything that can be held fixed is held fixed: the same 760
patients in the same order, the same site-grouped folds (taken from `fold_of` as written by
the MIL, verified equal to GroupKFold(5) by site), the same ridge (alpha 100, chunked
closed-form), the same probe filter, the same median imputation, the same R10 genomic-context
scheme. Only the embedding changes.

THE CIRCULARITY, STATED UP FRONT. The methylation-supervised embedding was trained on six
aggregate targets, one of which ('global') is the mean over every probe in this scan. So this
is NOT a clean test of "morphology predicts the methylome" -- the representation was built to
encode methylation. Two things bound the problem:

  * It is not patient-level leakage. Z is out-of-fold: each patient's embedding comes from a
    model that never saw that patient, and the ridge folds are the SAME folds, so they nest.
  * The live worry is a scalar shortcut -- that the embedding only learned each patient's
    overall methylation level, and individual CpGs score above null merely by correlating
    with that level.

Arm A_globalmean is the control for exactly that: it predicts every CpG from the patient's
OBSERVED global mean methylation, i.e. ground truth for the dominant axis. If the embedding
does not beat that baseline, the honest reading is 'global axis only' and the verdict says so.

Arms (identical folds, identical ridge, only the predictor differs):
  A_meth        512-d methylation-supervised embedding      <- the question
  A_sub         512-d subtype-supervised embedding          <- reproduces m13
  A_sublabel    1-d LUAD/LUSC label                         <- reproduces m13
  A_globalmean  1-d observed global mean methylation        <- circularity control
  A_null        A_meth with the WSI-methylation pairing permuted

m13's null used two independent permutations (one for the observed matrix, one for the ridge
target), so it scored Y[p1] against a model fit on Y[p2]. That is still a valid null -- both
are decoupled from Z -- but it is not the single-permutation null it reads as. This run uses
one permutation and ALSO recomputes m13's two-permutation variant, so the effect on the
threshold is visible rather than assumed.

Run:  python3 m14_percpg_meth_supervised.py
"""
import csv
import gzip
import json
import os

import numpy as np
from scipy.stats import rankdata, spearmanr
from sklearn.model_selection import GroupKFold

from audit_kit import Run, sha256

R08 = "/Users/rezanehzati/quantara-staging/r08/"
ARM = "/Users/rezanehzati/Projects/quantara/nsclc-rwpr-study/armD-meth/"
R15 = "/Users/rezanehzati/Projects/quantara/results/R15-wsi-methylation/evidence/"
MANIFEST = "/tmp/manifest450k.csv.gz"

SEED = 20260817
CHUNK = 40000
RIDGE_ALPHA = 100.0

# from R15 / m13, run 20260814T211255Z-R15wsimeth-3f8bece
EXPECT = {
    "matched": 760,
    "probes": 399579,
    "m13_n_above_null_subtype_emb": 110212,
    "m13_n_above_null_subtype_label": 220925,
    "m13_median_rho_subtype_emb": 0.006,
    "r10_context": {"Island": 0.54, "OpenSea": 1.44},
}
# pre-declared decision thresholds
MARGIN = 1.25          # A_meth must exceed A_sub by this factor to call supervision the cause
REPRO_MIN_CORR = 0.999  # elementwise agreement required against m13's saved rho vectors


def per_cpg_rho(Y, P):
    """Spearman per column, vectorised. Exactly equivalent to looping spearmanr:
    average-rank transform each column, then Pearson. Gate G3 checks that against m13."""
    Ry = rankdata(Y, axis=0).astype(np.float32)
    Rp = rankdata(P, axis=0).astype(np.float32)
    Ry -= Ry.mean(0)
    Rp -= Rp.mean(0)
    num = (Ry * Rp).sum(0)
    den = np.sqrt((Ry * Ry).sum(0) * (Rp * Rp).sum(0))
    del Ry, Rp
    return np.nan_to_num(num / np.maximum(den, 1e-12)).astype(np.float32)


def ridge_cv(X, Yv, splits):
    """Out-of-fold ridge, chunked over CpGs. Identical to m13's."""
    X = np.atleast_2d(np.asarray(X, np.float64))
    if X.shape[0] != Yv.shape[0]:
        X = X.T
    pred = np.full(Yv.shape, np.nan, np.float32)
    for tr, te in splits:
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-8
        A, B_ = (X[tr] - mu) / sd, (X[te] - mu) / sd
        G = A.T @ A
        for s0 in range(0, Yv.shape[1], CHUNK):
            Yb = Yv[tr, s0:s0 + CHUNK].astype(np.float64)
            m = Yb.mean(0)
            W = np.linalg.solve(G + RIDGE_ALPHA * np.eye(G.shape[0]), A.T @ (Yb - m))
            pred[te, s0:s0 + CHUNK] = (B_ @ W + m).astype(np.float32)
    return pred


def main():
    run = Run("R17percpgmeth")
    cfg_hash = run.start(
        {"question": "does a methylation-supervised morphology representation retain "
                     "genome-wide methylation structure, where a subtype-supervised one "
                     "did not (R15 Finding 2)?",
         "design": "five arms, identical site-grouped folds and identical ridge; only the "
                   "predictor differs. Folds taken from the MIL's own fold_of so the ridge "
                   "folds nest inside the representation's folds.",
         "arms": {"A_meth": "512-d methylation-supervised embedding",
                  "A_sub": "512-d subtype-supervised embedding (reproduces m13)",
                  "A_sublabel": "1-d LUAD/LUSC label (reproduces m13)",
                  "A_globalmean": "1-d OBSERVED global mean methylation (circularity control)",
                  "A_null": "A_meth with WSI-methylation pairing permuted, single permutation"},
         "circularity": "the 'global' training target is the mean over every scanned probe, so "
                        "this is not a clean 'morphology predicts methylome' test. Z is "
                        "out-of-fold so there is no patient-level leakage; A_globalmean bounds "
                        "the scalar-shortcut explanation.",
         "seed": SEED, "ridge_alpha": RIDGE_ALPHA, "chunk": CHUNK,
         "margin_factor": MARGIN, "repro_min_corr": REPRO_MIN_CORR,
         "expectations": EXPECT,
         "verdict_rules": [
             "n_meth > MARGIN*n_sub AND n_meth > n_sublabel -> SUPERVISION_EXPLAINS_ARM1_NEGATIVE",
             "abs(n_meth-n_sub) <= 0.25*n_sub -> ARM1_NEGATIVE_ROBUST_TO_SUPERVISION",
             "otherwise -> INTERMEDIATE",
             "secondary, independent of the above: n_meth > n_globalmean -> "
             "CARRIES_CPG_SPECIFIC_STRUCTURE_BEYOND_GLOBAL_AXIS, else GLOBAL_AXIS_ONLY"],
         "status": "peer-review audit copy"},
        [ARM + "mil_meth_output.npz", ARM + "mil_output.npz", ARM + "mil_input.json",
         R08 + "beta_450k.npy", R08 + "samples.tsv", R15 + "per_cpg_rho_grouped.npy",
         R15 + "per_cpg_rho_subtype_only.npy", MANIFEST])

    # ---------------- cohort, embeddings, folds (no outcome read yet) ----------------
    dm = np.load(ARM + "mil_meth_output.npz", allow_pickle=True)
    ds = np.load(ARM + "mil_output.npz", allow_pickle=True)
    pm = [str(x) for x in dm["pids"]]
    ps = [str(x) for x in ds["pids"]]
    run.gate("G0_cohort", EXPECT["matched"], len(pm), len(pm) == EXPECT["matched"])
    run.gate("G0b_same_order", "both arms cover the same patients in the same order",
             pm == ps, pm == ps,
             "the two embeddings are compared elementwise, so order equality is load-bearing")

    Z_meth = dm["Z"].astype(np.float64)
    Z_sub = ds["Z"].astype(np.float64)
    fold = dm["fold_of"].astype(int)
    site = np.array([p[5:7] for p in pm], object)

    derived = np.empty(len(pm), int)
    for k, (_, te) in enumerate(GroupKFold(n_splits=5).split(Z_meth, groups=site)):
        derived[te] = k
    run.gate("G1_folds_nest",
             "MIL fold_of == GroupKFold(5) by site, and identical across both arms",
             {"meth_matches_groupkfold": bool((fold == derived).all()),
              "arms_identical": bool((fold == ds["fold_of"].astype(int)).all())},
             bool((fold == derived).all()) and bool((fold == ds["fold_of"].astype(int)).all()),
             "using the MIL's own folds makes the ridge folds nest inside the "
             "representation's folds; re-deriving them would only be safe if they agree")

    splits = [(np.where(fold != k)[0], np.where(fold == k)[0]) for k in range(5)]
    for tr, te in splits:
        assert not (set(site[tr]) & set(site[te]))
    run.gate("G2_sites_disjoint", "train/test sites disjoint in every fold", "verified", True)

    lab = {p["pid"]: p for p in json.load(open(ARM + "mil_input.json"))["patients"]}
    sub = np.array([lab[p]["subtype"] for p in pm], float)

    # ---------------- methylation matrix, exactly m13's construction ----------------
    smeta = list(csv.DictReader(open(R08 + "samples.tsv"), delimiter="\t"))
    seen = {}
    for r in smeta:
        if r["sample_type"] == "Solid Tissue Normal":
            continue
        if r["case"] not in seen or r["file_id"] < seen[r["case"]]["file_id"]:
            seen[r["case"]] = r
    M = np.load(R08 + "beta_450k.npy", mmap_mode="r")
    probes = open(R08 + "probes.txt").read().split()
    is_t = np.array([r["sample_type"] != "Solid Tissue Normal" for r in smeta])
    nanf = np.isnan(np.asarray(M[:, np.where(is_t)[0][:400]], np.float32)).mean(1)
    good = nanf <= 0.05
    Y = np.asarray(M[:, [int(seen[p]["col"]) for p in pm]], np.float32)[good]
    pnames = [p for p, k in zip(probes, good) if k]
    del M
    Y = np.where(np.isnan(Y), np.nanmedian(Y, axis=1, keepdims=True), Y).T
    run.gate("G3_probes", EXPECT["probes"], int(Y.shape[1]), int(Y.shape[1]) == EXPECT["probes"],
             "same probe filter as m13; a different count would break comparability")

    # ---------------- reproduce m13 before trusting the reimplementation ----------------
    rho_sub = per_cpg_rho(Y, ridge_cv(Z_sub, Y, splits))
    rho_sublabel = per_cpg_rho(Y, ridge_cv(sub.reshape(-1, 1), Y, splits))
    m13_sub = np.load(R15 + "per_cpg_rho_grouped.npy")
    m13_lab = np.load(R15 + "per_cpg_rho_subtype_only.npy")
    c_sub = float(np.corrcoef(rho_sub, m13_sub)[0, 1])
    c_lab = float(np.corrcoef(rho_sublabel, m13_lab)[0, 1])
    run.gate("G4_reproduces_m13",
             f"vectorised Spearman reproduces m13's per-CpG vectors, r > {REPRO_MIN_CORR}",
             {"subtype_embedding_r": round(c_sub, 6), "subtype_label_r": round(c_lab, 6)},
             c_sub > REPRO_MIN_CORR and c_lab > REPRO_MIN_CORR,
             "validates the rank-transform shortcut and the fold_of substitution against the "
             "loop-over-spearmanr implementation that produced the published numbers")

    # spot-check the vectorised Spearman against scipy on real columns
    idx = np.random.default_rng(SEED).choice(Y.shape[1], 200, replace=False)
    P_sub = ridge_cv(Z_sub, Y[:, idx], splits)
    ref = np.array([spearmanr(Y[:, j], P_sub[:, i]).statistic for i, j in enumerate(idx)])
    mad = float(np.abs(ref - per_cpg_rho(Y[:, idx], P_sub)).max())
    run.gate("G5_spearman_exact", "max abs deviation from scipy spearmanr < 1e-5",
             round(mad, 9), mad < 1e-5)

    # ---------------- the arms ----------------
    rho_meth = per_cpg_rho(Y, ridge_cv(Z_meth, Y, splits))
    bulk = Y.mean(1)
    rho_global = per_cpg_rho(Y, ridge_cv(bulk.reshape(-1, 1), Y, splits))

    rng = np.random.default_rng(SEED)
    perm = rng.permutation(len(pm))
    rho_null = per_cpg_rho(Y[perm], ridge_cv(Z_meth, Y[perm], splits))
    # m13's two-permutation variant, recomputed so the threshold difference is visible
    p1, p2 = rng.permutation(len(pm)), rng.permutation(len(pm))
    rho_null_m13style = per_cpg_rho(Y[p1], ridge_cv(Z_meth, Y[p2], splits))

    for nm, v in (("meth", rho_meth), ("sub", rho_sub), ("sublabel", rho_sublabel),
                  ("globalmean", rho_global), ("null", rho_null)):
        run.gate(f"G6_finite_{nm}", "no NaN/inf in the rho vector",
                 int((~np.isfinite(v)).sum()), bool(np.isfinite(v).all()))

    thr = float(np.quantile(rho_null, 0.95))
    thr_m13 = float(np.quantile(rho_null_m13style, 0.95))
    run.gate("G7_null_centred", "single-permutation null median near zero",
             round(float(np.median(rho_null)), 4), abs(float(np.median(rho_null))) < 0.05)

    n = {k: int((v > thr).sum()) for k, v in
         (("meth", rho_meth), ("sub", rho_sub), ("sublabel", rho_sublabel),
          ("globalmean", rho_global))}
    run.log("null_thresholds", single_perm_p95=round(thr, 5),
            m13_two_perm_p95=round(thr_m13, 5),
            single_perm_median=round(float(np.median(rho_null)), 5),
            m13_style_median=round(float(np.median(rho_null_m13style)), 5))
    run.log("HEADLINE_per_cpg", threshold=round(thr, 5),
            **{f"n_above_{k}": v for k, v in n.items()},
            **{f"frac_{k}": round(v / len(rho_meth), 4) for k, v in n.items()},
            median_rho_meth=round(float(np.median(rho_meth)), 4),
            median_rho_sub=round(float(np.median(rho_sub)), 4),
            max_rho_meth=round(float(rho_meth.max()), 4))

    # m13's own counts used m13's own threshold; reproduce them on m13's saved null
    m13_thr = float(np.quantile(np.load(R15 + "per_cpg_rho_null.npy"), 0.95))
    run.gate("G8_reproduces_m13_counts",
             f"m13's published counts recovered on m13's own null threshold "
             f"({EXPECT['m13_n_above_null_subtype_emb']}, "
             f"{EXPECT['m13_n_above_null_subtype_label']}), within 1%",
             {"subtype_emb": int((m13_sub > m13_thr).sum()),
              "subtype_label": int((m13_lab > m13_thr).sum()), "thr": round(m13_thr, 5)},
             abs(int((m13_sub > m13_thr).sum()) - EXPECT["m13_n_above_null_subtype_emb"])
             <= 0.01 * EXPECT["m13_n_above_null_subtype_emb"]
             and abs(int((m13_lab > m13_thr).sum()) - EXPECT["m13_n_above_null_subtype_label"])
             <= 0.01 * EXPECT["m13_n_above_null_subtype_label"])

    # ---------------- genomic context, R10's scheme ----------------
    need = set(pnames)
    ann = {}
    with gzip.open(MANIFEST, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("IlmnID"):
                break
        for line in fh:
            p = line.split(",")
            if len(p) > 25 and p[0].startswith("cg") and p[0] in need:
                ann[p[0]] = p[25].strip() or "OpenSea"
    isl = np.array([ann.get(p, "OpenSea") for p in pnames], object)
    ctx = {}
    for arm, v in (("A_meth", rho_meth), ("A_sub", rho_sub)):
        top = v > thr
        d = {}
        for c in ("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"):
            inc = isl == c
            if not inc.sum():
                continue
            o = float((inc & top).sum() / max(top.sum(), 1))
            e = float(inc.sum() / len(isl))
            d[c] = {"n_probes": int(inc.sum()), "frac_of_predictable": round(o, 4),
                    "frac_of_array": round(e, 4), "enrichment": round(o / e, 3)}
        ctx[arm] = d
        run.log("context", arm=arm,
                **{c: d[c]["enrichment"] for c in d})
    spread = {a: round(max(d[c]["enrichment"] for c in d)
                       - min(d[c]["enrichment"] for c in d), 3) for a, d in ctx.items()}
    run.log("context_spread", note="m13/Arm-1 was flat at 0.97-1.01 (spread 0.04); "
            "R10 from the assay itself was structured (Island 0.54, OpenSea 1.44, spread 0.90)",
            **spread)

    # ---------------- mechanical verdict ----------------
    if n["meth"] > MARGIN * n["sub"] and n["meth"] > n["sublabel"]:
        verdict = "SUPERVISION_EXPLAINS_ARM1_NEGATIVE"
    elif abs(n["meth"] - n["sub"]) <= 0.25 * n["sub"]:
        verdict = "ARM1_NEGATIVE_ROBUST_TO_SUPERVISION"
    else:
        verdict = "INTERMEDIATE"
    secondary = ("CARRIES_CPG_SPECIFIC_STRUCTURE_BEYOND_GLOBAL_AXIS"
                 if n["meth"] > n["globalmean"] else "GLOBAL_AXIS_ONLY")
    run.log("VERDICT", primary=verdict, secondary=secondary)

    results = {
        "status": "peer-review audit copy",
        "question": "is R15's genome-wide null an artefact of subtype supervision?",
        "VERDICT": {"primary": verdict, "secondary": secondary,
                    "rules_were_pre_declared_in": "config.yaml, hash " + cfg_hash},
        "cohort": {"n": len(pm), "sites": len(set(site)), "probes": int(Y.shape[1]),
                   "lusc": int(sub.sum()), "luad": int((1 - sub).sum())},
        "null": {"single_permutation_p95": thr, "single_permutation_median":
                 float(np.median(rho_null)), "m13_two_permutation_p95": thr_m13,
                 "m13_two_permutation_median": float(np.median(rho_null_m13style)),
                 "note": "m13 drew two independent permutations; this run uses one and "
                         "reports both so the threshold effect is visible"},
        "arms": {k: {"n_above_null": n[k], "frac_above_null": round(n[k] / len(rho_meth), 4),
                     "median_rho": float(np.median(v)), "max_rho": float(v.max())}
                 for k, v in (("meth", rho_meth), ("sub", rho_sub),
                              ("sublabel", rho_sublabel), ("globalmean", rho_global))},
        "m13_reference": {"n_above_null_subtype_emb": EXPECT["m13_n_above_null_subtype_emb"],
                          "n_above_null_subtype_label": EXPECT["m13_n_above_null_subtype_label"],
                          "reproduced_here": {"subtype_emb": int((m13_sub > m13_thr).sum()),
                                              "subtype_label": int((m13_lab > m13_thr).sum())}},
        "genomic_context": ctx,
        "context_spread": spread,
        "r10_reference_context": EXPECT["r10_context"],
        "circularity_control": {
            "baseline": "observed global mean methylation as a 1-d predictor",
            "n_above_null": n["globalmean"],
            "reading": secondary,
            "caveat": "each CpG contributes 1/%d to that mean, so self-inclusion is "
                      "negligible but not literally zero" % Y.shape[1]},
        "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir},
    }
    for nm, v in (("per_cpg_rho_meth_supervised.npy", rho_meth),
                  ("per_cpg_rho_subtype_supervised.npy", rho_sub),
                  ("per_cpg_rho_subtype_label.npy", rho_sublabel),
                  ("per_cpg_rho_globalmean.npy", rho_global),
                  ("per_cpg_rho_null_single_perm.npy", rho_null),
                  ("per_cpg_rho_null_m13style.npy", rho_null_m13style)):
        np.save(os.path.join(run.dir, nm), v.astype(np.float32))
    with open(os.path.join(run.dir, "probes_scanned.txt"), "w") as fh:
        fh.write("\n".join(pnames) + "\n")
    run.write("results.json", results)
    run.finalize()
    print("\nVERDICT:", verdict, "|", secondary)
    print("RUN DIR:", run.dir)


if __name__ == "__main__":
    main()
