#!/usr/bin/env python3
"""R07 — does DNA methylation predict EGFR-TKI response? GSE147377, n=79.

This is the only cohort in our holdings pairing genome-wide methylation with a treatment
response endpoint, and it is the question the withdrawn lung manuscript was actually about.
It is also small: 40 partial responders, 29 progressors, 10 stable. Declared underpowered
in advance.

A structural risk must be handled first. The beta matrix columns are named EXP_1..EXP_79,
not GSM accessions, so the mapping to phenotype rests on assuming series-matrix column order.
A mis-mapping would invalidate everything. The positive controls therefore double as a
mapping test:

  CONTROL A (named, directional): cg05575921 in AHRR is the best-established smoking
      methylation biomarker; smokers show LOWER methylation. If the mapping is wrong this
      association disappears.
  CONTROL B (structural): sex is near-perfectly predictable from methylation via X-linked
      probes. With a correct mapping some probes separate the sexes almost completely.

If neither control holds, the mapping is wrong and no response result is reported.
"""
import csv, gzip, json, os, sys
import numpy as np
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_kit import Run, bh

D = "/Users/rezanehzati/quantara-staging/staged/labels/lung/GSE147377/"
BETA = D + "suppl/GSE147377_AverageBeta_Matrix.csv.gz"
MATRIX = D + "matrix/GSE147377_series_matrix.txt.gz"
SEED = 20260807
AHRR = "cg05575921"
EXPECT = {"n_samples": 79, "pr": 40, "pd": 29, "sd": 10,
          "control_a_direction": "smokers LOWER at cg05575921"}


def phenotype():
    """Sample-order phenotype from the series matrix."""
    fields = {}
    gsm = None
    for line in gzip.open(MATRIX, "rt", errors="ignore"):
        if line.startswith("!Sample_geo_accession"):
            gsm = [v.strip().strip('"') for v in line.split("\t")[1:]]
        if line.startswith("!Sample_characteristics_ch1"):
            vals = [v.strip().strip('"') for v in line.split("\t")[1:]]
            key = vals[0].split(":")[0].strip() if ":" in vals[0] else "field"
            fields[key] = [v.split(": ", 1)[-1].strip() if ": " in v else "" for v in vals]
        if line.startswith("!series_matrix_table_begin"):
            break
    return gsm, fields


def main():
    run = Run("R07tki")
    cfg = {
        "report": "R07",
        "question": "does genome-wide methylation predict EGFR-TKI response?",
        "cohort": "GSE147377, 79 advanced LUAD, all TKI-treated, Illumina 450k",
        "endpoint": "partial response (PR) vs progressive disease (PD); SD excluded as ambiguous",
        "declared_limitation": "underpowered by design: ~69 patients against ~450k probes",
        "structural_risk": "beta matrix columns are EXP_1..EXP_79, not GSM ids; phenotype "
                           "mapping assumes series-matrix column order",
        "control_A": "cg05575921 (AHRR) vs smoking; smokers expected LOWER",
        "control_B": "sex separability from X-linked probes",
        "control_rule": "if neither control holds, the mapping is wrong and no response "
                        "result is reported",
        "model": "elastic-net logistic on the top-variance probes, selection inside folds",
        "cv": "10x repeated 5-fold nested CV", "seed": SEED, "expectations": EXPECT,
    }
    cfg_hash = run.start(cfg, [BETA, MATRIX])

    gsm, fields = phenotype()
    run.gate("G0_samples", EXPECT["n_samples"], len(gsm), len(gsm) == EXPECT["n_samples"])
    resp = fields.get("tki_response", [])
    import collections
    rc = collections.Counter(resp)
    run.gate("G1_response_labels", {"PR": EXPECT["pr"], "PD": EXPECT["pd"], "SD": EXPECT["sd"]},
             dict(rc), rc.get("PR") == EXPECT["pr"] and rc.get("PD") == EXPECT["pd"])

    # ---- load beta values (AVG_Beta columns only) ----
    fh = gzip.open(BETA, "rt", errors="ignore")
    hdr = fh.readline().rstrip("\n").split(",")
    bidx = [i for i, c in enumerate(hdr) if c.endswith(".AVG_Beta")]
    run.gate("G2_beta_columns", EXPECT["n_samples"], len(bidx), len(bidx) == EXPECT["n_samples"],
             "one AVG_Beta column per sample, in EXP_n order")
    probes, rows = [], []
    ahrr = None
    for line in fh:
        p = line.rstrip("\n").split(",")
        vals = np.array([float(p[i]) if p[i] not in ("", "NA", "NaN") else np.nan
                         for i in bidx])
        if p[0] == AHRR:
            ahrr = vals
        probes.append(p[0])
        rows.append(vals)
    B = np.array(rows, dtype=np.float32)
    run.log("beta_loaded", probes=len(probes), samples=B.shape[1])
    run.gate("G3_ahrr_present", f"{AHRR} present", ahrr is not None, ahrr is not None)

    # ---- CONTROL A: AHRR vs smoking ----
    smoke = np.array(fields.get("smoking", []), object)
    ok = np.isin(smoke, ["Never-smoker", "smoker"]) & ~np.isnan(ahrr)
    a = ahrr[ok & (smoke == "smoker")]
    b = ahrr[ok & (smoke == "Never-smoker")]
    U, p_a = mannwhitneyu(a, b)
    lower_in_smokers = bool(np.median(a) < np.median(b))
    ctrl_a = bool(p_a < 0.05 and lower_in_smokers)
    run.log("control_A_AHRR", p=f"{p_a:.4g}", median_smoker=round(float(np.median(a)), 4),
            median_never=round(float(np.median(b)), 4), n=int(ok.sum()),
            direction_correct=lower_in_smokers, passes=ctrl_a,
            note="direction is the mapping evidence; significance is power-limited "
                 "(this cohort is 55/79 never-smokers)")

    # ---- CONTROL B: sex separability ----
    sex = np.array(fields.get("gender", []), object)
    oks = np.isin(sex, ["Male", "Female"])
    male = sex == "Male"
    # best single-probe separation between sexes
    from sklearn.metrics import roc_auc_score
    var = np.nanvar(B, axis=1)
    cand = np.argsort(var)[-20000:]
    best, bestp = 0.5, None
    for i in cand:
        v = B[i][oks]
        if np.isnan(v).any():
            continue
        auc = roc_auc_score(male[oks].astype(int), v)
        auc = max(auc, 1 - auc)
        if auc > best:
            best, bestp = auc, probes[i]
    ctrl_b = bool(best > 0.95)
    run.log("control_B_sex", best_auc=round(best, 4), probe=bestp, passes=ctrl_b)
    # PRE-SPECIFIED RULE, applied as written: halt only if NEITHER control holds.
    run.gate("G4_mapping_controls", "at least one of {AHRR-smoking, sex-separability} holds",
             {"A_AHRR": {"p": round(float(p_a), 4), "direction_correct": lower_in_smokers,
                         "passes": ctrl_a},
              "B_sex": {"best_auc": round(best, 4), "passes": ctrl_b}},
             ctrl_a or ctrl_b,
             "config rule: 'if neither control holds, the mapping is wrong'. An earlier "
             "implementation gated on control A alone and halted incorrectly.")

    # ---- PRIMARY: PR vs PD ----
    sel = np.isin(resp, ["PR", "PD"])
    y = (np.array(resp, object)[sel] == "PD").astype(int)   # 1 = progression
    X = B[:, sel].T.astype(np.float64)
    keep = ~np.isnan(X).any(axis=0)
    X = X[:, keep]
    pnames = [p for p, k in zip(probes, keep) if k]
    run.log("primary_cohort", n=int(sel.sum()), progressors=int(y.sum()),
            probes_complete=int(keep.sum()))

    from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV, StratifiedKFold
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.feature_selection import SelectKBest, f_classif
    from sklearn.linear_model import LogisticRegression
    pipe = Pipeline([("sel", SelectKBest(f_classif, k=2000)), ("sc", StandardScaler()),
                     ("clf", LogisticRegression(penalty="elasticnet", solver="saga",
                                                max_iter=4000, tol=1e-3, random_state=SEED))])
    grid = {"clf__C": [0.01, 0.1, 1.0], "clf__l1_ratio": [0.5]}
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=SEED)
    oof = np.full((10, len(y)), np.nan)
    for k, (tr, te) in enumerate(rskf.split(X, y)):
        rep = k // 5
        gs = GridSearchCV(pipe, grid, scoring="roc_auc",
                          cv=StratifiedKFold(3, shuffle=True, random_state=SEED + rep),
                          n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[rep, te] = gs.predict_proba(X[te])[:, 1]
        if k % 10 == 0:
            run.log("cv_progress", fold=k + 1, of=50)
    pred = oof.mean(0)
    auc = float(roc_auc_score(y, pred))
    per_rep = [float(roc_auc_score(y, oof[r])) for r in range(10)]

    # permutation null
    rng = np.random.default_rng(SEED)
    null = []
    fixed = Pipeline([("sel", SelectKBest(f_classif, k=2000)), ("sc", StandardScaler()),
                      ("clf", LogisticRegression(penalty="elasticnet", solver="saga", C=0.1,
                                                 l1_ratio=0.5, max_iter=2000, tol=1e-3,
                                                 random_state=SEED))])
    for i in range(100):
        yp = rng.permutation(y)
        o = np.full(len(y), np.nan)
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED + i).split(X, yp):
            fixed.fit(X[tr], yp[tr])
            o[te] = fixed.predict_proba(X[te])[:, 1]
        null.append(roc_auc_score(yp, o))
        if i % 25 == 0:
            run.log("perm_progress", done=i + 1, of=100)
    null = np.array(null)
    pval = float((np.sum(null >= auc) + 1) / 101)
    run.gate("G6_negative_control", "permutation null contains 0.50",
             {"null_mean": round(float(null.mean()), 4),
              "p2.5": round(float(np.percentile(null, 2.5)), 4),
              "p97.5": round(float(np.percentile(null, 97.5)), 4)},
             bool(np.percentile(null, 2.5) <= 0.5 <= np.percentile(null, 97.5)),
             f"100 permutations; observed p={pval:.4f}")

    np.savetxt(os.path.join(run.dir, "oof_tki_response.csv"),
               np.column_stack([y, pred]), delimiter=",",
               header="progressive_disease,mean_oof_probability", comments="")
    results = {
        "control_A_AHRR_smoking": {"p": float(p_a), "median_smoker": float(np.median(a)),
                                   "median_never": float(np.median(b)),
                                   "lower_in_smokers": lower_in_smokers, "recovered": ctrl_a},
        "control_B_sex": {"best_probe_auc": best, "probe": bestp, "recovered": ctrl_b},
        "primary_tki_response": {"n": int(sel.sum()), "n_progressors": int(y.sum()),
                                 "auroc": auc,
                                 "per_repeat_mean": float(np.mean(per_rep)),
                                 "per_repeat_sd": float(np.std(per_rep, ddof=1)),
                                 "permutation_p": pval,
                                 "null_mean": float(null.mean())},
        "_decision": {"interpretation":
                      "methylation predicts TKI response" if pval < 0.05 else
                      "no evidence that methylation predicts TKI response at this sample size",
                      "prespecified_in_config_sha256": cfg_hash},
        "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir,
                        "probes_used": int(keep.sum())},
    }
    run.write("results.json", results)
    run.write("permutation_null.json", {"n_perm": 100, "null": null.tolist(), "observed": auc})
    run.log("decision", **results["_decision"])
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
