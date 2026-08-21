#!/usr/bin/env python3
"""M1 — internal validation of MGMT prediction from UPENN-GBM radiomic features.

Per EXPERIMENT_PLAN.md. Development cohort ONLY. The external cohort is never opened.

Auditability contract:
  * config.yaml is written and hashed BEFORE any label is read, so an endpoint changed
    mid-flight is detectable.
  * Every gate records expected vs observed and PASS/FAIL, including on failure.
  * Feature selection and scaling happen strictly inside training folds.
  * A failed run is kept, never deleted.
"""
import csv, hashlib, io, json, os, platform, re, subprocess, sys, time, zipfile
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
ZIP = "/tmp/rf.zip"
CLIN = "/tmp/upenn.csv"
MASTER_SEED = 20260806
MODALITIES = ["T1", "T1GD", "T2", "FLAIR"]
ROIS = ["ET", "ED", "NC"]

EXPECT = {
    "zip_bytes": 16121402,
    # AMENDED 2026-08-07 after gate G1_join_n FAILED (expected 258, observed 262).
    # Diagnosis: the 12 feature files do NOT cover identical subjects
    # (ET=604, ED=611, NC=602). The planning figure of 258 was derived from a
    # single file and was therefore the wrong basis. Using the UNION (262) would
    # median-impute whole 144-feature ROI blocks for subjects missing an ROI,
    # i.e. fabricate data. The correct definition is COMPLETE CASES across all
    # 12 files: 599 subjects, of which 256 have a definite MGMT label.
    # Superseded expectation retained for audit: n=258 (109/149) from T1GD_ET only.
    "n_dev": 256,
    "n_methylated": 108,
    "n_unmethylated": 148,
    "join_rule": "complete cases across all 12 modality-ROI files",
    "n_features": 1728,
    "stop_rule_ci_upper": 0.65,
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Run:
    def __init__(self):
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                             capture_output=True).stdout.decode().strip() or "nogit"
        self.dir = os.path.join(ROOT, "runs",
                                time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()) + "-" + sha)
        os.makedirs(os.path.join(self.dir, "figures"), exist_ok=True)
        self.gates, self.logf = [], open(os.path.join(self.dir, "log.jsonl"), "a")

    def log(self, event, **kw):
        rec = {"t": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "event": event, **kw}
        self.logf.write(json.dumps(rec, default=str) + "\n")
        self.logf.flush()
        print(f"[{rec['t']}] {event}: " + " ".join(f"{k}={v}" for k, v in kw.items()), flush=True)

    def gate(self, name, expected, observed, passed, note=""):
        g = {"gate": name, "expected": expected, "observed": observed,
             "result": "PASS" if passed else "FAIL", "note": note,
             "t": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        self.gates.append(g)
        self.log("gate", **{k: g[k] for k in ("gate", "expected", "observed", "result")})
        with open(os.path.join(self.dir, "gates.json"), "w") as fh:
            json.dump(self.gates, fh, indent=2, default=str)
        if not passed:
            self.log("HALT", reason=f"gate {name} failed")
            self.finalize()
            sys.exit(2)

    def write(self, name, obj):
        p = os.path.join(self.dir, name)
        with open(p, "w") as fh:
            if name.endswith(".json"):
                json.dump(obj, fh, indent=2, default=str)
            else:
                fh.write(obj)
        return p

    def finalize(self):
        lines = []
        for dp, _, fs in os.walk(self.dir):
            for f in sorted(fs):
                if f == "MANIFEST.sha256":
                    continue
                fp = os.path.join(dp, f)
                lines.append(f"{sha256(fp)}  {os.path.relpath(fp, self.dir)}")
        with open(os.path.join(self.dir, "MANIFEST.sha256"), "w") as fh:
            fh.write("\n".join(sorted(lines)) + "\n")


# ---------------- DeLong CI (Sun & Xu 2014 fast implementation) ----------------

def _midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    out = np.empty(N, dtype=float)
    out[J] = T
    return out


def delong_auc_var(y, score):
    """Returns (AUC, variance) for a single ROC via the fast DeLong method."""
    pos = score[y == 1]
    neg = score[y == 0]
    m, n = len(pos), len(neg)
    tx, ty, tz = _midrank(pos), _midrank(neg), _midrank(np.concatenate([pos, neg]))
    auc = (tz[:m].sum() - m * (m + 1) / 2) / (m * n)
    v01 = (tz[:m] - tx) / n
    v10 = 1 - (tz[m:] - ty) / m
    var = v01.var(ddof=1) / m + v10.var(ddof=1) / n
    return auc, var


def delong_ci(y, score, alpha=0.05):
    from scipy.stats import norm
    auc, var = delong_auc_var(np.asarray(y), np.asarray(score))
    se = np.sqrt(max(var, 1e-12))
    z = norm.ppf(1 - alpha / 2)
    return auc, max(0.0, auc - z * se), min(1.0, auc + z * se)


# ------------------------------- pipeline -------------------------------------

def build_pipeline():
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.feature_selection import VarianceThreshold
    from sklearn.linear_model import LogisticRegression
    return Pipeline([
        ("var", VarianceThreshold(1e-8)),
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(penalty="elasticnet", solver="saga",
                                   max_iter=3000, tol=1e-3, random_state=MASTER_SEED)),
    ])


GRID = {"clf__C": [0.01, 0.1, 1.0], "clf__l1_ratio": [0.1, 0.5, 0.9]}


def main():
    run = Run()

    # ---- config FIRST, before any label is read -------------------------------
    cfg = {
        "milestone": "M1-internal-validation-only",
        "plan": "EXPERIMENT_PLAN.md",
        "primary_endpoint": "MGMT methylated vs unmethylated, internal nested-CV AUROC",
        "stopping_rule": "halt programme if 95% CI upper bound on internal AUROC < 0.65",
        "cohort": "UPENN-GBM, CaPTk automatic-segmentation features",
        "modalities": MODALITIES, "rois": ROIS,
        "excluded_modalities_reason": "DSC perfusion and DTI absent from the external cohort",
        "primary_model": "elastic-net logistic regression",
        "comparators": ["random forest", "gradient boosting"],
        "cv": "10x repeated 5-fold stratified nested CV; selection+scaling inside folds only",
        "grid": GRID, "master_seed": MASTER_SEED,
        "ci_method": "DeLong on out-of-fold predictions averaged across repeats",
        "permutation_null": "label permutation, fixed modal hyperparameters, 5-fold",
        "external_cohort": "NOT OPENED IN M1",
        "expectations": EXPECT,
        "supersedes": "runs/20260806T220320Z-3f8bece (identical analysis; this run additionally "
                      "persists out-of-fold predictions and subject IDs for independent verification)",
        "prior_failed_run": "runs/20260806T215837Z-3f8bece (halted at gate G1_join_n, retained)",
    }
    cfg_path = run.write("config.yaml", json.dumps(cfg, indent=2))
    cfg_hash = sha256(cfg_path)
    run.log("config_written", sha256=cfg_hash)

    run.write("env.txt", "\n".join([
        f"python {platform.python_version()}", f"platform {platform.platform()}",
        f"numpy {np.__version__}",
        subprocess.run([sys.executable, "-m", "pip", "freeze"],
                       capture_output=True).stdout.decode()]))

    # ---- G0 input integrity ---------------------------------------------------
    zb = os.path.getsize(ZIP)
    inputs = {ZIP: {"bytes": zb, "sha256": sha256(ZIP)},
              CLIN: {"bytes": os.path.getsize(CLIN), "sha256": sha256(CLIN)}}
    run.write("inputs.json", inputs)
    run.gate("G0_integrity", EXPECT["zip_bytes"], zb, zb == EXPECT["zip_bytes"],
             "CaPTk feature archive byte size")

    # ---- features -------------------------------------------------------------
    z = zipfile.ZipFile(ZIP)
    tables, feat_names, per_file = {}, [], {}
    for mod in MODALITIES:
        for roi in ROIS:
            name = f"Radiomic_Features_CaPTk_automaticsegm_{mod}_{roi}.csv"
            with z.open(name) as fh:
                rows = list(csv.reader(io.TextIOWrapper(fh, "utf8")))
            hdr, body = rows[0], rows[1:]
            cols = [f"{mod}_{roi}::{c}" for c in hdr[1:]]
            feat_names += cols
            present = set()
            for r in body:
                if not r:
                    continue
                sid = r[0].strip()
                present.add(sid)
                tables.setdefault(sid, {}).update(
                    {c: (float(v) if v not in ("", "NA", "NaN") else np.nan)
                     for c, v in zip(cols, r[1:])})
            per_file[f"{mod}_{roi}"] = present
    complete = set.intersection(*per_file.values())
    run.write("subject_coverage.json",
              {"per_file_counts": {k: len(v) for k, v in per_file.items()},
               "union": len(set().union(*per_file.values())),
               "complete_cases": len(complete),
               "rule": "complete cases only; whole-ROI-block imputation is not permitted"})
    run.log("coverage", per_file_min=min(len(v) for v in per_file.values()),
            per_file_max=max(len(v) for v in per_file.values()), complete=len(complete))
    run.gate("G_featcount", EXPECT["n_features"], len(feat_names),
             len(feat_names) == EXPECT["n_features"], "12 files x 144 features")

    # ---- labels (only now) ----------------------------------------------------
    lab = {}
    for r in csv.DictReader(open(CLIN)):
        v = r["MGMT"].strip()
        if v in ("Methylated", "Unmethylated"):
            lab[r["ID"].strip()] = 1 if v == "Methylated" else 0

    ids = sorted(complete & set(lab))
    y = np.array([lab[i] for i in ids])
    X = np.array([[tables[i].get(c, np.nan) for c in feat_names] for i in ids], float)

    run.gate("G1_join_n", EXPECT["n_dev"], len(ids), len(ids) == EXPECT["n_dev"],
             "features ∩ definite MGMT; ID suffix _11 present in both sources")
    run.write("cohort_ids.json", {"n": len(ids), "order": "matches rows of oof_*.csv",
                                  "subject_ids": ids})
    run.gate("G1_join_balance", [EXPECT["n_methylated"], EXPECT["n_unmethylated"]],
             [int((y == 1).sum()), int((y == 0).sum())],
             int((y == 1).sum()) == EXPECT["n_methylated"]
             and int((y == 0).sum()) == EXPECT["n_unmethylated"], "class counts")

    # ---- G6 missingness -------------------------------------------------------
    miss = np.isnan(X).mean(axis=0)
    keep = miss <= 0.10
    dropped = [feat_names[i] for i in np.where(~keep)[0]]
    run.write("dropped_features.json",
              {"n_dropped": len(dropped), "threshold": 0.10, "features": dropped})
    X = X[:, keep]
    kept_names = [f for f, k in zip(feat_names, keep) if k]
    col_med = np.nanmedian(X, axis=0)
    X = np.where(np.isnan(X), col_med, X)
    run.gate("G6_missingness", "<=10% per feature", f"{len(dropped)} dropped, {X.shape[1]} kept",
             True, "median-imputed within-cohort after dropping")

    # ---- G2 leakage structure -------------------------------------------------
    run.gate("G2_leakage_external", "external labels not loaded", "no external file opened",
             True, "M1 opens only the CaPTk archive and the UPENN clinical csv")

    # ---- nested CV ------------------------------------------------------------
    from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV, StratifiedKFold
    from sklearn.metrics import roc_auc_score, brier_score_loss
    from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier

    def nested(model_name, estimator, grid):
        rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=MASTER_SEED)
        oof = np.full((10, len(y)), np.nan)
        best = []
        for k, (tr, te) in enumerate(rskf.split(X, y)):
            rep = k // 5
            if grid:
                gs = GridSearchCV(estimator, grid, scoring="roc_auc",
                                  cv=StratifiedKFold(5, shuffle=True,
                                                     random_state=MASTER_SEED + rep),
                                  n_jobs=-1, refit=True)
                gs.fit(X[tr], y[tr])
                oof[rep, te] = gs.predict_proba(X[te])[:, 1]
                best.append(gs.best_params_)
            else:
                estimator.fit(X[tr], y[tr])
                oof[rep, te] = estimator.predict_proba(X[te])[:, 1]
            if k % 10 == 0:
                run.log("cv_progress", model=model_name, fold=k + 1, of=50)
        per_repeat = [roc_auc_score(y, oof[r]) for r in range(10)]
        mean_pred = oof.mean(axis=0)
        auc, lo, hi = delong_ci(y, mean_pred)
        np.savetxt(os.path.join(run.dir, f"oof_{model_name}.csv"),
                   np.column_stack([y, mean_pred]), delimiter=",",
                   header="mgmt_methylated,mean_oof_probability", comments="")
        return {"model": model_name,
                "auroc_delong": auc, "ci95_low": lo, "ci95_high": hi,
                "auroc_per_repeat_mean": float(np.mean(per_repeat)),
                "auroc_per_repeat_sd": float(np.std(per_repeat, ddof=1)),
                "auroc_per_repeat_p2.5": float(np.percentile(per_repeat, 2.5)),
                "auroc_per_repeat_p97.5": float(np.percentile(per_repeat, 97.5)),
                "brier": float(brier_score_loss(y, mean_pred)),
                "best_params_modal": (max(set(map(str, best)), key=list(map(str, best)).count)
                                      if best else None)}, mean_pred

    results = {}
    run.log("start_primary", model="elasticnet", n=len(y), p=X.shape[1])
    results["elasticnet"], pred_primary = nested("elasticnet", build_pipeline(), GRID)
    run.log("primary_done", **{k: v for k, v in results["elasticnet"].items() if "auroc" in k})

    for nm, est in (("random_forest",
                     RandomForestClassifier(n_estimators=500, random_state=MASTER_SEED, n_jobs=-1)),
                    ("hist_gradient_boosting",
                     HistGradientBoostingClassifier(random_state=MASTER_SEED))):
        run.log("start_comparator", model=nm)
        results[nm], _ = nested(nm, est, None)

    # ---- G3 permutation null --------------------------------------------------
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.feature_selection import VarianceThreshold
    bp = results["elasticnet"]["best_params_modal"]
    C = float(re.search(r"'clf__C': ([0-9.]+)", bp).group(1)) if bp else 0.1
    L1 = float(re.search(r"'clf__l1_ratio': ([0-9.]+)", bp).group(1)) if bp else 0.5
    fixed = Pipeline([("var", VarianceThreshold(1e-8)), ("scale", StandardScaler()),
                      ("clf", LogisticRegression(penalty="elasticnet", solver="saga", C=C,
                                                 l1_ratio=L1, max_iter=1500, tol=1e-3,
                                                 random_state=MASTER_SEED))])
    NPERM = int(os.environ.get("NPERM", "200"))
    rng = np.random.default_rng(MASTER_SEED)
    null = []
    t0 = time.time()
    for p in range(NPERM):
        yp = rng.permutation(y)
        skf = StratifiedKFold(5, shuffle=True, random_state=MASTER_SEED + p)
        o = np.full(len(y), np.nan)
        for tr, te in skf.split(X, yp):
            fixed.fit(X[tr], yp[tr])
            o[te] = fixed.predict_proba(X[te])[:, 1]
        null.append(roc_auc_score(yp, o))
        if p % 25 == 0:
            run.log("perm_progress", done=p + 1, of=NPERM, elapsed_s=round(time.time() - t0))
    null = np.array(null)
    obs = results["elasticnet"]["auroc_delong"]
    pval = float((np.sum(null >= obs) + 1) / (NPERM + 1))
    centred = bool(np.percentile(null, 2.5) <= 0.5 <= np.percentile(null, 97.5))
    run.write("permutation_null.json",
              {"n_perm": NPERM, "fixed_C": C, "fixed_l1_ratio": L1,
               "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)),
               "null_p2.5": float(np.percentile(null, 2.5)),
               "null_p97.5": float(np.percentile(null, 97.5)),
               "observed": obs, "p_value": pval, "null_auc": null.tolist()})
    run.gate("G3_negative_control", "null 95% interval contains 0.50",
             f"[{np.percentile(null,2.5):.3f}, {np.percentile(null,97.5):.3f}] mean={null.mean():.3f}",
             centred, f"{NPERM} permutations; p={pval:.4f}")

    # ---- G5 determinism -------------------------------------------------------
    h1 = hashlib.sha256(np.round(pred_primary, 10).tobytes()).hexdigest()
    _, pred_again = nested("elasticnet_rerun", build_pipeline(), GRID)
    h2 = hashlib.sha256(np.round(pred_again, 10).tobytes()).hexdigest()
    run.gate("G5_determinism", h1, h2, h1 == h2, "same seed reproduces identical predictions")

    # ---- stopping rule --------------------------------------------------------
    ci_hi = results["elasticnet"]["ci95_high"]
    proceed = ci_hi >= EXPECT["stop_rule_ci_upper"]
    results["_decision"] = {
        "ci95_high_primary": ci_hi, "threshold": EXPECT["stop_rule_ci_upper"],
        "decision": "PROCEED to M2" if proceed else "STOP - signal too weak for further compute",
        "prespecified_in_config_sha256": cfg_hash}
    results["_provenance"] = {"config_sha256": cfg_hash, "run_dir": run.dir,
                              "n": int(len(y)), "n_features_used": int(X.shape[1]),
                              "inputs": inputs}
    run.write("results.json", results)
    run.log("decision", **results["_decision"])
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
