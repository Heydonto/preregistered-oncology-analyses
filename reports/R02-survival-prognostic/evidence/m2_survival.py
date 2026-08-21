#!/usr/bin/env python3
"""R02 — does radiomics predict overall survival in UPENN-GBM? With positive controls.

Purpose is twofold:
  1. Test the prognostic endpoint, which unlike grade and IDH has ground truth here.
  2. Provide the POSITIVE CONTROLS that R01 lacked. Age and MGMT status are established
     prognostic factors in glioblastoma. If this pipeline recovers them and does not
     recover a radiomic signal, the R01 null is materially strengthened. If it fails to
     recover them, the machinery is suspect and R01 must be reconsidered.

Cohort note: every subject with a survival time in this collection is deceased, so the
data are UNCENSORED. Regression on log survival time is therefore valid and Cox is not
required for the primary; Cox is still used for the low-dimensional controls.
"""
import csv, hashlib, io, json, os, platform, re, subprocess, sys, time, zipfile
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
ZIP, CLIN = "/tmp/rf.zip", "/tmp/upenn.csv"
MASTER_SEED = 20260807
MODALITIES, ROIS = ["T1", "T1GD", "T2", "FLAIR"], ["ET", "ED", "NC"]

EXPECT = {
    "zip_bytes": 16121402,
    "n_complete_features": 599,
    "n_with_survival": 574,
    "all_deceased": True,
    "n_mgmt_definite": 247,
    "signal_rule_cindex_ci_low": 0.55,
}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Run:
    def __init__(self, tag):
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                             capture_output=True).stdout.decode().strip() or "nogit"
        self.dir = os.path.join(ROOT, "runs",
                               time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()) + f"-{tag}-{sha}")
        os.makedirs(self.dir, exist_ok=True)
        self.gates, self.logf = [], open(os.path.join(self.dir, "log.jsonl"), "a")

    def log(self, e, **kw):
        r = {"t": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "event": e, **kw}
        self.logf.write(json.dumps(r, default=str) + "\n"); self.logf.flush()
        print(f"[{r['t']}] {e}: " + " ".join(f"{k}={v}" for k, v in kw.items()), flush=True)

    def gate(self, name, expected, observed, ok, note=""):
        g = {"gate": name, "expected": expected, "observed": observed,
             "result": "PASS" if ok else "FAIL", "note": note,
             "t": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        self.gates.append(g); self.log("gate", gate=name, expected=expected,
                                       observed=observed, result=g["result"])
        with open(os.path.join(self.dir, "gates.json"), "w") as fh:
            json.dump(self.gates, fh, indent=2, default=str)
        if not ok:
            self.log("HALT", reason=f"gate {name} failed"); self.finalize(); sys.exit(2)

    def write(self, n, o):
        p = os.path.join(self.dir, n)
        with open(p, "w") as fh:
            json.dump(o, fh, indent=2, default=str) if n.endswith(".json") else fh.write(o)
        return p

    def finalize(self):
        L = []
        for dp, _, fs in os.walk(self.dir):
            for f in sorted(fs):
                if f != "MANIFEST.sha256":
                    fp = os.path.join(dp, f)
                    L.append(f"{sha256(fp)}  {os.path.relpath(fp, self.dir)}")
        with open(os.path.join(self.dir, "MANIFEST.sha256"), "w") as fh:
            fh.write("\n".join(sorted(L)) + "\n")


def cindex(time_obs, score_higher_is_longer):
    from lifelines.utils import concordance_index
    return float(concordance_index(time_obs, score_higher_is_longer))


def boot_ci(fn, n, rng, B=2000):
    vals = []
    for _ in range(B):
        i = rng.integers(0, n, n)
        try:
            vals.append(fn(i))
        except Exception:
            pass
    v = np.array(vals)
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def main():
    run = Run("R02surv")
    cfg = {
        "report": "R02", "milestone": "prognostic endpoint + positive controls",
        "primary": "radiomics -> overall survival, nested-CV Harrell C-index",
        "positive_controls": ["age -> survival", "MGMT status -> survival"],
        "signal_rule": "radiomics judged prognostic only if C-index 95% CI lower bound > 0.55",
        "control_rule": "if BOTH positive controls fail to recover a known association, "
                        "the pipeline is suspect and R01 must be reconsidered",
        "censoring": "none - all subjects with a survival time are deceased; "
                     "log-time regression is valid",
        "model": "elastic-net regression on log10(survival days), selection inside folds",
        "cv": "10x repeated 5-fold nested CV", "master_seed": MASTER_SEED,
        "expectations": EXPECT, "relates_to": "R01 (MGMT classification, null)",
    }
    cfg_hash = sha256(run.write("config.yaml", json.dumps(cfg, indent=2)))
    run.log("config_written", sha256=cfg_hash)
    run.write("env.txt", subprocess.run([sys.executable, "-m", "pip", "freeze"],
                                        capture_output=True).stdout.decode())

    zb = os.path.getsize(ZIP)
    run.write("inputs.json", {ZIP: {"bytes": zb, "sha256": sha256(ZIP)},
                              CLIN: {"bytes": os.path.getsize(CLIN), "sha256": sha256(CLIN)}})
    run.gate("G0_integrity", EXPECT["zip_bytes"], zb, zb == EXPECT["zip_bytes"])

    # ---- features: complete cases across all 12 files (same rule as R01) ------
    z = zipfile.ZipFile(ZIP)
    tables, feat_names, per_file = {}, [], {}
    for mod in MODALITIES:
        for roi in ROIS:
            with z.open(f"Radiomic_Features_CaPTk_automaticsegm_{mod}_{roi}.csv") as fh:
                rows = list(csv.reader(io.TextIOWrapper(fh, "utf8")))
            cols = [f"{mod}_{roi}::{c}" for c in rows[0][1:]]
            feat_names += cols
            pres = set()
            for r in rows[1:]:
                if not r:
                    continue
                pres.add(r[0].strip())
                tables.setdefault(r[0].strip(), {}).update(
                    {c: (float(v) if v not in ("", "NA", "NaN") else np.nan)
                     for c, v in zip(cols, r[1:])})
            per_file[f"{mod}_{roi}"] = pres
    complete = set.intersection(*per_file.values())
    run.gate("G1_complete_features", EXPECT["n_complete_features"], len(complete),
             len(complete) == EXPECT["n_complete_features"], "same complete-case rule as R01")

    # ---- outcome -------------------------------------------------------------
    def num(v):
        try:
            return float(v)
        except Exception:
            return None
    clin = {}
    for r in csv.DictReader(open(CLIN)):
        clin[r["ID"].strip()] = (num(r["Survival_from_surgery_days_UPDATED"]),
                                 r["Survival_Status"].strip(),
                                 num(r["Age_at_scan_years"]), r["MGMT"].strip())
    ids = sorted(i for i in complete if i in clin and clin[i][0] is not None and clin[i][0] > 0)
    run.gate("G2_survival_n", EXPECT["n_with_survival"], len(ids),
             len(ids) == EXPECT["n_with_survival"], "complete features AND survival time > 0")
    status = {clin[i][1] for i in ids}
    run.gate("G3_no_censoring", "all Deceased", sorted(status), status == {"Deceased"},
             "uncensored cohort validates log-time regression")

    t = np.array([clin[i][0] for i in ids], float)
    logt = np.log10(t)
    age = np.array([clin[i][2] for i in ids], float)
    X = np.array([[tables[i].get(c, np.nan) for c in feat_names] for i in ids], float)
    miss = np.isnan(X).mean(0); keep = miss <= 0.10
    X = X[:, keep]
    X = np.where(np.isnan(X), np.nanmedian(X, 0), X)
    run.gate("G4_features", 1728, int(keep.sum()), int(keep.sum()) == 1728,
             "no feature exceeded 10% missingness")
    run.write("cohort.json", {"n": len(ids), "subject_ids": ids,
                              "median_survival_days": float(np.median(t)),
                              "iqr": [float(np.percentile(t, 25)), float(np.percentile(t, 75))]})

    # ---- POSITIVE CONTROL 1: age --------------------------------------------
    rng = np.random.default_rng(MASTER_SEED)
    c_age = cindex(t, -age)  # older -> shorter survival, so negate
    lo, hi = boot_ci(lambda i: cindex(t[i], -age[i]), len(t), rng)
    from scipy.stats import spearmanr
    rho_age, p_age = spearmanr(age, t)
    run.log("control_age", cindex=round(c_age, 4), ci=[round(lo, 4), round(hi, 4)],
            spearman=round(float(rho_age), 4), p=f"{p_age:.2e}")

    # ---- POSITIVE CONTROL 2: MGMT -------------------------------------------
    mg_ids = [k for k, i in enumerate(ids) if clin[ids[k]][3] in ("Methylated", "Unmethylated")]
    mg = np.array([1 if clin[ids[k]][3] == "Methylated" else 0 for k in mg_ids])
    tm = t[mg_ids]
    run.gate("G5_mgmt_subset", EXPECT["n_mgmt_definite"], len(mg_ids),
             len(mg_ids) == EXPECT["n_mgmt_definite"], "definite MGMT within the survival cohort")
    from lifelines.statistics import logrank_test
    lr = logrank_test(tm[mg == 1], tm[mg == 0],
                      np.ones((mg == 1).sum()), np.ones((mg == 0).sum()))
    med_m, med_u = float(np.median(tm[mg == 1])), float(np.median(tm[mg == 0]))
    run.log("control_mgmt", n=len(mg_ids), median_methylated=med_m, median_unmethylated=med_u,
            logrank_p=f"{lr.p_value:.4g}")

    controls = {
        "age": {"cindex": c_age, "ci95": [lo, hi], "spearman_rho": float(rho_age),
                "spearman_p": float(p_age), "direction": "older = shorter survival",
                "recovered": bool(lo > 0.55)},
        "mgmt": {"n": len(mg_ids), "median_days_methylated": med_m,
                 "median_days_unmethylated": med_u, "logrank_p": float(lr.p_value),
                 "direction_expected": "methylated = longer survival",
                 "direction_observed": "methylated longer" if med_m > med_u else "methylated shorter",
                 "recovered": bool(lr.p_value < 0.05 and med_m > med_u)},
    }
    run.gate("G6_positive_controls",
             "at least one known prognostic association recovered",
             {"age_recovered": controls["age"]["recovered"],
              "mgmt_recovered": controls["mgmt"]["recovered"]},
             controls["age"]["recovered"] or controls["mgmt"]["recovered"],
             "if both fail, the pipeline cannot be trusted to detect real signal")

    # ---- PRIMARY: radiomics -> survival -------------------------------------
    from sklearn.model_selection import RepeatedKFold, GridSearchCV, KFold
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.feature_selection import VarianceThreshold
    from sklearn.linear_model import ElasticNet

    pipe = Pipeline([("var", VarianceThreshold(1e-8)), ("scale", StandardScaler()),
                     ("reg", ElasticNet(max_iter=5000, random_state=MASTER_SEED))])
    grid = {"reg__alpha": [0.01, 0.1, 1.0], "reg__l1_ratio": [0.1, 0.5, 0.9]}
    rkf = RepeatedKFold(n_splits=5, n_repeats=10, random_state=MASTER_SEED)
    oof = np.full((10, len(ids)), np.nan)
    for k, (tr, te) in enumerate(rkf.split(X)):
        rep = k // 5
        gs = GridSearchCV(pipe, grid, scoring="neg_mean_squared_error",
                          cv=KFold(5, shuffle=True, random_state=MASTER_SEED + rep), n_jobs=-1)
        gs.fit(X[tr], logt[tr])
        oof[rep, te] = gs.predict(X[te])
        if k % 10 == 0:
            run.log("cv_progress", fold=k + 1, of=50)
    pred = oof.mean(0)
    c_rad = cindex(t, pred)  # higher predicted log-time -> longer survival
    rng2 = np.random.default_rng(MASTER_SEED + 1)
    rlo, rhi = boot_ci(lambda i: cindex(t[i], pred[i]), len(t), rng2)
    rho_r, p_r = spearmanr(pred, logt)
    np.savetxt(os.path.join(run.dir, "oof_survival.csv"),
               np.column_stack([t, logt, pred]), delimiter=",",
               header="survival_days,log10_days,predicted_log10_days", comments="")

    signal = bool(rlo > EXPECT["signal_rule_cindex_ci_low"])
    results = {
        "primary_radiomics": {"cindex": c_rad, "ci95": [rlo, rhi],
                              "spearman_rho_pred_vs_obs": float(rho_r),
                              "spearman_p": float(p_r),
                              "prognostic_by_prespecified_rule": signal},
        "positive_controls": controls,
        "_decision": {"rule": f"C-index CI lower > {EXPECT['signal_rule_cindex_ci_low']}",
                      "observed_ci_low": rlo,
                      "decision": "radiomics prognostic" if signal
                                  else "radiomics NOT prognostic by pre-specified rule",
                      "prespecified_in_config_sha256": cfg_hash},
        "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir, "n": len(ids),
                        "n_features": int(X.shape[1])},
    }
    run.write("results.json", results)
    run.log("decision", **results["_decision"])
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
