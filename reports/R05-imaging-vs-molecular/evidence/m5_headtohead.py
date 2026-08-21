#!/usr/bin/env python3
"""R05 — imaging versus molecular markers, head to head in the SAME patients.

R04's audit identified this as the analysis that settles the comparison. R02 measured
imaging on UPENN-GBM and R04 measured molecular markers on TCGA/MSK, so the two were matched
on disease but not on population. UPENN carries IDH1 and MGMT for the same patients that
have radiomic features and survival, so the comparison can be made within one cohort.

Pre-specified questions:
  Q1 POSITIVE CONTROL - is MGMT prognostic in this subset? (R02 already showed p=5.8e-8)
  Q2 C-index of molecular markers alone (IDH1 + MGMT), cross-validated
  Q3 C-index of imaging alone, on the same patients, cross-validated
  Q4 C-index of imaging + molecular combined; does imaging add beyond molecular?
  Decision rule: imaging adds independent prognostic value only if the likelihood-ratio test
  for adding the radiomic score to a molecular Cox model gives p < 0.05.
"""
import csv, io, json, os, sys, zipfile
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from lifelines.statistics import logrank_test
from scipy.stats import chi2
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import ElasticNet

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_kit import Run, sha256, boot_ci

ZIP, CLIN = "/tmp/rf.zip", "/tmp/upenn.csv"
SEED = 20260807
MODS, ROIS = ["T1", "T1GD", "T2", "FLAIR"], ["ET", "ED", "NC"]
EXPECT = {"zip_bytes": 16121402, "n_both_markers": 291, "lr_alpha": 0.05}


def main():
    run = Run("R05head2head")
    cfg = {
        "report": "R05",
        "question": "imaging vs molecular markers, head to head in the same UPENN patients",
        "motivation": "R04 audit: R02 (imaging) and R04 (molecular) used different cohorts",
        "Q1_positive_control": "MGMT prognostic within this subset",
        "Q2": "molecular alone (IDH1 + MGMT), cross-validated C-index",
        "Q3": "imaging alone, same patients, cross-validated C-index",
        "Q4": "combined; LR test for imaging adding beyond molecular",
        "decision_rule": "imaging adds independent value only if LR p < 0.05",
        "cv": "5-fold, out-of-fold predictions, selection inside folds",
        "seed": SEED, "expectations": EXPECT,
        "relates_to": ["R02 imaging C=0.602 (n=574)", "R04 molecular within-GBM C=0.541-0.561"],
    }
    cfg_hash = run.start(cfg, [ZIP, CLIN])
    run.gate("G0_integrity", EXPECT["zip_bytes"], os.path.getsize(ZIP),
             os.path.getsize(ZIP) == EXPECT["zip_bytes"])

    # features, complete cases (same rule as R01/R02)
    z = zipfile.ZipFile(ZIP)
    tab, names, per = {}, [], {}
    for m in MODS:
        for roi in ROIS:
            with z.open(f"Radiomic_Features_CaPTk_automaticsegm_{m}_{roi}.csv") as fh:
                rows = list(csv.reader(io.TextIOWrapper(fh, "utf8")))
            cols = [f"{m}_{roi}::{c}" for c in rows[0][1:]]
            names += cols
            pres = set()
            for r in rows[1:]:
                if not r:
                    continue
                pres.add(r[0].strip())
                tab.setdefault(r[0].strip(), {}).update(
                    {c: (float(v) if v not in ("", "NA", "NaN") else np.nan)
                     for c, v in zip(cols, r[1:])})
            per[f"{m}_{roi}"] = pres
    complete = set.intersection(*per.values())

    def num(v):
        try:
            return float(v)
        except Exception:
            return None
    clin = {r["ID"].strip(): r for r in csv.DictReader(open(CLIN))}
    ids = sorted(i for i in complete if i in clin
                 and num(clin[i]["Survival_from_surgery_days_UPDATED"])
                 and num(clin[i]["Survival_from_surgery_days_UPDATED"]) > 0
                 and clin[i]["IDH1"].strip() in ("Wildtype", "Mutated")
                 and clin[i]["MGMT"].strip() in ("Methylated", "Unmethylated"))
    run.gate("G1_head2head_n", "<= 291 (subset with both markers)", len(ids),
             0 < len(ids) <= EXPECT["n_both_markers"],
             "features + survival + definite IDH1 + definite MGMT")
    run.log("cohort", n=len(ids))

    T = np.array([num(clin[i]["Survival_from_surgery_days_UPDATED"]) for i in ids])
    E = np.ones(len(ids), int)   # uncensored cohort, as established in R02
    idh = np.array([1 if clin[i]["IDH1"].strip() == "Mutated" else 0 for i in ids])
    mgmt = np.array([1 if clin[i]["MGMT"].strip() == "Methylated" else 0 for i in ids])
    age = np.array([num(clin[i]["Age_at_scan_years"]) or np.nan for i in ids])
    X = np.array([[tab[i].get(c, np.nan) for c in names] for i in ids], float)
    X = np.where(np.isnan(X), np.nanmedian(X, 0), X)
    run.gate("G2_uncensored", "all deceased", sorted({clin[i]["Survival_Status"].strip() for i in ids}),
             {clin[i]["Survival_Status"].strip() for i in ids} == {"Deceased"})

    # Q1 positive control
    lr = logrank_test(T[mgmt == 1], T[mgmt == 0], E[mgmt == 1], E[mgmt == 0])
    run.gate("G3_positive_control_MGMT", "MGMT prognostic p<0.05",
             {"p": f"{lr.p_value:.3e}", "median_meth": float(np.median(T[mgmt == 1])),
              "median_unmeth": float(np.median(T[mgmt == 0]))},
             lr.p_value < 0.05, "R02 found p=5.8e-8 on the larger subset")

    # Q3 imaging, cross-validated on THESE patients
    logt = np.log10(T)
    pipe = Pipeline([("var", VarianceThreshold(1e-8)), ("scale", StandardScaler()),
                     ("reg", ElasticNet(max_iter=20000, random_state=SEED))])
    grid = {"reg__alpha": [0.01, 0.1, 1.0], "reg__l1_ratio": [0.1, 0.5, 0.9]}
    kf = KFold(5, shuffle=True, random_state=SEED)
    rad = np.full(len(ids), np.nan)
    for tr, te in kf.split(X):
        gs = GridSearchCV(pipe, grid, scoring="neg_mean_squared_error",
                          cv=KFold(5, shuffle=True, random_state=SEED), n_jobs=-1)
        gs.fit(X[tr], logt[tr])
        rad[te] = gs.predict(X[te])
    c_img = float(concordance_index(T, rad, E))

    # Q2 molecular, cross-validated the same way
    def cv_cox(df, cols):
        pred = np.full(len(df), np.nan)
        for tr, te in KFold(5, shuffle=True, random_state=SEED).split(df):
            m = CoxPHFitter().fit(df.iloc[tr][["T", "E"] + cols], "T", "E")
            pred[te] = -m.predict_partial_hazard(df.iloc[te][cols]).values
        return float(concordance_index(df["T"], pred, df["E"])), pred

    def z(v):
        return (v - np.nanmean(v)) / np.nanstd(v)
    df = pd.DataFrame({"T": T, "E": E, "idh": idh, "mgmt": mgmt,
                       "rad": z(rad), "age": z(age)})
    c_mol, _ = cv_cox(df, ["idh", "mgmt"])
    c_comb, _ = cv_cox(df, ["idh", "mgmt", "rad"])
    c_all, _ = cv_cox(df, ["idh", "mgmt", "rad", "age"])

    # Q4 LR test: does imaging add beyond molecular?
    m1 = CoxPHFitter().fit(df[["T", "E", "idh", "mgmt"]], "T", "E")
    m2 = CoxPHFitter().fit(df[["T", "E", "idh", "mgmt", "rad"]], "T", "E")
    lr_stat = 2 * (m2.log_likelihood_ - m1.log_likelihood_)
    lr_p = float(chi2.sf(lr_stat, 1))
    # and the reverse: does molecular add beyond imaging?
    m3 = CoxPHFitter().fit(df[["T", "E", "rad"]], "T", "E")
    lr_rev = 2 * (m2.log_likelihood_ - m3.log_likelihood_)
    lr_rev_p = float(chi2.sf(lr_rev, 2))

    adds = lr_p < EXPECT["lr_alpha"]
    results = {
        "cohort": {"n": len(ids), "idh_mutant": int(idh.sum()),
                   "mgmt_methylated": int(mgmt.sum()), "all_deceased": True},
        "Q1_positive_control_MGMT": {"p": float(lr.p_value),
                                     "median_days_methylated": float(np.median(T[mgmt == 1])),
                                     "median_days_unmethylated": float(np.median(T[mgmt == 0]))},
        "Q2_molecular_alone_cindex": c_mol,
        "Q3_imaging_alone_cindex": c_img,
        "Q4_combined_cindex": c_comb,
        "Q4b_combined_plus_age_cindex": c_all,
        "LR_imaging_beyond_molecular": {"chi2": float(lr_stat), "df": 1, "p": lr_p},
        "LR_molecular_beyond_imaging": {"chi2": float(lr_rev), "df": 2, "p": lr_rev_p},
        "cox_terms": {k: {"HR": float(np.exp(m2.params_[k])),
                          "p": float(m2.summary.loc[k, "p"])} for k in ("idh", "mgmt", "rad")},
        "_decision": {"rule": f"imaging adds independent value if LR p < {EXPECT['lr_alpha']}",
                      "observed_p": lr_p,
                      "decision": "imaging ADDS independent prognostic value beyond IDH1+MGMT"
                                  if adds else
                                  "imaging does NOT add beyond IDH1+MGMT in this subset",
                      "prespecified_in_config_sha256": cfg_hash},
        "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir,
                        "n_features": int(X.shape[1])},
    }
    np.savetxt(os.path.join(run.dir, "head2head.csv"),
               np.column_stack([T, idh, mgmt, rad, age]), delimiter=",",
               header="survival_days,idh_mutated,mgmt_methylated,radiomic_score,age", comments="")
    run.write("results.json", results)
    run.log("decision", **results["_decision"])
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
