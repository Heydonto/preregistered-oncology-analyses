#!/usr/bin/env python3
"""R08 — the lung chromatin arm with methylation. TCGA-LUAD + LUSC, 450k arrays.

R06 showed KEAP1 and SMARCA4 mutations are prognostic in NSCLC using mutation data alone.
This report adds the methylation layer the response letter proposed, on a cohort large enough
to carry it: 843 primary tumours and 74 matched-type normals across 833 patients.

Four pre-registered positive controls, in descending order of expected effect. This cohort was
chosen partly because it can power the AHRR control that R07 could not (93 lifelong
non-smokers against 256 current smokers, versus R07's 24 smokers):

  PC1 tumour vs normal   thousands of differentially methylated probes; the largest effect
                         available and a direct test that the matrix is correctly assembled
  PC2 sex                near-perfect separation from X-linked probes; also verifies the
                         sample-to-clinical join, the failure mode R07 had to rule out
  PC3 AHRR cg05575921    lower methylation in smokers; the control R07 could not power
  PC4 LUAD vs LUSC       histology is known to be separable from methylation

Gate: PC1 and PC2 must both hold. They test assembly and join correctness, so a failure means
the data are wrong rather than the biology surprising. PC3 and PC4 are reported either way.

Primary questions, fixed before any outcome is read:
  Q1 do KEAP1 / SMARCA4 / KMT2D mutations associate with differential methylation?
  Q2 is there a methylation signature prognostic for overall survival, nested-CV honest?
  Q3 does methylation add prognostic value beyond mutation status, stage and age?

Q3 is the question that matters for the arm: R06 established the mutation signal, so the
methylation layer has to earn its place on top of it, not merely correlate with survival.
"""
import csv, json, os, subprocess, sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from lifelines.statistics import logrank_test
from scipy.stats import chi2
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.linear_model import ElasticNet

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_kit import Run, bh, boot_ci

GCS = "gs://heydonto-quantara-lungcdx/data-request-2026-08/lung/tcga_methylation_450k"
LOCAL = "/Users/rezanehzati/quantara-staging/r08"
CLIN = LOCAL + "/clinical_merged.tsv"
GENES = {9817: "KEAP1", 6597: "SMARCA4", 8085: "KMT2D", 6794: "STK11", 7157: "TP53"}
CHROM = ["KEAP1", "SMARCA4", "KMT2D"]
SILENT = {"Silent", "Intron", "3'UTR", "5'UTR", "5'Flank", "IGR", "RNA"}
AHRR = "cg05575921"
SEED = 20260807
EXPECT = {"files": 919, "probes": 486427, "tumours": 843, "normals": 74,
          "max_probe_nan_frac": 0.05}


def sh(*a):
    return subprocess.run(a, capture_output=True, text=True).stdout


def load_matrix(run):
    for f in ("beta_450k.npy", "probes.txt", "samples.tsv", "ACQUISITION.txt",
              "MANIFEST_450k.tsv"):
        if not os.path.exists(f"{LOCAL}/{f}"):
            run.log("fetch", file=f)
            subprocess.run(["gcloud", "storage", "cp", f"{GCS}/{f}", f"{LOCAL}/{f}"], check=True)
    M = np.load(f"{LOCAL}/beta_450k.npy", mmap_mode="r")
    probes = open(f"{LOCAL}/probes.txt").read().split()
    samples = list(csv.DictReader(open(f"{LOCAL}/samples.tsv"), delimiter="\t"))
    return M, probes, samples


def mutations():
    per = {g: set() for g in GENES.values()}
    n = 0
    for study in ("luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"):
        body = json.dumps({"sampleListId": f"{study}_all", "entrezGeneIds": list(GENES)})
        out = subprocess.run(
            ["curl", "-s", "--retry", "3", "-X", "POST",
             f"https://www.cbioportal.org/api/molecular-profiles/{study}_mutations/"
             f"mutations/fetch?projection=DETAILED",
             "-H", "Content-Type: application/json", "-d", body],
            capture_output=True, timeout=300).stdout
        d = json.loads(out)
        n += len(d)
        for x in d:
            g = GENES.get(x.get("entrezGeneId"))
            if g and x.get("mutationType") not in SILENT:
                pid = x.get("patientId")
                if pid:
                    per[g].add(pid)
    return per, n


def main():
    run = Run("R08tcgameth")
    cfg = {
        "report": "R08",
        "question": "does DNA methylation add to the NSCLC chromatin/prognostic picture "
                    "established in R06?",
        "cohort": "TCGA-LUAD + TCGA-LUSC, Illumina 450k, SeSAMe level-3 betas from GDC",
        "clinical_source": "survival/sex/age/stage from the cBioPortal PanCancer Atlas (the curated TCGA-CDR endpoint, which is the published recommendation for TCGA survival analysis); smoking status from the GDC cases endpoint, which carries it and TCGA-CDR does not. GDC demographic.gender is empty for this project - the field is now sex_at_birth - so the curated SEX field is used instead",
        "acquisition": "919 open-access files, MD5-verified, assembled on a transient VM",
        "positive_controls": {
            "PC1": "tumour vs normal differential methylation (largest available effect)",
            "PC2": "sex from X-linked probes (also verifies the sample-clinical join)",
            "PC3": f"{AHRR} lower in smokers (the control R07 could not power)",
            "PC4": "LUAD vs LUSC separable"},
        "gate_rule": "PC1 and PC2 must both hold - they test assembly and join correctness",
        "Q1": "differential methylation by KEAP1/SMARCA4/KMT2D mutation status, BH-corrected",
        "Q2": "prognostic methylation signature for OS, nested CV",
        "Q3": "does methylation add beyond mutation status + stage + age (LR test)",
        "decision_rule_Q3": "methylation adds only if LR p < 0.05 AND the seed-stability check "
                            "in the audit clears 0.05 in a majority of 8 seeds",
        "probe_filter": f"drop probes with >{EXPECT['max_probe_nan_frac']:.0%} NA across "
                        "tumours; SeSAMe masks unreliable probes as NA",
        "platform_restriction": "450k ONLY. The 311 HM27 and 53 EPICv2 files in this project "
                                "were excluded before download: mixing array generations "
                                "would confound every comparison",
        "replicate_rule": "4 samples have 2 files; the lowest file_id is kept, deterministically",
        "relates_to": ["R06 KEAP1 HR 2.07, SMARCA4 HR 1.73 (mutation only)",
                       "R07 AHRR control underpowered at n=79"],
        "seed": SEED, "expectations": EXPECT,
    }
    cfg_hash = run.start(cfg, [CLIN])

    M, probes, samples = load_matrix(run)
    acq = dict(l.split("=", 1) for l in
               open(f"{LOCAL}/ACQUISITION.txt").read().strip().split("\n"))
    run.gate("G0_files", EXPECT["files"], int(acq["files"]), int(acq["files"]) == EXPECT["files"],
             "all MD5-verified on the VM; fewer is a failure, not a smaller matrix")
    run.gate("G1_shape", (EXPECT["probes"], EXPECT["files"]), tuple(M.shape),
             tuple(M.shape) == (EXPECT["probes"], EXPECT["files"]))
    run.log("acquisition", **acq)

    # ---- deterministic replicate resolution + column bookkeeping ----
    by_sample = {}
    for s in samples:
        k = s["sample"]
        if k not in by_sample or s["file_id"] < by_sample[k]["file_id"]:
            by_sample[k] = s
    keep_cols = sorted(int(s["col"]) for s in by_sample.values())
    dropped = len(samples) - len(keep_cols)
    run.gate("G2_replicates", 4, dropped, dropped == 4,
             "4 samples carry 2 files; lowest file_id kept")
    smeta = [samples[c] for c in keep_cols]
    is_tum = np.array([s["sample_type"] != "Solid Tissue Normal" for s in smeta])
    run.gate("G3_tissue_split", {"tumour": EXPECT["tumours"], "normal": EXPECT["normals"]},
             {"tumour": int(is_tum.sum()), "normal": int((~is_tum).sum())},
             int((~is_tum).sum()) == EXPECT["normals"])

    B = np.asarray(M[:, keep_cols], np.float32)
    del M

    # ---- probe filter, computed on tumours only so normals cannot drive inclusion ----
    nan_frac = np.isnan(B[:, is_tum]).mean(1)
    good = nan_frac <= EXPECT["max_probe_nan_frac"]
    run.log("probe_filter", kept=int(good.sum()), dropped=int((~good).sum()),
            frac_kept=round(float(good.mean()), 4))
    run.gate("G4_probes_kept", "> 300,000 probes survive the NA filter", int(good.sum()),
             int(good.sum()) > 300_000)
    B = B[good]
    pnames = [p for p, k in zip(probes, good) if k]
    pidx = {p: i for i, p in enumerate(pnames)}

    # ---- PC1 tumour vs normal ----
    rng = np.random.default_rng(SEED)
    sub = rng.choice(len(pnames), 20000, replace=False)
    ps = []
    for i in sub:
        t, n = B[i, is_tum], B[i, ~is_tum]
        t, n = t[~np.isnan(t)], n[~np.isnan(n)]
        ps.append(mannwhitneyu(t, n).pvalue if len(t) > 20 and len(n) > 20 else 1.0)
    q = bh(ps)
    n_dm = int((q < 0.05).sum())
    run.gate("PC1_tumour_vs_normal", "> 2,000 of 20,000 probes differential at q<0.05", n_dm,
             n_dm > 2000, f"{int(is_tum.sum())} tumours vs {int((~is_tum).sum())} normals")

    # ---- clinical join ----
    clin = {r["case"]: r for r in csv.DictReader(open(CLIN), delimiter="\t")}
    cases = [s["case"] for s in smeta]
    matched = sum(1 for c in cases if c in clin)
    run.gate("G5_clinical_join", "> 95% of samples match a clinical record",
             f"{matched}/{len(cases)}", matched / len(cases) > 0.95)

    sex = np.array([clin.get(c, {}).get("sex", "") for c in cases], object)
    # ---- PC2 sex ----
    ok = np.isin(sex, ["Male", "Female"]) & is_tum
    male = (sex == "Male")[ok]
    var = np.nanvar(B[:, ok], axis=1)
    cand = np.argsort(var)[-20000:]
    best, bestp = 0.5, None
    for i in cand:
        v = B[i, ok]
        m = ~np.isnan(v)
        if m.sum() < ok.sum() * 0.9 or len(set(male[m])) < 2:
            continue
        a = roc_auc_score(male[m].astype(int), v[m])
        a = max(a, 1 - a)
        if a > best:
            best, bestp = a, pnames[i]
    run.gate("PC2_sex", "a probe separates sexes at AUROC>0.95",
             {"auroc": round(best, 4), "probe": bestp, "n": int(ok.sum())}, best > 0.95,
             "also verifies the sample-to-clinical join")

    # ---- PC3 AHRR vs smoking (the control R07 could not power) ----
    smoke = np.array([clin.get(c, {}).get("smoke_status", "") for c in cases], object)
    never = np.array(["Non-Smoker" in s for s in smoke]) & is_tum
    curr = np.array([s == "Current Smoker" for s in smoke]) & is_tum
    pc3 = {"probe_present": AHRR in pidx}
    if AHRR in pidx:
        a, b = B[pidx[AHRR], curr], B[pidx[AHRR], never]
        a, b = a[~np.isnan(a)], b[~np.isnan(b)]
        u = mannwhitneyu(a, b)
        pc3.update({"n_current": int(len(a)), "n_never": int(len(b)),
                    "median_current": float(np.median(a)), "median_never": float(np.median(b)),
                    "lower_in_smokers": bool(np.median(a) < np.median(b)),
                    "p": float(u.pvalue),
                    "holds": bool(u.pvalue < 0.05 and np.median(a) < np.median(b))})
    run.log("PC3_AHRR_smoking", **pc3)

    # ---- PC4 LUAD vs LUSC ----
    proj = np.array([s["project"] for s in smeta], object)
    lu = is_tum & np.isin(proj, ["TCGA-LUAD", "TCGA-LUSC"])
    ylu = (proj[lu] == "TCGA-LUSC").astype(int)
    bestlu = 0.5
    for i in cand[-5000:]:
        v = B[i, lu]
        m = ~np.isnan(v)
        if m.sum() < lu.sum() * 0.9:
            continue
        a = roc_auc_score(ylu[m], v[m])
        bestlu = max(bestlu, max(a, 1 - a))
    run.log("PC4_luad_vs_lusc", best_probe_auroc=round(float(bestlu), 4),
            n_luad=int((ylu == 0).sum()), n_lusc=int(ylu.sum()),
            holds=bool(bestlu > 0.80))

    # ---- mutation join ----
    per, nrec = mutations()
    run.gate("G6_mutation_api", "> 500 records", nrec, nrec > 500)
    mut = {g: np.array([c in per[g] for c in cases]) for g in GENES.values()}
    run.log("mutation_counts", **{g: int((mut[g] & is_tum).sum()) for g in GENES.values()})

    # ---- Q1 differential methylation by chromatin-gene mutation ----
    q1 = {}
    for g in CHROM:
        m1 = mut[g] & is_tum
        m0 = (~mut[g]) & is_tum
        if m1.sum() < 20:
            q1[g] = {"skipped": "fewer than 20 mutated tumours", "n": int(m1.sum())}
            continue
        ps = []
        for i in sub:
            a, b = B[i, m1], B[i, m0]
            a, b = a[~np.isnan(a)], b[~np.isnan(b)]
            ps.append(mannwhitneyu(a, b).pvalue if len(a) > 10 and len(b) > 10 else 1.0)
        qq = bh(ps)
        q1[g] = {"n_mutated": int(m1.sum()), "n_wt": int(m0.sum()),
                 "probes_tested": len(sub), "n_sig_q05": int((qq < 0.05).sum()),
                 "frac_sig": float((qq < 0.05).mean()),
                 "min_q": float(qq.min())}
        run.log("Q1", gene=g, **{k: v for k, v in q1[g].items() if k != "probes_tested"})

    # ---- survival frame ----
    rows = []
    for j, s in enumerate(smeta):
        if not is_tum[j]:
            continue
        c = clin.get(s["case"])
        if not c:
            continue
        dead = c["os_status"].startswith("1")
        try:
            t = float(c["os_months"])
        except (TypeError, ValueError):
            continue
        if t <= 0:
            continue
        age = float(c["age"]) if c["age"] else np.nan
        # cBioPortal writes these uppercase ("STAGE IIIA"); an earlier case-sensitive parse
        # silently made this covariate constant, which surfaced as a singular Cox matrix.
        # G8 below now gates on the parsed fraction so it cannot pass unnoticed.
        st = (c["stage"] or "").upper().replace("STAGE", "").strip()
        late = 1 if st.startswith(("II", "III", "IV")) else 0
        rows.append({"col": j, "T": t, "E": int(dead), "age": age, "late_stage": late,
                     "lusc": int(s["project"] == "TCGA-LUSC"),
                     **{g: int(mut[g][j]) for g in GENES.values()}})
    df = pd.DataFrame(rows)
    run.gate("G7_survival_n", "> 700 tumours with usable survival", len(df), len(df) > 700,
             "OS months > 0 and status known (TCGA-CDR)")
    run.log("survival_cohort", n=len(df), events=int(df["E"].sum()))
    late_frac = float(df["late_stage"].mean())
    run.gate("G8_stage_parsed", "late-stage fraction in 0.20-0.60", round(late_frac, 4),
             0.20 <= late_frac <= 0.60,
             "guards against a stage-string parse failure making the covariate constant")
    for cvname in ("KEAP1", "SMARCA4", "late_stage", "lusc"):
        run.gate(f"G9_variance_{cvname}", "covariate is not constant",
                 int(df[cvname].nunique()), df[cvname].nunique() > 1)

    cols = df["col"].to_numpy()
    Xm = B[:, cols].T.astype(np.float64)
    Xm = np.where(np.isnan(Xm), np.nanmedian(Xm, 0), Xm)
    T, E = df["T"].to_numpy(), df["E"].to_numpy()

    # ---- Q2 prognostic methylation signature, nested CV ----
    logt = np.log10(T)
    pipe = Pipeline([("sel", SelectKBest(f_regression, k=5000)), ("sc", StandardScaler()),
                     ("reg", ElasticNet(max_iter=5000, random_state=SEED))])
    grid = {"reg__alpha": [0.01, 0.1, 1.0], "reg__l1_ratio": [0.5]}
    sig = np.full(len(df), np.nan)
    for tr, te in KFold(5, shuffle=True, random_state=SEED).split(Xm):
        gs = GridSearchCV(pipe, grid, scoring="neg_mean_squared_error",
                          cv=KFold(3, shuffle=True, random_state=SEED), n_jobs=-1)
        gs.fit(Xm[tr], logt[tr])
        sig[te] = gs.predict(Xm[te])
        run.log("cv_fold_done", n_test=len(te))
    c_meth = float(concordance_index(T, sig, E))
    lo, hi = boot_ci(lambda ix: concordance_index(T[ix], sig[ix], E[ix]), len(T), SEED, B=2000)
    run.log("Q2", cindex=round(c_meth, 4), ci95=[round(lo, 4), round(hi, 4)])

    # ---- Q3 does methylation add beyond mutation + stage + age? ----
    def zs(v):
        return (v - np.nanmean(v)) / np.nanstd(v)
    d2 = pd.DataFrame({"T": T, "E": E, "age": zs(df["age"].to_numpy()),
                       "late_stage": df["late_stage"].to_numpy(),
                       "lusc": df["lusc"].to_numpy(),
                       "KEAP1": df["KEAP1"].to_numpy(), "SMARCA4": df["SMARCA4"].to_numpy(),
                       "meth": zs(sig)}).dropna()
    base = ["age", "late_stage", "lusc", "KEAP1", "SMARCA4"]
    m1 = CoxPHFitter().fit(d2[["T", "E"] + base], "T", "E")
    m2 = CoxPHFitter().fit(d2[["T", "E"] + base + ["meth"]], "T", "E")
    lr = 2 * (m2.log_likelihood_ - m1.log_likelihood_)
    lr_p = float(chi2.sf(lr, 1))

    def cv_cox(frame, cc):
        """Rare binary covariates (SMARCA4 has 52 carriers) can be constant inside a fold,
        which makes the Cox design singular. Such columns are dropped for that fold only and
        the drop is logged, rather than the fold being silently skipped."""
        pred = np.full(len(frame), np.nan)
        for k, (tr, te) in enumerate(KFold(5, shuffle=True, random_state=SEED).split(frame)):
            use = [c for c in cc if frame.iloc[tr][c].nunique() > 1]
            if len(use) < len(cc):
                run.log("cv_cox_dropped_constant", fold=k,
                        dropped=[c for c in cc if c not in use])
            mm = CoxPHFitter(penalizer=0.01).fit(frame.iloc[tr][["T", "E"] + use], "T", "E")
            pred[te] = -mm.predict_partial_hazard(frame.iloc[te][use]).values
        return float(concordance_index(frame["T"], pred, frame["E"]))

    c_base = cv_cox(d2, base)
    c_full = cv_cox(d2, base + ["meth"])
    adds = lr_p < 0.05
    run.log("Q3", lr_chi2=round(lr, 3), lr_p=lr_p, c_base=round(c_base, 4),
            c_full=round(c_full, 4))

    results = {
        "acquisition": acq,
        "cohort": {"tumours": int(is_tum.sum()), "normals": int((~is_tum).sum()),
                   "probes_kept": int(good.sum()), "survival_n": len(df),
                   "events": int(df["E"].sum())},
        "PC1_tumour_vs_normal": {"probes_tested": len(sub), "n_sig_q05": n_dm,
                                 "frac_sig": n_dm / len(sub)},
        "PC2_sex": {"auroc": best, "probe": bestp},
        "PC3_AHRR_smoking": pc3,
        "PC4_luad_vs_lusc": {"best_probe_auroc": float(bestlu)},
        "Q1_differential_methylation": q1,
        "Q2_methylation_signature": {"cindex": c_meth, "ci95": [lo, hi]},
        "Q3_incremental": {
            "covariates_base": base, "lr_chi2": float(lr), "lr_p": lr_p,
            "cindex_base": c_base, "cindex_with_methylation": c_full,
            "cox_terms": {k: {"HR": float(np.exp(m2.params_[k])),
                              "p": float(m2.summary.loc[k, "p"])}
                          for k in base + ["meth"]}},
        "_decision": {
            "rule": "methylation adds only if LR p<0.05 and the audit's seed check agrees",
            "observed_lr_p": lr_p,
            "provisional": "methylation adds beyond mutation+stage+age" if adds
                           else "methylation does NOT add beyond mutation+stage+age",
            "note": "provisional pending the audit seed-stability check, per config",
            "prespecified_in_config_sha256": cfg_hash},
        "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir},
    }
    np.savetxt(os.path.join(run.dir, "oof_methylation_signature.csv"),
               np.column_stack([T, E, sig]), delimiter=",",
               header="os_months,event,methylation_score", comments="")
    run.write("results.json", results)
    run.write("survival_frame.csv", d2.to_csv(index=False))
    run.log("decision", **results["_decision"])
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
