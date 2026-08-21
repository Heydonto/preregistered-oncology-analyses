#!/usr/bin/env python3
"""Arm B, methylation track — 450k -> EPIC cross-array transfer of the R10 KEAP1-phenotype
model (plus a TCGA-LUAD OS Cox head) onto the sealed anti-PD-1 EPIC cohort (GSE115246).

Three stages, run in order. The stage boundary between 1 and 3 IS the seal:

  stage1  TCGA-LUAD only. Writes config.yaml (the whole plan, including everything stage 3
          will do) and hashes it. Reproduces R10's KEAP1 classifier as a positive-control
          gate, builds the probe-overlap map from the two official Illumina manifests,
          fits the deployable models (KEAP1 full-array, KEAP1 overlap-restricted, sex
          control, elastic-net OS Cox). NO EPIC file is opened.
  stage2  appends the unsealing event (with stage1's config hash) to
          gs://.../nsclc-rwpr-study/_provenance/unsealing_events.jsonl. NO EPIC file is
          opened.
  stage3  reads the anti-PD-1 EPIC cohort and evaluates zero-shot.

The pre-declared expectation for the response endpoint is an HONEST INCONCLUSIVE at n=81:
81 patients cannot resolve a modest AUROC. The deliverable is the measured
minimum-detectable AUROC alongside whatever the point estimate turns out to be, plus a
transfer positive control (sex predicted zero-shot across arrays) that proves the sample
mapping and the cross-array pipeline are sound even when the response endpoint is not
resolvable.
"""
import csv, gzip, io, json, os, subprocess, sys, time
import numpy as np
from scipy.stats import mannwhitneyu, norm
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, GridSearchCV, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from audit_kit import Run, bh, boot_ci, sha256   # noqa: E402

# ---------------------------------------------------------------- paths / constants
R08 = "/Users/rezanehzati/quantara-staging/r08"                       # TCGA 450k matrix (R10 inputs)
MAN450 = "/tmp/manifest450k.csv.gz"                                   # Illumina 450k v1.1 (GPL13534)
MANEPIC = "/tmp/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz"       # Illumina EPIC v1.0 (GPL21145)
EPIC_STAGE = "/Users/rezanehzati/quantara-staging/armB-meth"
SERIES_MATRIX = f"{EPIC_STAGE}/GSE115246_series_matrix.txt.gz"
SIGNALS = f"{EPIC_STAGE}/GSE115246_methylated_unmethylated_signal_intensities.txt.gz"
WORK = os.path.join(HERE, "work")
GCS_PROV = "gs://heydonto-quantara-lungcdx/nsclc-rwpr-study/_provenance/unsealing_events.jsonl"

GENES = {9817: "KEAP1", 6597: "SMARCA4", 8085: "KMT2D", 6794: "STK11", 7157: "TP53"}
SILENT = {"Silent", "Intron", "3'UTR", "5'UTR", "5'Flank", "IGR", "RNA"}
SEED = 20260810
R10_SEED = 20260807          # R10's seed — kept for the reproduction leg
KBEST = 5000                 # R10's SelectKBest k
COX_K = 500                  # univariate-Cox screen width for the OS head
N_PERM = 1000                # permutation shuffles (task floor: 100)
EXPECT = {"probes": 392489, "luad_tumours": 471, "keap1_mutant_luad": 84,
          "R10_keap1_auroc": 0.9098375784422296, "R10_sex_auroc": 1.0,
          "gate_repro_auroc": 0.88}


def cfg_dict():
    """The pre-registered plan. Hashed before any EPIC byte is read."""
    return {
        "study": "NSCLC RWPR — Arm B, methylation track",
        "arm": "B-meth",
        "question": "does the R10 KEAP1-methylation phenotype (and a TCGA-LUAD OS Cox head) "
                    "transfer 450k->EPIC and carry any signal for anti-PD-1 response in the "
                    "sealed 81-patient EPIC cohort?",
        "pre_declared_expectation": "HONEST INCONCLUSIVE on the response endpoint at n=81. "
                                    "A calibrated null with a measured minimum-detectable "
                                    "AUROC is the success condition; a manufactured positive "
                                    "is the failure condition.",
        "training_cohort": "TCGA-LUAD primary tumours, Illumina 450k, n=471, "
                           "identical inputs/filters to R10 (same beta matrix, same "
                           "392,489-probe filter, same cbioportal mutation call)",
        "holdout_cohort": "GSE115246 — 81 advanced NSCLC, anti-PD-1, Illumina EPIC "
                          "(read-only from gs://heydonto-quantara-lungcdx/data-request-2026-08/"
                          "labels/lung/GSE115246/)",
        "models": {
            "M1_keap1_full": "R10 pipeline verbatim: SelectKBest(f_classif,k=5000) -> "
                             "StandardScaler -> LogisticRegression(elasticnet, l1_ratio=0.5, "
                             "saga, max_iter=3000, tol=1e-3), C from inner 3-fold "
                             "GridSearchCV over [0.01,0.1,1.0]. Feature space = all 392,489 "
                             "retained 450k probes. Nested 5-fold CV AUROC is the R10 "
                             "reproduction number; a full-cohort refit is the deployable model.",
            "M1_keap1_overlap": "identical pipeline, feature space restricted a priori to "
                                "450k-probes-also-on-EPIC (official Illumina manifest "
                                "intersection). This is the transferable variant; both are "
                                "reported.",
            "M2_os_cox": "elastic-net Cox on the overlap feature space. Screen: vectorised "
                         "Cox score test at beta=0 (top 500 probes), inside each CV fold. "
                         "Fit: lifelines CoxPHFitter(penalizer=0.1, l1_ratio=0.5) on "
                         "z-scored probes. Honest C-index by 5-fold CV; full-cohort refit "
                         "is the deployable score.",
            "PC_sex": "same classifier pipeline, target = sex, overlap feature space. "
                      "Used as the CROSS-ARRAY transfer positive control: applied zero-shot "
                      "to the EPIC cohort it must recover sex at AUROC>=0.95, which proves "
                      "sample mapping and that 450k->EPIC transfer works at all."},
        "probe_mapping": {
            "source": "official Illumina manifests: HumanMethylation450 v1.1 (GEO GPL13534 "
                      "supplementary, sha256 pinned in inputs.json) and MethylationEPIC "
                      "15073387 v1.0 (GPL21145, shipped inside GSE115246_RAW.tar)",
            "rule": "a probe is transferable iff its IlmnID appears in both manifests AND is "
                    "present in the EPIC deposit's own probe index",
            "report": "manifest overlap size; how many of each model's 5,000 selected and "
                      "non-zero-coefficient probes survive; empirical coverage in the deposit"},
        "harmonisation": {
            "T1_naive": "apply the TCGA-fitted StandardScaler unchanged to EPIC betas "
                        "(measures raw cross-array shift; expected to be biased)",
            "T2_cohort_z": "PRIMARY. Replace the scaler's per-probe mean/sd with the EPIC "
                           "cohort's own per-probe mean/sd, then apply the frozen model "
                           "coefficients. Label-free, standard cross-cohort z-harmonisation; "
                           "removes array-level offset and scale without touching outcomes.",
            "missing_probe_rule": "a model probe absent from the EPIC deposit is imputed at "
                                  "the TCGA training median (T1) / at 0 after z-scoring (T2), "
                                  "and counted in the coverage report"},
        "endpoints": {
            "E1_primary": "KEAP1-phenotype score -> anti-PD-1 response (responder vs "
                          "non-responder as deposited). AUROC + 2000-replicate bootstrap 95% CI.",
            "E2": "OS/PFS association of the Cox score, IF the deposit carries survival: "
                  "Cox HR per SD with CI, Harrell C, log-rank on the median split.",
            "E3": f"measured power: {N_PERM} label permutations -> null AUROC distribution, "
                  "one-sided p = (#null>=obs + 1)/(N+1); minimum detectable AUROC at 80% "
                  "power via the R08 convention mde = 0.5 + (z.975 + z.80) * se with se "
                  "from the observed bootstrap CI width; plus the non-parametric 95th "
                  "percentile of the permutation null as a second floor."},
        "gates": {
            "G1_R10_reproduction": "nested-CV KEAP1 AUROC on TCGA >= 0.88 (R10 reported 0.9098)",
            "G2_R10_sex_control": "sex AUROC on TCGA > 0.95 through the same pipeline",
            "G3_preregistered_before_unseal": "config sha256 appended to unsealing_events.jsonl "
                                              "before any EPIC byte is read",
            "G4_probe_overlap_reported": "reporting gate — coverage numbers must be present",
            "G5_epic_sex_transfer": "sex recovered zero-shot on EPIC at AUROC >= 0.95 IF sex "
                                    "metadata ships; if it does not ship, gate is recorded "
                                    "NOT_TESTABLE with the reason, never silently dropped",
            "G6_permutation_null_centred": "permutation null mean AUROC within 0.5 +/- 0.05",
            "G7_no_silent_drops": "every deposited sample is accounted for: used, or excluded "
                                  "with a per-sample reason"},
        "sample_qc": "exclude a sample if >5% of the deployed model's probes are missing, or "
                     "if its mean beta over model probes is >5 SD from the cohort mean. Every "
                     "exclusion logged per sample with the reason.",
        "seed": SEED, "r10_seed": R10_SEED, "n_permutations": N_PERM,
        "expectations": EXPECT,
        "hard_constraints": ["data-request-2026-08/ is READ-ONLY",
                             "nsclc-rwpr-study/sealed/ never opened",
                             "CPTAC labels out of scope"],
        "relates_to": ["R10 KEAP1 phenotype AUROC 0.910 (CI 0.869-0.946), sex control 1.000",
                       "R08 methylation adds nothing prognostically in TCGA lung (C=0.486)",
                       "R07 methylation->TKI response inconclusive, MDE 0.685 at n=69"],
    }


# ---------------------------------------------------------------- TCGA loading (R10 verbatim)
def load_tcga(run):
    M = np.load(f"{R08}/beta_450k.npy", mmap_mode="r")
    probes = open(f"{R08}/probes.txt").read().split()
    samples = list(csv.DictReader(open(f"{R08}/samples.tsv"), delimiter="\t"))
    by = {}
    for s in samples:
        if s["sample"] not in by or s["file_id"] < by[s["sample"]]["file_id"]:
            by[s["sample"]] = s
    cols = sorted(int(s["col"]) for s in by.values())
    smeta = [samples[c] for c in cols]
    is_tum = np.array([s["sample_type"] != "Solid Tissue Normal" for s in smeta])
    B = np.asarray(M[:, cols], np.float32)
    del M
    good = np.isnan(B[:, is_tum]).mean(1) <= 0.05
    B = B[good]
    pnames = [p for p, k in zip(probes, good) if k]
    run.gate("G0_probes", EXPECT["probes"], len(pnames), len(pnames) == EXPECT["probes"],
             "identical filter to R08/R10")
    clin = {r["case"]: r for r in csv.DictReader(open(f"{R08}/clinical_merged.tsv"),
                                                delimiter="\t")}
    cases = [s["case"] for s in smeta]
    proj = np.array([s["project"] for s in smeta], object)
    luad = is_tum & (proj == "TCGA-LUAD")
    run.gate("G1_luad", EXPECT["luad_tumours"], int(luad.sum()),
             int(luad.sum()) == EXPECT["luad_tumours"])

    per = {g: set() for g in GENES.values()}
    for study in ("luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"):
        d = json.loads(subprocess.run(
            ["curl", "-s", "--retry", "3", "-X", "POST",
             f"https://www.cbioportal.org/api/molecular-profiles/{study}_mutations/"
             f"mutations/fetch?projection=DETAILED", "-H", "Content-Type: application/json",
             "-d", json.dumps({"sampleListId": f"{study}_all",
                               "entrezGeneIds": list(GENES)})],
            capture_output=True, timeout=300).stdout)
        for x in d:
            g = GENES.get(x.get("entrezGeneId"))
            if g and x.get("mutationType") not in SILENT and x.get("patientId"):
                per[g].add(x["patientId"])
    keap1 = np.array([c in per["KEAP1"] for c in cases])
    run.gate("G2_keap1_luad", EXPECT["keap1_mutant_luad"], int((keap1 & luad).sum()),
             int((keap1 & luad).sum()) == EXPECT["keap1_mutant_luad"])
    sex = np.array([clin.get(c, {}).get("sex", "") for c in cases], object)
    os_m = np.array([float(clin.get(c, {}).get("os_months") or "nan") for c in cases])
    os_e = np.array([str(clin.get(c, {}).get("os_status", "")).startswith("1")
                     for c in cases])
    return B, pnames, cases, luad, keap1, sex, os_m, os_e


def manifest_probes(path):
    s = set()
    with gzip.open(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("IlmnID"):
                break
        for line in fh:
            n = line.split(",", 2)[0]
            if n.startswith(("cg", "ch.", "rs")):
                s.add(n)
    return s


# ---------------------------------------------------------------- classifier (R10 pipeline)
def make_pipe(seed):
    return Pipeline([("sel", SelectKBest(f_classif, k=KBEST)), ("sc", StandardScaler()),
                     ("clf", LogisticRegression(penalty="elasticnet", solver="saga",
                                                l1_ratio=0.5, max_iter=3000, tol=1e-3,
                                                random_state=seed))])


def nested_cv(X, y, seed):
    oof = np.full(len(y), np.nan)
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(X, y):
        gs = GridSearchCV(make_pipe(seed), {"clf__C": [0.01, 0.1, 1.0]}, scoring="roc_auc",
                          cv=StratifiedKFold(3, shuffle=True, random_state=seed), n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.predict_proba(X[te])[:, 1]
    auc = float(roc_auc_score(y, oof))
    lo, hi = boot_ci(lambda ix: roc_auc_score(y[ix], oof[ix]) if len(set(y[ix])) == 2
                     else np.nan, len(y), seed, B=2000)
    return {"auroc": auc, "ci95": [lo, hi], "n": int(len(y)), "positives": int(y.sum())}, oof


def fit_full(X, y, seed):
    gs = GridSearchCV(make_pipe(seed), {"clf__C": [0.01, 0.1, 1.0]}, scoring="roc_auc",
                      cv=StratifiedKFold(3, shuffle=True, random_state=seed), n_jobs=-1)
    gs.fit(X, y)
    best = gs.best_estimator_
    sel = best.named_steps["sel"].get_support(indices=True)
    sc = best.named_steps["sc"]
    clf = best.named_steps["clf"]
    return {"selected_idx": sel.tolist(), "mean": sc.mean_.tolist(), "scale": sc.scale_.tolist(),
            "coef": clf.coef_[0].tolist(), "intercept": float(clf.intercept_[0]),
            "C": float(gs.best_params_["clf__C"]),
            "n_selected": int(len(sel)), "n_nonzero": int((clf.coef_[0] != 0).sum())}


def impute_median(X):
    """columns = probes; NaN -> column median. Returns (X, medians)."""
    med = np.nanmedian(X, 0)
    med = np.where(np.isnan(med), 0.5, med)
    return np.where(np.isnan(X), med, X), med


# ---------------------------------------------------------------- Cox
def cox_score_z(X, t, e, chunk=20000):
    """Vectorised Cox score-test z at beta=0 for every column of X independently.
    X: (n, p) no NaN. t: times. e: event bool. Breslow handling of the risk set.
    Chunked over probes so peak memory stays ~ (n x chunk)."""
    o = np.argsort(t, kind="stable")
    es = np.asarray(e, bool)[o]
    n = X.shape[0]
    m = np.arange(n, 0, -1).astype(float)[:, None]
    out = np.empty(X.shape[1])
    for a in range(0, X.shape[1], chunk):
        Xs = X[o, a:a + chunk]
        c1 = np.cumsum(Xs[::-1], 0)[::-1]
        c2 = np.cumsum((Xs ** 2)[::-1], 0)[::-1]
        mean = c1 / m
        var = c2 / m - mean ** 2
        U = (Xs[es] - mean[es]).sum(0)
        I = var[es].sum(0)
        out[a:a + chunk] = U / np.sqrt(np.maximum(I, 1e-12))
    return out


def fit_cox(Xtr, ttr, etr, k=COX_K):
    """screen by |score z| then penalised Cox on z-scored probes. Returns model dict."""
    import pandas as pd
    from lifelines import CoxPHFitter
    z = np.abs(cox_score_z(Xtr, ttr, etr))
    idx = np.argsort(-z)[:k]
    Z = Xtr[:, idx]
    mu, sd = Z.mean(0), Z.std(0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    df = pd.DataFrame((Z - mu) / sd, columns=[f"p{i}" for i in range(len(idx))])
    df["T"] = np.maximum(ttr, 1e-3)
    df["E"] = etr.astype(int)
    cph = CoxPHFitter(penalizer=0.1, l1_ratio=0.5)
    cph.fit(df, "T", "E")
    coef = cph.params_.values
    return {"idx": idx.tolist(), "mean": mu.tolist(), "scale": sd.tolist(),
            "coef": coef.tolist(), "n_nonzero": int((np.abs(coef) > 1e-8).sum())}


def cox_score(model, X):
    Z = X[:, np.asarray(model["idx"])]
    Z = (Z - np.asarray(model["mean"])) / np.asarray(model["scale"])
    return Z @ np.asarray(model["coef"])


def cindex(t, e, risk):
    """Harrell's C for a risk score (higher = worse)."""
    t, e, risk = np.asarray(t, float), np.asarray(e, bool), np.asarray(risk, float)
    num = den = 0.0
    for i in range(len(t)):
        if not e[i]:
            continue
        cmp = t > t[i]
        den += cmp.sum()
        num += (risk[cmp] < risk[i]).sum() + 0.5 * (risk[cmp] == risk[i]).sum()
    return num / den if den else float("nan")


# ---------------------------------------------------------------- STAGE 1
def stage1():
    run = Run("armBmeth-s1")
    cfg = cfg_dict()
    cfg_hash = run.start(cfg, [f"{R08}/clinical_merged.tsv", f"{R08}/beta_450k.npy",
                               f"{R08}/probes.txt", f"{R08}/samples.tsv", MAN450, MANEPIC])
    run.log("stage1_start", config_sha256=cfg_hash, note="TCGA only; no EPIC byte read")

    B, pnames, cases, luad, keap1, sex, os_m, os_e = load_tcga(run)

    # ---- probe overlap from the two official manifests
    p450 = manifest_probes(MAN450)
    pepic = manifest_probes(MANEPIC)
    ov = p450 & pepic
    pidx = {p: i for i, p in enumerate(pnames)}
    ov_in_model = [p for p in pnames if p in ov]
    ov_cols = np.array([pidx[p] for p in ov_in_model])
    overlap_stats = {
        "manifest_450k_loci": len(p450), "manifest_epic_loci": len(pepic),
        "manifest_overlap_loci": len(ov),
        "r08_retained_probes": len(pnames),
        "r08_retained_on_epic": len(ov_in_model),
        "r08_retained_on_epic_frac": round(len(ov_in_model) / len(pnames), 4)}
    run.log("probe_overlap", **overlap_stats)

    # ---- reproduction leg: R10 verbatim (full array, R10 seed)
    Xall, med_all = impute_median(B[:, luad].T.astype(np.float64))
    y_k = keap1[luad].astype(int)
    r_full, oof_full = nested_cv(Xall, y_k, R10_SEED)
    run.log("repro_keap1_full", **r_full)
    run.gate("G1_R10_reproduction",
             f">= {EXPECT['gate_repro_auroc']} (R10 reported {EXPECT['R10_keap1_auroc']:.4f})",
             round(r_full["auroc"], 6), r_full["auroc"] >= EXPECT["gate_repro_auroc"],
             "same inputs, same config, same pipeline as R10 m11")
    m_keap1_full = fit_full(Xall, y_k, SEED)

    # ---- overlap-restricted KEAP1 variant (transferable), study seed
    Xov = Xall[:, ov_cols]
    del Xall
    r_ov, oof_ov = nested_cv(Xov, y_k, SEED)
    run.log("keap1_overlap_nestedcv", **r_ov)
    m_keap1_ov = fit_full(Xov, y_k, SEED)

    # ---- sex control: R10 reproduction (full array) then the transferable overlap model
    sexmask = luad & np.isin(sex, ["Male", "Female"])
    Xsex, _ = impute_median(B[:, sexmask].T.astype(np.float64))
    del B
    y_s = (sex[sexmask] == "Male").astype(int)
    r_sex, oof_sex = nested_cv(Xsex, y_s, R10_SEED)
    run.log("repro_sex_full", **r_sex)
    run.gate("G2_R10_sex_control", "> 0.95", round(r_sex["auroc"], 6), r_sex["auroc"] > 0.95)
    Xsex_ov = Xsex[:, ov_cols]
    del Xsex
    r_sex_ov, _ = nested_cv(Xsex_ov, y_s, SEED)
    run.log("sex_overlap_nestedcv", **r_sex_ov)
    m_sex_ov = fit_full(Xsex_ov, y_s, SEED)
    del Xsex_ov

    for nm, m, space in (("M1_keap1_full", m_keap1_full, pnames),
                         ("M1_keap1_overlap", m_keap1_ov, ov_in_model),
                         ("PC_sex_overlap", m_sex_ov, ov_in_model)):
        sel = [space[i] for i in m["selected_idx"]]
        nz = [p for p, c in zip(sel, m["coef"]) if c != 0]
        m["selected_probes"] = sel
        m["nonzero_probes"] = nz
        m["n_selected_on_epic"] = int(sum(p in ov for p in sel))
        m["n_nonzero_on_epic"] = int(sum(p in ov for p in nz))
        run.log("model_fit", model=nm, C=m["C"], n_selected=m["n_selected"],
                n_nonzero=m["n_nonzero"], n_selected_on_epic=m["n_selected_on_epic"],
                n_nonzero_on_epic=m["n_nonzero_on_epic"])

    # ---- OS Cox head on the overlap space
    lu = np.where(luad)[0]
    ok = np.isfinite(os_m[lu]) & (os_m[lu] >= 0)
    tt, ee = os_m[lu][ok], os_e[lu][ok]
    Xc = Xov[ok]
    run.log("cox_cohort", n=int(ok.sum()), events=int(ee.sum()),
            dropped_no_os=int((~ok).sum()))
    cv_risk = np.full(len(tt), np.nan)
    for tr, te in KFold(5, shuffle=True, random_state=SEED).split(Xc):
        mm = fit_cox(Xc[tr], tt[tr], ee[tr])
        cv_risk[te] = cox_score(mm, Xc[te])
    c_cv = cindex(tt, ee, cv_risk)
    lo, hi = boot_ci(lambda ix: cindex(tt[ix], ee[ix], cv_risk[ix]), len(tt), SEED, B=1000)
    run.log("cox_cv", cindex=round(float(c_cv), 4), ci95=[round(lo, 4), round(hi, 4)],
            n=int(len(tt)), events=int(ee.sum()))
    m_cox = fit_cox(Xc, tt, ee)
    m_cox["probes"] = [ov_in_model[i] for i in m_cox["idx"]]
    m_cox["n_nonzero_on_epic"] = int(sum(p in ov for p, c in zip(m_cox["probes"], m_cox["coef"])
                                        if abs(c) > 1e-8))
    tr_risk = cox_score(m_cox, Xc)
    run.log("cox_full_fit", n_nonzero=m_cox["n_nonzero"],
            train_cindex_optimistic=round(float(cindex(tt, ee, tr_risk)), 4))

    os.makedirs(WORK, exist_ok=True)
    np.savez_compressed(f"{WORK}/tcga_medians.npz", full=med_all, overlap=med_all[ov_cols])
    art = {"overlap_stats": overlap_stats,
           "overlap_probes": ov_in_model,
           "all_probes_n": len(pnames),
           "models": {"M1_keap1_full": m_keap1_full, "M1_keap1_overlap": m_keap1_ov,
                      "PC_sex_overlap": m_sex_ov, "M2_os_cox": m_cox},
           "tcga_results": {"repro_keap1_full_nestedcv": r_full,
                            "repro_sex_full_nestedcv": r_sex,
                            "keap1_overlap_nestedcv": r_ov,
                            "sex_overlap_nestedcv": r_sex_ov,
                            "os_cox_cv_cindex": float(c_cv),
                            "os_cox_cv_cindex_ci95": [lo, hi],
                            "os_cox_n": int(len(tt)), "os_cox_events": int(ee.sum())},
           "config_sha256": cfg_hash, "run_dir": run.dir}
    json.dump(art, open(f"{WORK}/stage1_artifacts.json", "w"))
    np.savetxt(f"{WORK}/oof_keap1_repro_full.csv",
               np.column_stack([y_k, oof_full]), delimiter=",",
               header="keap1_mutant,oof_probability", comments="")
    np.savetxt(f"{WORK}/oof_keap1_overlap.csv",
               np.column_stack([y_k, oof_ov]), delimiter=",",
               header="keap1_mutant,oof_probability", comments="")
    np.savetxt(f"{WORK}/oof_sex_repro_full.csv",
               np.column_stack([y_s, oof_sex]), delimiter=",",
               header="male,oof_probability", comments="")
    np.savetxt(f"{WORK}/oof_os_cox_tcga.csv",
               np.column_stack([tt, ee.astype(int), cv_risk]), delimiter=",",
               header="os_months,event,cv_risk_score", comments="")
    run.write("stage1_results.json",
              {**{k: v for k, v in art.items() if k not in ("overlap_probes", "models")},
               "model_summaries": {
                   nm: {k: v for k, v in m.items()
                        if k not in ("selected_idx", "mean", "scale", "coef",
                                     "selected_probes", "nonzero_probes", "idx", "probes")}
                   for nm, m in art["models"].items()}})
    run.finalize()
    print("\nSTAGE1 RUN DIR:", run.dir)
    print("CONFIG SHA256:", cfg_hash)


# ---------------------------------------------------------------- STAGE 2 — unsealing event
def stage2():
    art = json.load(open(f"{WORK}/stage1_artifacts.json"))
    cfg_hash = art["config_sha256"]
    ev = {
        "cohort": "Anti-PD-1 EPIC methylation (81 pts)",
        "event_number": 2,
        "arm": "B-meth",
        "geo_series": "GSE115246",
        "gcs_prefix": "gs://heydonto-quantara-lungcdx/data-request-2026-08/labels/lung/"
                      "GSE115246/",
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config_sha256": cfg_hash,
        "config_path": "nsclc-rwpr-study/armB-meth/evidence/config.yaml",
        "models_frozen_before_unseal": sorted(art["models"].keys()),
        "endpoints_preregistered": ["anti-PD-1 response AUROC (primary)",
                                    "OS/PFS association of the Cox score (if deposited)",
                                    "measured minimum-detectable AUROC by permutation"],
        "gates_preregistered": ["G1_R10_reproduction", "G2_R10_sex_control",
                                "G3_preregistered_before_unseal", "G4_probe_overlap_reported",
                                "G5_epic_sex_transfer", "G6_permutation_null_centred",
                                "G7_no_silent_drops"],
        "power_statement_preregistered": "n=81 is expected to be unable to resolve a modest "
                                        "AUROC; the minimum-detectable AUROC at 80% power is "
                                        "reported alongside the point estimate and an "
                                        "inconclusive verdict is a pre-declared success mode",
        "tcga_training_numbers_at_seal": art["tcga_results"],
        "probe_overlap_at_seal": art["overlap_stats"],
        "access_mode": "READ-ONLY on data-request-2026-08/",
        "operator": "Reza Nehzati (agent-executed)",
        "sealed_scopes_untouched": ["nsclc-rwpr-study/sealed/"],
        "declared_prehash_reads": "GCS object listings and byte counts under the GSE115246 "
                                  "prefix; suppl/filelist.txt; the Illumina EPIC v1.0 "
                                  "manifest shipped inside GSE115246_RAW.tar (array design, "
                                  "not sample data). No sample metadata, no beta values and "
                                  "no outcome field was read before this event.",
        "note": "event_number is the protocol-assigned index (Arm B methylation track), not a "
                "line number in this file. Event #1 is the CPTAC/Arm-B-WSI exposure.",
    }
    # Append-only, race-safe: read the current object with its generation, append our line,
    # write back under an if-generation-match precondition. Arm B-WSI logs event #1 to the
    # same file from a parallel track, so a plain read-modify-write could silently drop it.
    tmp = f"{WORK}/unsealing_events.jsonl"
    for attempt in range(6):
        st = subprocess.run(["gsutil", "stat", GCS_PROV], capture_output=True, text=True)
        gen = "0"
        existing = ""
        if st.returncode == 0:
            for line in st.stdout.splitlines():
                if "Generation:" in line:
                    gen = line.split(":", 1)[1].strip()
            existing = subprocess.run(["gsutil", "cat", GCS_PROV],
                                      capture_output=True, text=True).stdout
        if existing and not existing.endswith("\n"):
            existing += "\n"
        if '"arm": "B-meth"' in existing and cfg_hash in existing:
            print("event #2 for arm B-meth already present; not duplicating")
            break
        open(tmp, "w").write(existing + json.dumps(ev) + "\n")
        r = subprocess.run(["gsutil", "-h", f"x-goog-if-generation-match:{gen}",
                            "cp", tmp, GCS_PROV], capture_output=True, text=True)
        if r.returncode == 0:
            print(f"appended under generation precondition {gen} (attempt {attempt + 1})")
            break
        print("precondition failed, retrying:", r.stderr.strip()[:200])
        time.sleep(2)
    else:
        raise SystemExit("could not append the unsealing event")
    print(json.dumps(ev, indent=2))
    json.dump(ev, open(f"{WORK}/unsealing_event_2.json", "w"), indent=2)
    print("\n--- unsealing_events.jsonl now reads ---")
    print(subprocess.run(["gsutil", "cat", GCS_PROV], capture_output=True, text=True).stdout)


# ---------------------------------------------------------------- STAGE 3 — the EPIC cohort
def parse_series_matrix(path, keep_probes):
    """Returns (meta_lines dict, sample_ids list, probe list, beta matrix probes x samples)."""
    meta, ids, rows, prows = {}, None, [], []
    with gzip.open(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!"):
                k, _, v = line.rstrip("\n").partition("\t")
                meta.setdefault(k, []).append([x.strip().strip('"') for x in v.split("\t")])
        hdr = next(fh).rstrip("\n").split("\t")
        ids = [h.strip().strip('"') for h in hdr[1:]]
        for line in fh:
            if line.startswith("!series_matrix_table_end"):
                break
            p, _, rest = line.rstrip("\n").partition("\t")
            p = p.strip().strip('"')
            if p not in keep_probes:
                continue
            # fields may be empty ("" = not reported), so parse positionally, never by
            # whitespace collapsing — that silently shifts a whole row.
            f = rest.split("\t")
            if len(f) < len(ids):
                f += [""] * (len(ids) - len(f))
            v = np.empty(len(ids))
            for j in range(len(ids)):
                t = f[j].strip().strip('"')
                try:
                    v[j] = float(t)
                except ValueError:
                    v[j] = np.nan
            prows.append(p)
            rows.append(v)
    Bm = np.vstack(rows) if rows else np.zeros((0, len(ids)))
    return meta, ids, prows, Bm


def apply_logistic(model, Xs, mode):
    """Xs: samples x (the model's own selected probes, in model order). NaN allowed."""
    mu = np.asarray(model["mean"]); sd = np.asarray(model["scale"])
    if mode == "T1_naive":
        Xs = np.where(np.isnan(Xs), mu, Xs)     # TCGA training mean of that probe
        Z = (Xs - mu) / sd
    else:                                       # T2_cohort_z  (primary)
        cm = np.nanmean(Xs, 0); cs = np.nanstd(Xs, 0)
        cm = np.where(np.isnan(cm), mu, cm)
        cs = np.where((~np.isfinite(cs)) | (cs < 1e-8), 1.0, cs)
        Z = (Xs - cm) / cs
        Z = np.where(np.isnan(Z), 0.0, Z)       # imputed at the cohort mean
        Z = Z * (sd / sd)                       # no-op; scale already unit in both spaces
    lp = Z @ np.asarray(model["coef"]) + model["intercept"]
    return lp, 1.0 / (1.0 + np.exp(-lp))


def apply_cox(model, Xs, mode):
    """Xs: samples x (model['probes'] in order)."""
    mu = np.asarray(model["mean"]); sd = np.asarray(model["scale"])
    if mode == "T1_naive":
        Xs = np.where(np.isnan(Xs), mu, Xs)
        Z = (Xs - mu) / sd
    else:
        cm = np.nanmean(Xs, 0); cs = np.nanstd(Xs, 0)
        cm = np.where(np.isnan(cm), mu, cm)
        cs = np.where((~np.isfinite(cs)) | (cs < 1e-8), 1.0, cs)
        Z = np.where(np.isnan((Xs - cm) / cs), 0.0, (Xs - cm) / cs)
    return Z @ np.asarray(model["coef"])


def build_X(probe_list, cache_probes, cache_beta):
    """samples x len(probe_list); NaN where the deposit has no such probe."""
    pos = {p: i for i, p in enumerate(cache_probes)}
    n = cache_beta.shape[1]
    X = np.full((n, len(probe_list)), np.nan)
    hit = 0
    for j, p in enumerate(probe_list):
        i = pos.get(p)
        if i is not None:
            X[:, j] = cache_beta[i]
            hit += 1
    return X, hit


def auroc_ci_perm(y, score, seed, n_perm=N_PERM, boot=2000):
    y = np.asarray(y, int); score = np.asarray(score, float)
    auc = float(roc_auc_score(y, score))
    lo, hi = boot_ci(lambda ix: roc_auc_score(y[ix], score[ix]) if len(set(y[ix])) == 2
                     else np.nan, len(y), seed, B=boot)
    rng = np.random.default_rng(seed)
    null = [float(roc_auc_score(rng.permutation(y), score)) for _ in range(n_perm)]
    null = np.array(null)
    # two-sided-aware: report the one-sided p for the observed direction and for |AUROC-0.5|
    p_one = float((np.sum(null >= auc) + 1) / (n_perm + 1))
    p_two = float((np.sum(np.abs(null - 0.5) >= abs(auc - 0.5)) + 1) / (n_perm + 1))
    se = (hi - lo) / (2 * 1.96)
    mde = 0.5 + (norm.ppf(0.975) + norm.ppf(0.80)) * se
    return {"auroc": auc, "ci95": [lo, hi], "n": int(len(y)), "positives": int(y.sum()),
            "perm_null_mean": float(null.mean()), "perm_null_sd": float(null.std()),
            "perm_null_p95": float(np.percentile(null, 95)),
            "perm_p_one_sided": p_one, "perm_p_two_sided": p_two,
            "n_permutations": int(n_perm),
            "se_from_bootstrap": float(se),
            "min_detectable_auroc_80pct_power": float(mde)}, null


def stage3():
    """Post-unsealing. Parses the EPIC deposit, caches the model-relevant betas, dumps the
    complete deposited metadata verbatim. Interpretation of the metadata vocabulary happens
    in b2_eval.py so the mapping step is separately auditable."""
    run = Run("armBmeth-s3-parse")
    art = json.load(open(f"{WORK}/stage1_artifacts.json"))
    cfg_hash = art["config_sha256"]
    ev_local = json.load(open(f"{WORK}/unsealing_event_2.json"))
    remote = subprocess.run(["gsutil", "cat", GCS_PROV], capture_output=True, text=True).stdout
    logged = (cfg_hash in remote) and ('"event_number": 2' in remote)
    run.start(cfg_dict(), [SERIES_MATRIX, MANEPIC, f"{WORK}/stage1_artifacts.json"])
    run.gate("G3_preregistered_before_unseal",
             "config sha256 present in gs://.../unsealing_events.jsonl before label read",
             {"config_sha256": cfg_hash, "found_in_remote": logged,
              "event_ts": ev_local["timestamp_utc"]}, logged)

    models = art["models"]
    need = set(art["overlap_probes"]) | set(models["M1_keap1_full"]["selected_probes"]) \
        | set(models["M2_os_cox"]["probes"])
    meta, ids, prows, Bm = parse_series_matrix(SERIES_MATRIX, need)
    run.log("epic_deposit", n_samples=len(ids), probes_requested=len(need),
            probes_found=len(prows), meta_keys=sorted(meta.keys()))
    json.dump(meta, open(f"{run.dir}/epic_metadata_raw.json", "w"), indent=2)
    json.dump(meta, open(f"{WORK}/epic_metadata_dump.json", "w"), indent=2)
    np.savez_compressed(f"{WORK}/epic_betas.npz", beta=Bm.astype(np.float32),
                        probes=np.array(prows), samples=np.array(ids))
    summ = {"n_samples": len(ids), "probes_found": len(prows),
            "probes_requested": len(need), "meta_keys": sorted(meta.keys()),
            "value_range": [float(np.nanmin(Bm)), float(np.nanmax(Bm))],
            "nan_frac": float(np.isnan(Bm).mean())}
    json.dump(summ, open(f"{WORK}/stage3_parse_summary.json", "w"), indent=2)
    run.write("stage3_parse.json", summ)
    run.finalize()
    print("\nSTAGE3 RUN DIR:", run.dir)
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    s = sys.argv[1] if len(sys.argv) > 1 else "stage1"
    {"stage1": stage1, "stage2": stage2, "stage3": stage3}[s]()
