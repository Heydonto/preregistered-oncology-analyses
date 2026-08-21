#!/usr/bin/env python3
"""R10 audit. Four adversarial checks, each aimed at the way a result could be spurious.

A. Q3's NRF2 enrichment could be an artefact of gene size. Random sets were matched on the
   NUMBER of genes, not on probes-per-gene. If NRF2 genes carry more probes and probe count
   correlates with the significant fraction, the enrichment is a size effect. Re-run the null
   matched on probe count.
B. Q4's overlap needs a chance baseline. "12.4% shared" means nothing without the fraction two
   independent signatures of these sizes would share anyway.
C. Q1's classifier could be reading smoking rather than KEAP1. Retrain after DELETING every
   probe shared with the smoking signature; if AUROC holds, the signal is its own.
D. Q1's AUROC could be a fold-split artefact (this is what changed R05's conclusion). Re-run
   under 5 seeds.
"""
import csv, gzip, json, os, subprocess, sys
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_kit import bh

LOCAL = "/Users/rezanehzati/quantara-staging/r08"
MANIFEST = "/tmp/manifest450k.csv.gz"
RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs")
GENES = {9817: "KEAP1", 6597: "SMARCA4", 8085: "KMT2D", 6794: "STK11", 7157: "TP53"}
SILENT = {"Silent", "Intron", "3'UTR", "5'UTR", "5'Flank", "IGR", "RNA"}
SEED = 20260807
NRF2 = {"NQO1", "HMOX1", "GCLC", "GCLM", "GSR", "TXN", "TXNRD1", "SLC7A11", "SRXN1", "PRDX1",
        "GPX2", "G6PD", "PGD", "TALDO1", "IDH1", "ME1", "AKR1C1", "AKR1C2", "AKR1B10",
        "ABCC1", "ABCC2", "ABCC3", "FTL", "FTH1", "CBR1", "GSTP1", "GSTM1", "GSTA1",
        "UGT1A1", "CAT", "SOD1", "EPHX1", "CES1", "OSGIN1", "TXNDC17",
        "KEAP1", "NFE2L2", "CUL3", "SQSTM1"}


def main():
    rd = sorted(d for d in os.listdir(RUNS) if "R10keap1pheno" in d)[-1]
    rd = os.path.join(RUNS, rd)
    res = json.load(open(rd + "/results.json"))
    print(f"auditing {os.path.basename(rd)}\n")
    out = {"run_dir": rd}

    M = np.load(f"{LOCAL}/beta_450k.npy", mmap_mode="r")
    probes = open(f"{LOCAL}/probes.txt").read().split()
    samples = list(csv.DictReader(open(f"{LOCAL}/samples.tsv"), delimiter="\t"))
    by = {}
    for s in samples:
        if s["sample"] not in by or s["file_id"] < by[s["sample"]]["file_id"]:
            by[s["sample"]] = s
    cols = sorted(int(s["col"]) for s in by.values())
    smeta = [samples[c] for c in cols]
    is_tum = np.array([s["sample_type"] != "Solid Tissue Normal" for s in smeta])
    B = np.asarray(M[:, cols], np.float32); del M
    good = np.isnan(B[:, is_tum]).mean(1) <= 0.05
    B = B[good]
    pnames = [p for p, k in zip(probes, good) if k]
    clin = {r["case"]: r for r in csv.DictReader(open(f"{LOCAL}/clinical_merged.tsv"),
                                                delimiter="\t")}
    cases = [s["case"] for s in smeta]
    proj = np.array([s["project"] for s in smeta], object)
    luad = is_tum & (proj == "TCGA-LUAD")
    per = set()
    for study in ("luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"):
        d = json.loads(subprocess.run(
            ["curl", "-s", "--retry", "3", "-X", "POST",
             f"https://www.cbioportal.org/api/molecular-profiles/{study}_mutations/"
             f"mutations/fetch?projection=DETAILED", "-H", "Content-Type: application/json",
             "-d", json.dumps({"sampleListId": f"{study}_all", "entrezGeneIds": [9817]})],
            capture_output=True, timeout=300).stdout)
        per |= {x["patientId"] for x in d
                if x.get("mutationType") not in SILENT and x.get("patientId")}
    keap1 = np.array([c in per for c in cases])
    smoke = np.array([clin.get(c, {}).get("smoke_status", "") for c in cases], object)

    def qvals(m1, m0):
        r = mannwhitneyu(B[:, m1], B[:, m0], axis=1, nan_policy="omit")
        return bh(np.nan_to_num(r.pvalue, nan=1.0))

    sig_k = qvals(keap1 & luad, (~keap1) & luad) < 0.05
    wt = (~keap1) & luad
    ever = np.array([("Smoker" in s) and ("Non-Smoker" not in s) for s in smoke])
    never = np.array(["Non-Smoker" in s for s in smoke])
    sig_s = qvals(wt & ever, wt & never) < 0.05

    # ---- A. NRF2 enrichment, probe-count-matched null ----
    print("A. Q3 NRF2 enrichment — is it a gene-size artefact?")
    ann = {}
    with gzip.open(MANIFEST, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("IlmnID"):
                break
        need = set(pnames)
        for line in fh:
            p = line.split(",")
            if len(p) < 26 or not p[0].startswith("cg") or p[0] not in need:
                continue
            ann[p[0]] = {g for g in p[21].split(";") if g}
    hits, tot = {}, {}
    for p, s in zip(pnames, sig_k):
        for g in ann.get(p, ()):
            tot[g] = tot.get(g, 0) + 1
            if s:
                hits[g] = hits.get(g, 0) + 1
    genes = [g for g in tot if tot[g] >= 5]
    frac = {g: hits.get(g, 0) / tot[g] for g in genes}
    npres = sorted(set(genes) & NRF2)
    obs = float(np.mean([frac[g] for g in npres]))
    rho, prho = spearmanr([tot[g] for g in genes], [frac[g] for g in genes])
    print(f"   probes-per-gene vs frac_sig: Spearman rho={rho:.4f} (p={prho:.2g})")
    # probe-count-matched sampling: for each NRF2 gene draw a gene with similar probe count
    rng = np.random.default_rng(SEED)
    bygroup = {}
    for g in genes:
        bygroup.setdefault(min(tot[g], 60) // 5, []).append(g)
    null_m = []
    for _ in range(1000):
        pick = []
        for g in npres:
            pool = bygroup[min(tot[g], 60) // 5]
            pick.append(pool[rng.integers(len(pool))])
        null_m.append(float(np.mean([frac[x] for x in pick])))
    null_m = np.array(null_m)
    p_m = float((np.sum(null_m >= obs) + 1) / 1001)
    print(f"   NRF2 observed {obs:.4f}")
    print(f"   unmatched null (report)      mean {res['Q3_nrf2_enrichment']['random_null_mean']:.4f}"
          f"  p={res['Q3_nrf2_enrichment']['p_empirical']:.4f}")
    print(f"   probe-count-matched null     mean {null_m.mean():.4f}  p={p_m:.4f}")
    print(f"   -> {'SURVIVES' if p_m < 0.05 else 'DOES NOT SURVIVE'} size matching")
    out["A_nrf2"] = {"observed": obs, "matched_null_mean": float(null_m.mean()),
                     "matched_p": p_m, "survives": bool(p_m < 0.05),
                     "spearman_probes_vs_frac": [float(rho), float(prho)],
                     "n_nrf2_genes": len(npres)}

    # ---- B. Q4 overlap against a chance baseline ----
    print("\nB. Q4 overlap — what would independent signatures share by chance?")
    n = len(pnames)
    exp = float(sig_k.mean() * sig_s.mean() * n)
    obs_ov = int((sig_k & sig_s).sum())
    maxpos = int(min(sig_k.sum(), sig_s.sum()))
    print(f"   KEAP1 signature {int(sig_k.sum()):,} probes; smoking signature {int(sig_s.sum()):,}")
    print(f"   observed overlap {obs_ov:,};  expected if independent {exp:,.0f}"
          f";  maximum possible {maxpos:,}")
    print(f"   observed/expected = {obs_ov/exp:.2f}x    "
          f"observed/maximum = {obs_ov/maxpos:.2f}")
    print(f"   -> smoking accounts for {obs_ov/int(sig_k.sum())*100:.1f}% of the KEAP1 "
          f"signature, only {obs_ov/exp:.2f}x chance")
    out["B_overlap"] = {"keap1_size": int(sig_k.sum()), "smoking_size": int(sig_s.sum()),
                        "observed": obs_ov, "expected_if_independent": exp,
                        "max_possible": maxpos, "obs_over_exp": obs_ov / exp,
                        "frac_of_keap1": obs_ov / int(sig_k.sum()),
                        "n_never_smokers_wt": int((wt & never).sum()),
                        "n_keap1_mutant": int((keap1 & luad).sum())}

    # ---- C. classifier with smoking-shared probes deleted ----
    print("\nC. Q1 classifier — does it survive deleting every smoking-shared probe?")

    def run_clf(mask_probes, y, mask_s, seed):
        X = B[mask_probes][:, mask_s].T.astype(np.float64)
        X = np.where(np.isnan(X), np.nanmedian(X, 0), X)
        yy = y[mask_s].astype(int)
        pipe = Pipeline([("sel", SelectKBest(f_classif, k=5000)), ("sc", StandardScaler()),
                         ("clf", LogisticRegression(penalty="elasticnet", solver="saga",
                                                    l1_ratio=0.5, max_iter=3000, tol=1e-3,
                                                    random_state=SEED))])
        oof = np.full(len(yy), np.nan)
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(X, yy):
            gs = GridSearchCV(pipe, {"clf__C": [0.01, 0.1, 1.0]}, scoring="roc_auc",
                              cv=StratifiedKFold(3, shuffle=True, random_state=SEED), n_jobs=-1)
            gs.fit(X[tr], yy[tr])
            oof[te] = gs.predict_proba(X[te])[:, 1]
        return float(roc_auc_score(yy, oof))

    keep = ~sig_s
    auc_nosmoke = run_clf(keep, keap1, luad, SEED)
    print(f"   all probes            AUROC {res['Q1_keap1_classifier']['auroc']:.4f}")
    print(f"   smoking probes removed AUROC {auc_nosmoke:.4f}  "
          f"({int(keep.sum()):,} of {n:,} probes retained)")
    out["C_no_smoking_probes"] = {"auroc_all": res["Q1_keap1_classifier"]["auroc"],
                                 "auroc_smoking_removed": auc_nosmoke,
                                 "probes_retained": int(keep.sum()),
                                 "holds": bool(auc_nosmoke > 0.85)}

    # ---- D. seed stability ----
    print("\nD. Q1 seed stability (this is what changed R05's conclusion)")
    aucs = [run_clf(np.ones(n, bool), keap1, luad, s) for s in range(SEED, SEED + 5)]
    print(f"   AUROC across 5 seeds: {[round(a,4) for a in aucs]}")
    print(f"   median {np.median(aucs):.4f}  range {min(aucs):.4f}-{max(aucs):.4f}")
    out["D_seed_stability"] = {"aucs": [float(a) for a in aucs],
                               "median": float(np.median(aucs)),
                               "range": [float(min(aucs)), float(max(aucs))]}

    json.dump(out, open(rd + "/audit_addendum.json", "w"), indent=2, default=str)
    print(f"\nwrote {rd}/audit_addendum.json")


if __name__ == "__main__":
    main()
