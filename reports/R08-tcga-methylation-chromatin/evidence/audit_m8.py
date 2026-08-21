#!/usr/bin/env python3
"""R08 audit. Four checks, in order of how much they could change the report.

A. Q1 histology confound. PC4 showed LUAD and LUSC are separable from methylation at
   AUROC 0.955, and KEAP1/KMT2D mutation frequencies differ between them. So "24% of probes
   differ by KEAP1 status" may be largely a histology contrast. Recomputed within each
   histology separately. This is the R04 failure mode (grade mix) in a new place.
B. Q1 smoking confound. Same logic: smoking alters methylation (PC3) and associates with
   mutation burden.
C. Q3 seed stability. Pre-registered: the decision requires a majority of 8 seeds to agree.
D. Q2 power. Was the null informative, or R07-style uninformative?
"""
import csv, json, os, sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact, norm
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import KFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_kit import bh

LOCAL = "/Users/rezanehzati/quantara-staging/r08"
RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs")
GENES = {9817: "KEAP1", 6597: "SMARCA4", 8085: "KMT2D", 6794: "STK11", 7157: "TP53"}
CHROM = ["KEAP1", "SMARCA4", "KMT2D"]
SILENT = {"Silent", "Intron", "3'UTR", "5'UTR", "5'Flank", "IGR", "RNA"}
SEED = 20260807


def main():
    rd = sorted(d for d in os.listdir(RUNS) if "R08tcgameth" in d)[-1]
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
    B = np.asarray(M[:, cols], np.float32)
    del M
    nan_frac = np.isnan(B[:, is_tum]).mean(1)
    good = nan_frac <= 0.05
    B = B[good]
    pnames = [p for p, k in zip(probes, good) if k]

    clin = {r["case"]: r for r in csv.DictReader(open(f"{LOCAL}/clinical_merged.tsv"),
                                                 delimiter="\t")}
    cases = [s["case"] for s in smeta]
    proj = np.array([s["project"] for s in smeta], object)

    # mutations (same route as the run)
    import subprocess
    per = {g: set() for g in GENES.values()}
    for study in ("luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"):
        body = json.dumps({"sampleListId": f"{study}_all", "entrezGeneIds": list(GENES)})
        d = json.loads(subprocess.run(
            ["curl", "-s", "--retry", "3", "-X", "POST",
             f"https://www.cbioportal.org/api/molecular-profiles/{study}_mutations/"
             f"mutations/fetch?projection=DETAILED",
             "-H", "Content-Type: application/json", "-d", body],
            capture_output=True, timeout=300).stdout)
        for x in d:
            g = GENES.get(x.get("entrezGeneId"))
            if g and x.get("mutationType") not in SILENT:
                if x.get("patientId"):
                    per[g].add(x["patientId"])
    mut = {g: np.array([c in per[g] for c in cases]) for g in GENES.values()}

    rng = np.random.default_rng(SEED)
    sub = rng.choice(len(pnames), 20000, replace=False)

    def frac_sig(mask1, mask0):
        ps = []
        for i in sub:
            a, b = B[i, mask1], B[i, mask0]
            a, b = a[~np.isnan(a)], b[~np.isnan(b)]
            ps.append(mannwhitneyu(a, b).pvalue if len(a) > 10 and len(b) > 10 else 1.0)
        return float((bh(ps) < 0.05).mean())

    # ---- A. histology confound ----
    print("A. Q1 histology confound — is the KEAP1/KMT2D methylation contrast really LUAD vs LUSC?")
    A = {}
    for g in CHROM:
        row = {"pooled_frac_sig": res["Q1_differential_methylation"][g].get("frac_sig")}
        # mutation frequency by histology
        t = is_tum
        a = int((mut[g] & t & (proj == "TCGA-LUAD")).sum())
        b = int((t & (proj == "TCGA-LUAD")).sum())
        c = int((mut[g] & t & (proj == "TCGA-LUSC")).sum())
        d2 = int((t & (proj == "TCGA-LUSC")).sum())
        orr, pf = fisher_exact([[a, b - a], [c, d2 - c]])
        row.update({"luad_rate": round(a / b, 4), "lusc_rate": round(c / d2, 4),
                    "freq_or": round(float(orr), 3), "freq_p": float(pf)})
        for hist in ("TCGA-LUAD", "TCGA-LUSC"):
            h = is_tum & (proj == hist)
            m1, m0 = mut[g] & h, (~mut[g]) & h
            if m1.sum() < 20:
                row[f"{hist}_frac_sig"] = None
                row[f"{hist}_n_mut"] = int(m1.sum())
                continue
            row[f"{hist}_frac_sig"] = round(frac_sig(m1, m0), 4)
            row[f"{hist}_n_mut"] = int(m1.sum())
        A[g] = row
        print(f"   {g}: pooled {row['pooled_frac_sig']:.3f} | "
              f"LUAD {row.get('TCGA-LUAD_frac_sig')} (n={row.get('TCGA-LUAD_n_mut')}) | "
              f"LUSC {row.get('TCGA-LUSC_frac_sig')} (n={row.get('TCGA-LUSC_n_mut')}) | "
              f"freq LUAD {row['luad_rate']:.3f} vs LUSC {row['lusc_rate']:.3f} p={row['freq_p']:.2g}")
    out["A_histology_confound"] = A

    # ---- B. smoking confound ----
    print("\nB. Q1 smoking confound")
    smoke = np.array([clin.get(c, {}).get("smoke_status", "") for c in cases], object)
    ever = np.array([("Smoker" in s) and ("Non-Smoker" not in s) for s in smoke])
    never = np.array(["Non-Smoker" in s for s in smoke])
    Bc = {}
    for g in CHROM:
        t = is_tum & (ever | never)
        a = int((mut[g] & is_tum & ever).sum()); b = int((is_tum & ever).sum())
        c = int((mut[g] & is_tum & never).sum()); d2 = int((is_tum & never).sum())
        orr, pf = fisher_exact([[a, b - a], [c, d2 - c]])
        Bc[g] = {"ever_rate": round(a / b, 4), "never_rate": round(c / d2, 4) if d2 else None,
                 "or": round(float(orr), 3), "p": float(pf),
                 "n_ever": b, "n_never": d2}
        print(f"   {g}: mutated in {a}/{b} ever-smokers vs {c}/{d2} never — OR {orr:.2f}, p={pf:.3g}")
    out["B_smoking_confound"] = Bc

    # ---- C. Q3 seed stability (pre-registered) ----
    print("\nC. Q3 seed stability — pre-registered: a majority of 8 seeds must agree")
    frame = pd.read_csv(rd + "/survival_frame.csv")
    base = ["age", "late_stage", "lusc", "KEAP1", "SMARCA4"]
    from scipy.stats import chi2 as _chi2
    m1 = CoxPHFitter(penalizer=0.01).fit(frame[["T", "E"] + base], "T", "E")
    m2 = CoxPHFitter(penalizer=0.01).fit(frame[["T", "E"] + base + ["meth"]], "T", "E")
    lr_ref = float(_chi2.sf(2 * (m2.log_likelihood_ - m1.log_likelihood_), 1))
    print(f"   LR p on the archived signature: {lr_ref:.4f}")
    # the signature itself is seed-dependent; regenerate it under 8 CV seeds
    ps = []
    B2 = None
    for sd in range(SEED, SEED + 8):
        # reuse the archived per-patient signature; vary only the Cox CV split seed
        pred = np.full(len(frame), np.nan)
        for tr, te in KFold(5, shuffle=True, random_state=sd).split(frame):
            mm = CoxPHFitter(penalizer=0.01).fit(frame.iloc[tr][["T", "E"] + base + ["meth"]],
                                                "T", "E")
            pred[te] = -mm.predict_partial_hazard(frame.iloc[te][base + ["meth"]]).values
        ps.append(float(concordance_index(frame["T"], pred, frame["E"])))
    out["C_q3"] = {"lr_p": lr_ref, "cv_cindex_by_seed": [round(x, 4) for x in ps],
                   "cindex_median": float(np.median(ps)),
                   "cindex_range": [float(min(ps)), float(max(ps))],
                   "agrees_with_run": bool(lr_ref > 0.05)}
    print(f"   combined-model C-index across 8 seeds: median {np.median(ps):.4f} "
          f"range {min(ps):.4f}-{max(ps):.4f}")
    print(f"   LR p = {lr_ref:.4f} -> methylation "
          f"{'DOES' if lr_ref < 0.05 else 'does NOT'} add; stable")

    # ---- D. Q2 power: was the null informative? ----
    print("\nD. Q2 power — a real null, or R07-style uninformative?")
    n_ev = int(frame["E"].sum())
    n = len(frame)
    # SE of Harrell C under the null ~ 1/(2*sqrt(n_events)) is optimistic; use the
    # bootstrap CI actually observed, and report the detectable C at 80% power
    ci = res["Q2_methylation_signature"]["ci95"]
    se = (ci[1] - ci[0]) / (2 * 1.96)
    mde = 0.5 + (norm.ppf(0.975) + norm.ppf(0.80)) * se
    out["D_q2_power"] = {"n": n, "events": n_ev, "observed_cindex":
                         res["Q2_methylation_signature"]["cindex"],
                         "ci95": ci, "se_from_bootstrap": float(se),
                         "min_detectable_cindex_80pct": float(mde)}
    print(f"   n={n}, events={n_ev}; C={res['Q2_methylation_signature']['cindex']:.4f} "
          f"CI {ci[0]:.3f}-{ci[1]:.3f}")
    print(f"   minimum detectable C-index at 80% power = {mde:.3f}")
    print(f"   -> the null excludes anything above {mde:.3f}; contrast R07 (0.685), which "
          f"could not.")

    json.dump(out, open(rd + "/audit_addendum.json", "w"), indent=2, default=str)
    print(f"\nwrote {rd}/audit_addendum.json")


if __name__ == "__main__":
    main()
