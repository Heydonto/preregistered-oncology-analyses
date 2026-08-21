#!/usr/bin/env python3
"""R10 — what IS the KEAP1 methylation phenotype? Characterising R08's one positive finding.

R08 established that 29.4% of tested probes are differentially methylated by KEAP1 mutation
status within lung adenocarcinoma. That is a statistical statement with no biological content:
it does not say how strong the effect is as a predictor, where in the genome it sits, which
genes are involved, or whether it is really a smoking signature wearing a KEAP1 badge.

This report answers those four questions on all 392,489 retained probes (R08 sampled 20,000).

Q1 HOW STRONG   nested-CV classifier of KEAP1 status from methylation, within LUAD. Converts
                "% of probes" into an AUROC with a confidence interval.
                Positive control: sex through the identical pipeline, expected near 1.0.
Q2 WHERE        are differential probes enriched at CpG islands and promoters, or in open sea?
                Tested against the array's own background composition.
Q3 WHICH GENES  a pre-declared canonical NRF2/KEAP1-pathway gene set, tested against 1,000
                size-matched random gene sets.
                THE PRIOR HERE IS WEAK AND IS DECLARED AS SUCH: KEAP1 loss activates NRF2 at
                the protein level by ending its degradation. That does not entail a promoter
                methylation change at NRF2 targets. A null result is informative, not a
                failure, and the unbiased top genes are reported regardless.
Q4 CONFOUNDED?  R08 declared the smoking confound unresolvable, because only 4 of 81
                never-smokers carry a KEAP1 mutation, so no never-smoker stratum can support a
                genome-wide test. This report attacks it from the other side: define the
                smoking signature in KEAP1-WILDTYPE tumours only, then ask how much of the
                KEAP1 signature is the same probes. Heavy overlap would mean the KEAP1
                phenotype is largely a smoking phenotype. This does not need never-smokers who
                carry the mutation.

Q4 is the point of the report. Q1-Q3 describe the phenotype; Q4 tests whether it is its own.
"""
import csv, gzip, json, os, subprocess, sys
import numpy as np
from scipy.stats import mannwhitneyu, fisher_exact
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_kit import Run, bh, boot_ci

LOCAL = "/Users/rezanehzati/quantara-staging/r08"
MANIFEST = "/tmp/manifest450k.csv.gz"
GENES = {9817: "KEAP1", 6597: "SMARCA4", 8085: "KMT2D", 6794: "STK11", 7157: "TP53"}
SILENT = {"Silent", "Intron", "3'UTR", "5'UTR", "5'Flank", "IGR", "RNA"}
SEED = 20260807
# Pre-declared canonical NRF2 (NFE2L2) targets and KEAP1-pathway members. Fixed before any
# enrichment was computed; hashed into config.yaml.
NRF2 = {"NQO1", "HMOX1", "GCLC", "GCLM", "GSR", "TXN", "TXNRD1", "SLC7A11", "SRXN1", "PRDX1",
        "GPX2", "G6PD", "PGD", "TALDO1", "IDH1", "ME1", "AKR1C1", "AKR1C2", "AKR1B10",
        "ABCC1", "ABCC2", "ABCC3", "FTL", "FTH1", "CBR1", "GSTP1", "GSTM1", "GSTA1",
        "UGT1A1", "CAT", "SOD1", "EPHX1", "CES1", "OSGIN1", "TXNDC17",
        "KEAP1", "NFE2L2", "CUL3", "SQSTM1"}
EXPECT = {"probes": 392489, "luad_tumours": 471, "keap1_mutant_luad": 84}


def load_annotation(probes_needed):
    """probe -> (genes, island_relation, refgene_group) from the Illumina v1.1 manifest."""
    ann = {}
    with gzip.open(MANIFEST, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("IlmnID"):
                break
        for line in fh:
            p = line.split(",")
            if len(p) < 26 or not p[0].startswith("cg"):
                continue
            if p[0] not in probes_needed:
                continue
            genes = {g for g in p[21].split(";") if g}
            ann[p[0]] = (genes, p[25].strip() or "OpenSea", p[23].split(";")[0] if p[23] else "")
    return ann


def main():
    run = Run("R10keap1pheno")
    cfg = {
        "report": "R10",
        "question": "what is the KEAP1 methylation phenotype R08 found, and is it its own?",
        "cohort": "TCGA-LUAD primary tumours only (LUSC excluded: R08 showed histology "
                  "confounds every pooled comparison)",
        "probes": "ALL 392,489 retained probes, not R08's 20,000 sample",
        "Q1": "nested-CV classifier of KEAP1 status; positive control = sex, same pipeline",
        "Q2": "CpG-island and gene-region enrichment against the array's own background",
        "Q3": "pre-declared NRF2/KEAP1-pathway gene set vs 1,000 size-matched random sets",
        "Q3_declared_weak_prior": "KEAP1 loss activates NRF2 post-translationally; this does "
                                  "NOT entail promoter methylation change at its targets. A "
                                  "null is informative and the unbiased top genes are reported "
                                  "either way",
        "Q4": "define the smoking signature in KEAP1-WILDTYPE tumours only, then measure its "
              "overlap with the KEAP1 signature. Attacks R08's declared-unresolvable smoking "
              "confound without needing never-smokers who carry the mutation",
        "gate_rule": "the sex positive control must reach AUROC>0.95 through the same "
                     "classifier pipeline, else the machinery is not trustworthy",
        "nrf2_gene_set": sorted(NRF2),
        "annotation": "Illumina HumanMethylation450 v1.1 manifest via GEO GPL13534 supplementary",
        "seed": SEED, "expectations": EXPECT,
        "relates_to": ["R08 KEAP1 29.4% of probes within LUAD, min q=1.7e-20",
                       "R08 limitation 3: smoking confound declared unresolvable"],
    }
    cfg_hash = run.start(cfg, [f"{LOCAL}/clinical_merged.tsv", MANIFEST])

    # ---- rebuild R08's matrix and filters exactly ----
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
    good = np.isnan(B[:, is_tum]).mean(1) <= 0.05
    B = B[good]
    pnames = [p for p, k in zip(probes, good) if k]
    run.gate("G0_probes", EXPECT["probes"], len(pnames), len(pnames) == EXPECT["probes"],
             "identical filter to R08, so the two reports describe the same probe set")

    clin = {r["case"]: r for r in csv.DictReader(open(f"{LOCAL}/clinical_merged.tsv"),
                                                delimiter="\t")}
    cases = [s["case"] for s in smeta]
    proj = np.array([s["project"] for s in smeta], object)
    luad = is_tum & (proj == "TCGA-LUAD")
    run.gate("G1_luad", EXPECT["luad_tumours"], int(luad.sum()), int(luad.sum()) == EXPECT["luad_tumours"])

    # mutations
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

    smoke = np.array([clin.get(c, {}).get("smoke_status", "") for c in cases], object)
    sex = np.array([clin.get(c, {}).get("sex", "") for c in cases], object)

    def diff_probes(mask1, mask0, label):
        """Mann-Whitney across all probes at once; returns BH q-values."""
        a, b = B[:, mask1], B[:, mask0]
        res = mannwhitneyu(a, b, axis=1, nan_policy="omit")
        q = bh(np.nan_to_num(res.pvalue, nan=1.0))
        run.log(label, n1=int(mask1.sum()), n0=int(mask0.sum()),
                n_sig=int((q < 0.05).sum()), frac=round(float((q < 0.05).mean()), 4))
        return q

    # ---- the KEAP1 signature, all probes, within LUAD ----
    q_keap1 = diff_probes(keap1 & luad, (~keap1) & luad, "signature_KEAP1")
    sig_keap1 = q_keap1 < 0.05
    run.gate("G3_signature_reproduces_R08", "frac significant within 5pp of R08's 0.294",
             round(float(sig_keap1.mean()), 4), abs(float(sig_keap1.mean()) - 0.294) < 0.05,
             "all 392,489 probes vs R08's 20,000 sample")

    # ---- Q1 classifier + sex positive control ----
    def classify(y, mask, label):
        X = B[:, mask].T.astype(np.float64)
        X = np.where(np.isnan(X), np.nanmedian(X, 0), X)
        yy = y[mask].astype(int)
        pipe = Pipeline([("sel", SelectKBest(f_classif, k=5000)), ("sc", StandardScaler()),
                         ("clf", LogisticRegression(penalty="elasticnet", solver="saga",
                                                    l1_ratio=0.5, max_iter=3000, tol=1e-3,
                                                    random_state=SEED))])
        oof = np.full(len(yy), np.nan)
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED).split(X, yy):
            gs = GridSearchCV(pipe, {"clf__C": [0.01, 0.1, 1.0]}, scoring="roc_auc",
                              cv=StratifiedKFold(3, shuffle=True, random_state=SEED), n_jobs=-1)
            gs.fit(X[tr], yy[tr])
            oof[te] = gs.predict_proba(X[te])[:, 1]
        auc = float(roc_auc_score(yy, oof))
        lo, hi = boot_ci(lambda ix: roc_auc_score(yy[ix], oof[ix]) if len(set(yy[ix])) == 2
                         else np.nan, len(yy), SEED, B=2000)
        run.log(label, auroc=round(auc, 4), ci95=[round(lo, 4), round(hi, 4)],
                n=len(yy), positives=int(yy.sum()))
        return {"auroc": auc, "ci95": [lo, hi], "n": len(yy), "positives": int(yy.sum())}, oof

    sexmask = luad & np.isin(sex, ["Male", "Female"])
    pc_sex, _ = classify((sex == "Male"), sexmask, "PC_sex_classifier")
    run.gate("G4_positive_control_sex", "AUROC > 0.95", round(pc_sex["auroc"], 4),
             pc_sex["auroc"] > 0.95, "same pipeline as the KEAP1 classifier")
    q1, oof_keap1 = classify(keap1, luad, "Q1_keap1_classifier")

    # ---- annotation-dependent analyses ----
    ann = load_annotation(set(pnames))
    run.gate("G5_annotation_coverage", "> 90% of probes annotated",
             f"{len(ann)}/{len(pnames)}", len(ann) / len(pnames) > 0.90)
    isl = np.array([ann.get(p, (set(), "OpenSea", ""))[1] for p in pnames], object)
    grp = np.array([ann.get(p, (set(), "", ""))[2] for p in pnames], object)

    # ---- Q2 where in the genome ----
    def enrich(cats, arr):
        out = {}
        for c in cats:
            inc = arr == c
            a = int((inc & sig_keap1).sum()); b = int((inc & ~sig_keap1).sum())
            cc = int((~inc & sig_keap1).sum()); d = int((~inc & ~sig_keap1).sum())
            orr, p = fisher_exact([[a, b], [cc, d]])
            out[c] = {"n_probes": int(inc.sum()), "n_sig": a,
                      "frac_sig": round(a / max(inc.sum(), 1), 4),
                      "odds_ratio": round(float(orr), 3), "p": float(p)}
        return out
    q2 = {"cpg_island_context": enrich(["Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf",
                                        "OpenSea"], isl),
          "gene_region": enrich(["TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR", ""],
                                grp),
          "overall_frac_sig": round(float(sig_keap1.mean()), 4)}
    for k, v in q2["cpg_island_context"].items():
        run.log("Q2_island", context=k, frac_sig=v["frac_sig"], odds_ratio=v["odds_ratio"])

    # ---- Q3 gene-level: NRF2 set vs random sets ----
    gene_hits, gene_tot = {}, {}
    for p, s in zip(pnames, sig_keap1):
        for g in ann.get(p, (set(),))[0]:
            gene_tot[g] = gene_tot.get(g, 0) + 1
            if s:
                gene_hits[g] = gene_hits.get(g, 0) + 1
    genes_all = [g for g in gene_tot if gene_tot[g] >= 5]
    frac = {g: gene_hits.get(g, 0) / gene_tot[g] for g in genes_all}
    nrf2_present = sorted(set(genes_all) & NRF2)
    obs = float(np.mean([frac[g] for g in nrf2_present])) if nrf2_present else None
    rng = np.random.default_rng(SEED)
    null = [float(np.mean([frac[g] for g in rng.choice(genes_all, len(nrf2_present),
                                                       replace=False)]))
            for _ in range(1000)]
    p_nrf2 = float((np.sum(np.array(null) >= obs) + 1) / 1001) if obs is not None else None
    top = sorted(genes_all, key=lambda g: (-frac[g], -gene_tot[g]))
    top = [g for g in top if gene_tot[g] >= 10][:30]
    q3 = {"genes_tested": len(genes_all), "nrf2_genes_on_array": nrf2_present,
          "nrf2_mean_frac_sig": obs, "random_null_mean": float(np.mean(null)),
          "random_null_p95": float(np.percentile(null, 95)), "p_empirical": p_nrf2,
          "enriched": bool(p_nrf2 is not None and p_nrf2 < 0.05),
          "top30_genes": [{"gene": g, "frac_sig": round(frac[g], 3), "probes": gene_tot[g]}
                          for g in top],
          "per_nrf2_gene": {g: {"frac_sig": round(frac[g], 3), "probes": gene_tot[g]}
                            for g in nrf2_present}}
    run.log("Q3_nrf2", observed=round(obs, 4) if obs else None,
            null_mean=round(float(np.mean(null)), 4), p=p_nrf2, n_genes=len(nrf2_present))

    # ---- Q4 is it a smoking signature? ----
    wt = (~keap1) & luad
    ever = np.array([("Smoker" in s) and ("Non-Smoker" not in s) for s in smoke])
    never = np.array(["Non-Smoker" in s for s in smoke])
    q4 = {"n_wt_ever": int((wt & ever).sum()), "n_wt_never": int((wt & never).sum())}
    if (wt & never).sum() >= 20:
        q_smoke = diff_probes(wt & ever, wt & never, "signature_smoking_in_KEAP1wt")
        sig_smoke = q_smoke < 0.05
        inter = int((sig_keap1 & sig_smoke).sum())
        orr, pf = fisher_exact([[inter, int((sig_keap1 & ~sig_smoke).sum())],
                                [int((~sig_keap1 & sig_smoke).sum()),
                                 int((~sig_keap1 & ~sig_smoke).sum())]])
        jac = inter / max(int((sig_keap1 | sig_smoke).sum()), 1)
        q4.update({
            "smoking_signature_size": int(sig_smoke.sum()),
            "keap1_signature_size": int(sig_keap1.sum()),
            "overlap": inter,
            "frac_of_keap1_shared": round(inter / max(int(sig_keap1.sum()), 1), 4),
            "jaccard": round(jac, 4), "odds_ratio": round(float(orr), 3), "p": float(pf),
            "verdict": None})
        q4["verdict"] = (
            "the KEAP1 signature is largely a smoking signature"
            if q4["frac_of_keap1_shared"] > 0.5 else
            "the KEAP1 signature is mostly distinct from the smoking signature")
        run.log("Q4_overlap", **{k: v for k, v in q4.items() if k != "verdict"})
        run.log("Q4_verdict", verdict=q4["verdict"])
    else:
        q4["skipped"] = f"only {q4['n_wt_never']} KEAP1-wildtype never-smokers"

    results = {
        "cohort": {"luad_tumours": int(luad.sum()),
                   "keap1_mutant": int((keap1 & luad).sum()),
                   "probes": len(pnames), "annotated": len(ann)},
        "signature": {"n_significant": int(sig_keap1.sum()),
                      "frac_significant": float(sig_keap1.mean()),
                      "R08_sampled_estimate": 0.2935},
        "PC_sex_classifier": pc_sex,
        "Q1_keap1_classifier": q1,
        "Q2_genomic_context": q2,
        "Q3_nrf2_enrichment": q3,
        "Q4_smoking_overlap": q4,
        "_decision": {
            "phenotype_strength_auroc": q1["auroc"],
            "nrf2_enriched": q3["enriched"],
            "distinct_from_smoking": (q4.get("frac_of_keap1_shared", 1) <= 0.5
                                      if "frac_of_keap1_shared" in q4 else None),
            "prespecified_in_config_sha256": cfg_hash},
        "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir},
    }
    np.savetxt(os.path.join(run.dir, "oof_keap1_classifier.csv"),
               np.column_stack([keap1[luad].astype(int), oof_keap1]), delimiter=",",
               header="keap1_mutant,oof_probability", comments="")
    run.write("results.json", results)
    run.write("keap1_signature_qvalues.csv",
              "probe,q_value,significant\n" + "\n".join(
                  f"{p},{q:.6g},{int(s)}" for p, q, s in
                  zip(pnames, q_keap1, sig_keap1) if s))
    run.log("decision", **results["_decision"])
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
