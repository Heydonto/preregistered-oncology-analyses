#!/usr/bin/env python3
"""R23 -- execute the hashed pre-registration on GSE68379 lung lines + GDSC2.

Plan: PREREG_PLAN.yaml, SHA-256 7e0322ba804aef2f84ecf559d2e348f146c09fab1086cb2148235a377cde2a66,
hashed before any IC50 was joined to any methylation value. This executor implements that plan and
nothing else. Where it departs, the departure is logged as a gate note.

Weaker than the Paper 2 protocol in one respect, stated rather than glossed: the plan author and
the executor are the same person, so this is the "hash before label" form and not the structural
isolation form.

Run:  python3 run_analysis.py
"""
import collections
import csv
import gzip
import hashlib
import io
import json
import os
import subprocess
import sys

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, "/Users/rezanehzati/Projects/quantara/paper-draft/02-revised-protocol-open-datasets")
from audit_kit import Run

STAGE = "/Users/rezanehzati/quantara-staging/drugresist/"
BETAS = "gs://heydonto-quantara-lungcdx/data-request-2026-08/labels/lung/GSE68379/processed/lung_betas.npz"
META = "gs://heydonto-quantara-lungcdx/data-request-2026-08/labels/lung/GSE68379/processed/lung_meta.csv"
PLAN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PREREG_PLAN.yaml")
PLAN_SHA = "7e0322ba804aef2f84ecf559d2e348f146c09fab1086cb2148235a377cde2a66"

SEED = 20260820
# DEPARTURE FROM THE HASHED PLAN, declared not silent. The plan specified 200 permutations.
# With a p-floor of 1/201 = 0.00498, Holm across 11 standard-of-care drugs cannot go below
# 0.0547 and across 32 chromatin drugs cannot go below 0.159 -- significance was arithmetically
# unreachable no matter what the data showed, and the first run duly reported 0 predictable
# drugs while Vorinostat sat at rho 0.72. Raised to 2000 (floor 0.0005; x32 = 0.016). The change
# can only make the test MORE able to reject, never less, and it is logged as gate G6.
N_PERM = 2000
N_VAR_PROBES = 5000
CV = 5
MIN_LINES = 100

SOC = ["Docetaxel", "Cisplatin", "Gemcitabine", "Paclitaxel", "Afatinib", "Gefitinib",
       "Erlotinib", "Osimertinib", "Olaparib", "Navitoclax", "Camptothecin"]
THESIS_NEG, THESIS_POS = -0.20, 0.20


def holm(pvals):
    idx = np.argsort(pvals)
    out = np.empty(len(pvals))
    run_max = 0.0
    for r, i in enumerate(idx):
        v = (len(pvals) - r) * pvals[i]
        run_max = max(run_max, v)
        out[i] = min(run_max, 1.0)
    return out


def fetch(uri, dest):
    if not os.path.exists(dest):
        subprocess.run(["gcloud", "storage", "cp", uri, dest], check=True, capture_output=True)
    return dest


def main():
    run = Run("R23drugresist")

    # ---- the plan must be the hashed one, checked before anything is computed ----
    got = hashlib.sha256(open(PLAN, "rb").read()).hexdigest()
    cfg = run.start(
        {"executes": "PREREG_PLAN.yaml", "declared_sha256": PLAN_SHA, "observed_sha256": got,
         "seed": SEED, "n_perm": N_PERM, "n_var_probes": N_VAR_PROBES, "cv": CV,
         "soc_agents": SOC, "thesis_thresholds": {"supported": THESIS_NEG,
                                                  "contradicted": THESIS_POS},
         "isolation": "hash-before-label only; plan author and executor are the same person",
         "status": "HELD FOR IP - not for publication pending patent"},
        [PLAN])
    run.gate("G0_plan_unmodified", PLAN_SHA, got, got == PLAN_SHA,
             "the plan executed must be the plan that was hashed before outcomes were joined")

    # ---- load ----
    d = np.load(fetch(BETAS, STAGE + "lung_betas.npz"), allow_pickle=True)
    betas = d["betas"].astype(np.float32)
    probes = np.array([str(x) for x in d["probes"]], object)
    gsms = [str(x) for x in d["gsm"]]
    meta = {r["gsm"]: r for r in csv.DictReader(open(fetch(META, STAGE + "lung_meta.csv")))}
    run.log("betas", shape=list(betas.shape), arrays=len(gsms), labelled=len(meta))

    ic = collections.defaultdict(dict)
    path_of, n_lines = {}, collections.Counter()
    with open(STAGE + "GDSC2_fitted.csv") as f:
        for r in csv.DictReader(f):
            ic[r["COSMIC_ID"].strip()][r["DRUG_NAME"]] = float(r["LN_IC50"])
            path_of[r["DRUG_NAME"]] = r["PATHWAY_NAME"]

    keep = [g for g in gsms if g in meta and meta[g]["cosmic_id"] in ic]
    run.gate("G1_overlap", "187 +/- 5 lines with both methylation and GDSC2", len(keep),
             abs(len(keep) - 187) <= 5)
    col = {g: i for i, g in enumerate(gsms)}
    X_all = betas[:, [col[g] for g in keep]].T          # lines x probes
    cos = [meta[g]["cosmic_id"] for g in keep]
    hist = [meta[g]["histology"] for g in keep]

    # ---- E4/E5 probe filtering ----
    frac_nan = np.isnan(X_all).mean(0)
    ok = frac_nan <= 0.05
    sexlike = np.array([False] * len(probes))           # 450k sex probes need the manifest;
    # the manifest is not staged locally, so E5 is applied via the chrX/chrY probe-name prefix
    # convention only where available. Recorded as a departure.
    keepp = ok & ~sexlike
    X = X_all[:, keepp]
    P = probes[keepp]
    X = np.where(np.isnan(X), np.nanmedian(X, axis=0, keepdims=True), X)
    run.gate("G2_probes", "300,000-460,000 probes after filtering", int(X.shape[1]),
             300000 <= X.shape[1] <= 460000,
             "E5 (sex chromosomes) could not be applied: the 450k manifest is not staged "
             "locally, so no probe was dropped on that criterion. Departure from the plan.")

    # ---- most-variable probe panel, chosen from PREDICTORS only, never from outcomes ----
    v = X.var(0)
    top = np.argsort(v)[::-1][:N_VAR_PROBES]
    Xv = StandardScaler().fit_transform(X[:, top])

    # ---- PC1: histology separability ----
    y_h = np.array([1 if "small_cell" in h else 0 for h in hist])
    skf = StratifiedKFold(CV, shuffle=True, random_state=SEED)
    oof = np.empty(len(y_h), int)
    for tr, te in skf.split(Xv, y_h):
        m = LogisticRegression(max_iter=5000).fit(Xv[tr], y_h[tr])
        oof[te] = m.predict(Xv[te])
    pc1 = float(balanced_accuracy_score(y_h, oof))
    run.gate("PC1_histology", "SCLC vs non-SCLC balanced accuracy > 0.80", round(pc1, 4),
             pc1 > 0.80, "plan: if this fails the beta matrix is not carrying biology; halt")

    # ---- H1 / H2: does methylation predict LN_IC50? ----
    rng = np.random.default_rng(SEED)

    def predict_drug(drug):
        """Out-of-fold ridge with a permutation null.

        Uses the kernel (dual) form: with n~150 lines and p=5000 probes, ridge predictions are
        H @ y where H depends only on X and the fold split. H is built once per drug, so each of
        the 2000 permutations costs one matrix-vector product instead of a refit. Without this the
        permutation count needed for a valid Holm correction would have taken hours.
        """
        idx = [i for i, c in enumerate(cos) if drug in ic[c]]
        if len(idx) < MIN_LINES:
            return None
        y = np.array([ic[cos[i]][drug] for i in idx])
        Xd = Xv[idx]
        K = Xd @ Xd.T
        n = len(y)
        H = np.zeros((n, n))
        for tr, te in KFold(CV, shuffle=True, random_state=SEED).split(Xd):
            A = K[np.ix_(tr, tr)] + 1000.0 * np.eye(len(tr))
            H[np.ix_(te, tr)] = K[np.ix_(te, tr)] @ np.linalg.inv(A)

        def rho_of(yy):
            yc = yy - yy.mean()
            return spearmanr(yy, H @ yc + yy.mean()).statistic

        obs = float(rho_of(y))
        null = np.array([rho_of(rng.permutation(y)) for _ in range(N_PERM)])
        p = (1 + int((null >= obs).sum())) / (N_PERM + 1)
        return {"drug": drug, "n": n, "rho": obs, "perm_p": p,
                "null_mean": float(null.mean()),
                "pathway": path_of.get(drug, "")}

    soc_res = [r for r in (predict_drug(d) for d in SOC) if r]
    chrom_drugs = sorted({d for d in path_of if "Chromatin" in path_of[d]})
    chrom_avail = [d for d in chrom_drugs
                   if sum(1 for c in cos if d in ic[c]) >= MIN_LINES]
    chrom_res = [r for r in (predict_drug(d) for d in chrom_avail) if r]
    for grp, res in (("H1_soc", soc_res), ("H2_chromatin", chrom_res)):
        if res:
            q = holm([r["perm_p"] for r in res])
            for r, qq in zip(res, q):
                r["holm_p"] = float(qq)
            n_sig = sum(1 for r in res if r["holm_p"] < 0.05)
            run.log(grp, n_drugs=len(res), n_predictable=n_sig,
                    best=max(res, key=lambda r: r["rho"])["drug"],
                    best_rho=round(max(r["rho"] for r in res), 4))

    run.gate("G6_multiplicity_floor",
             "the permutation floor must permit Holm significance across the larger family",
             {"n_perm": N_PERM, "floor": round(1 / (N_PERM + 1), 5),
              "n_chromatin": len(chrom_res),
              "min_reachable_holm": round(len(chrom_res) / (N_PERM + 1), 4)},
             len(chrom_res) / (N_PERM + 1) < 0.05,
             "the plan's 200 permutations made significance unreachable; this gate exists so the "
             "same mistake fails loudly next time")
    run.gate("G3_nulls_centred", "every permutation null within 0.05 of zero",
             round(float(np.max([abs(r["null_mean"]) for r in soc_res + chrom_res])), 4),
             all(abs(r["null_mean"]) < 0.05 for r in soc_res + chrom_res))
    run.gate("G4_min_lines", f"no retained drug below {MIN_LINES} lines",
             int(min(r["n"] for r in soc_res + chrom_res)),
             min(r["n"] for r in soc_res + chrom_res) >= MIN_LINES)

    # ---- H3 PRIMARY: are resistant lines differentially sensitive to chromatin drugs? ----
    soc_set = {r["drug"] for r in soc_res}
    chrom_set = {r["drug"] for r in chrom_res}
    run.gate("G5_disjoint_drug_sets", "resistance and chromatin indices use disjoint drugs",
             sorted(soc_set & chrom_set), not (soc_set & chrom_set))

    def index_over(drugs):
        z = {}
        for d in drugs:
            vals = {c: ic[c][d] for c in cos if d in ic[c]}
            if len(vals) < MIN_LINES:
                continue
            a = np.array(list(vals.values()))
            mu, sd = a.mean(), a.std(ddof=1)
            for c, val in vals.items():
                z.setdefault(c, []).append((val - mu) / sd)
        return {c: float(np.mean(v)) for c, v in z.items() if len(v) >= 3}

    ri, ci = index_over(soc_set), index_over(chrom_set)
    common = [c for c in cos if c in ri and c in ci]
    a = np.array([ri[c] for c in common])
    b = np.array([ci[c] for c in common])
    rho = float(spearmanr(a, b).statistic)
    nullv = np.array([spearmanr(a, rng.permutation(b)).statistic for _ in range(N_PERM)])
    p3 = (1 + int((np.abs(nullv) >= abs(rho)).sum())) / (N_PERM + 1)

    # ---- POST-HOC CONTROL, not in the hashed plan, labelled as such -------------------
    # Both indices are mean z-scored LN_IC50, so both load on the general drug-sensitivity axis:
    # a line resistant to nearly everything scores high on each. A positive correlation between
    # them may therefore say nothing about chromatin drugs specifically. This partials out a
    # general index built from every GDSC2 drug that is in NEITHER set. The plan did not declare
    # this and the result is reported as an observation, not as a pre-registered test.
    other = [d for d in path_of
             if d not in soc_set and d not in chrom_set
             and sum(1 for c in cos if d in ic[c]) >= MIN_LINES]
    gi = index_over(set(other))
    common3 = [c for c in common if c in gi]
    g = np.array([gi[c] for c in common3])
    a3 = np.array([ri[c] for c in common3])
    b3 = np.array([ci[c] for c in common3])

    def prank(x):
        from scipy.stats import rankdata
        r = rankdata(x).astype(float)
        return (r - r.mean()) / r.std(ddof=1)

    ra, rb, rg = prank(a3), prank(b3), prank(g)
    # residualise both on the general axis, then correlate
    res_a = ra - rg * (ra @ rg) / (rg @ rg)
    res_b = rb - rg * (rb @ rg) / (rg @ rg)
    partial = float(spearmanr(res_a, res_b).statistic)
    nullp = np.array([spearmanr(res_a, rng.permutation(res_b)).statistic
                      for _ in range(N_PERM)])
    p_partial = (1 + int((np.abs(nullp) >= abs(partial)).sum())) / (N_PERM + 1)
    run.log("H3_POSTHOC_general_sensitivity_control", n_other_drugs=len(other),
            n_lines=len(common3), raw_rho=round(rho, 4), partial_rho=round(partial, 4),
            perm_p=round(p_partial, 5),
            note="POST-HOC, not pre-declared; partials out an index built from drugs in neither set")

    if rho <= THESIS_NEG and p3 < 0.05:
        verdict = "THESIS_SUPPORTED_resistant_lines_more_sensitive_to_chromatin_drugs"
    elif rho >= THESIS_POS and p3 < 0.05:
        verdict = "THESIS_CONTRADICTED_resistant_lines_are_cross_resistant"
    else:
        verdict = "NULL_no_relationship_at_the_pre_declared_effect_size"
    run.log("H3_PRIMARY", n_lines=len(common), rho=round(rho, 4), perm_p=round(p3, 5),
            null_mean=round(float(nullv.mean()), 4), verdict=verdict)

    out = {"status": "HELD FOR IP - not for publication pending patent",
           "plan_sha256": got,
           "VERDICT": {"H3_primary": verdict, "rho": rho, "perm_p": p3,
                       "n_lines": len(common),
                       "rules_pre_declared_in": "PREREG_PLAN.yaml hashed " + PLAN_SHA[:16]},
           "cohort": {"lines_with_both": len(keep), "probes_after_filter": int(X.shape[1]),
                      "panel": N_VAR_PROBES,
                      "histology": dict(collections.Counter(hist).most_common())},
           "PC1_histology_balanced_accuracy": pc1,
           "H1_standard_of_care": sorted(soc_res, key=lambda r: -r["rho"]),
           "H2_chromatin": sorted(chrom_res, key=lambda r: -r["rho"]),
           "H3_indices": {"n_soc_drugs": len(soc_set), "n_chromatin_drugs": len(chrom_set),
                          "disjoint": True},
           "H3_POSTHOC_general_sensitivity_control": {
               "STATUS": "POST-HOC OBSERVATION, NOT PRE-DECLARED",
               "n_control_drugs": len(other), "n_lines": len(common3),
               "raw_rho": rho, "partial_rho": partial, "perm_p": p_partial,
               "reading": ("the raw association is largely the general drug-sensitivity axis"
                           if abs(partial) < 0.2 else
                           "the association survives adjustment for general sensitivity")},
           "_provenance": {"config_sha256": cfg, "run_dir": run.dir}}
    np.save(os.path.join(run.dir, "resistance_index.npy"), a)
    np.save(os.path.join(run.dir, "chromatin_index.npy"), b)
    run.write("results.json", out)
    run.finalize()
    print("\nH3 PRIMARY:", verdict)
    print(f"  rho={rho:+.4f} p={p3:.5f} on {len(common)} lines")
    print("RUN DIR:", run.dir)


if __name__ == "__main__":
    main()
