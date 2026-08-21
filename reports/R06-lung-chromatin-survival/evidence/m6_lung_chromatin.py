#!/usr/bin/env python3
"""R06 — chromatin-regulator mutations and survival in NSCLC, two independent cohorts.

This is the reframed lung arm: the response letter proposed abandoning the acquired-resistance
claim (which public data cannot power) in favour of chromatin dysregulation. Here it is tested
directly on mutation and survival data, with no imaging and no methylation required.

  luad_mskcc_2020     604 patients, OS and RFS complete
  nsclc_ctdx_msk_2022 1,127 patients, OS complete

Genes: KMT2D, SMARCA4, KEAP1 (chromatin/oxidative-stress regulators).
POSITIVE CONTROLS: TP53 and STK11, both established in NSCLC. If neither is recovered the
pipeline is untrustworthy and nothing else is reported.

Mutation data comes from the cBioPortal REST API, not the repository files: the two studies'
data_mutations.txt are 132-byte git-LFS pointer stubs whose LFS objects are unretrievable
upstream (media endpoint 404, S3 403). This is recorded as gate G1.
"""
import csv, json, os, subprocess, sys
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_kit import Run, bh

L = "/Users/rezanehzati/quantara-staging/staged/labels/lung/"
STUDIES = {"luad_mskcc_2020": 604, "nsclc_ctdx_msk_2022": 1127}
GENES = {8085: "KMT2D", 6597: "SMARCA4", 9817: "KEAP1", 6794: "STK11", 7157: "TP53"}
CHROM = ["KMT2D", "SMARCA4", "KEAP1"]
CONTROLS = ["TP53", "STK11"]
SILENT = {"Silent", "Intron", "3'UTR", "5'UTR", "5'Flank", "IGR", "RNA"}
SEED = 20260807


def api_mutations(study):
    body = json.dumps({"sampleListId": f"{study}_all", "entrezGeneIds": list(GENES)})
    out = subprocess.run(
        ["curl", "-s", "--retry", "3", "-X", "POST",
         f"https://www.cbioportal.org/api/molecular-profiles/{study}_mutations/"
         f"mutations/fetch?projection=DETAILED",
         "-H", "Content-Type: application/json", "-d", body],
        capture_output=True, timeout=300).stdout
    d = json.loads(out)
    per = {g: set() for g in GENES.values()}
    for x in d:
        g = GENES.get(x.get("entrezGeneId"))
        if g and x.get("mutationType") not in SILENT:
            pid = x.get("patientId") or x.get("sampleId")
            if pid:
                per[g].add(pid)
    return per, len(d)


def rd(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def cl(v):
    v = (v or "").strip()
    return "" if v in ("NA", "[Not Available]", "[Not Applicable]", "[Unknown]") else v


def main():
    run = Run("R06lungchrom")
    cfg = {
        "report": "R06",
        "question": "do chromatin-regulator mutations associate with overall survival in NSCLC?",
        "motivation": "the reframed lung arm proposed in RESPONSE_TO_COLLEAGUE; needs no "
                      "methylation download and no imaging",
        "cohorts": STUDIES,
        "genes_of_interest": CHROM,
        "positive_controls": CONTROLS,
        "control_rule": "if neither TP53 nor STK11 is recovered in the larger cohort, halt",
        "primary": "per-gene Cox on overall survival, each cohort separately",
        "key_test": "replication - does a gene reach significance in BOTH cohorts?",
        "multiplicity": "Benjamini-Hochberg across the 5 genes within each cohort",
        "data_source_note": "mutations via cBioPortal REST API; the repository data_mutations.txt "
                            "files are 132-byte LFS pointer stubs with unretrievable objects",
        "seed": SEED,
    }
    cfg_hash = run.start(cfg, [L + s + "/data_clinical_patient.txt" for s in STUDIES])

    allres, cohort_meta = {}, {}
    for study, expect_n in STUDIES.items():
        pat = {r["PATIENT_ID"]: r for r in rd(L + study + "/data_clinical_patient.txt")}
        run.gate(f"G0_{study}_patients", expect_n, len(pat), len(pat) == expect_n)

        per, nrec = api_mutations(study)
        run.gate(f"G1_{study}_api_mutations", "> 100 records from REST API", nrec, nrec > 100,
                 "repository files are LFS stubs; API is the working route")

        rows = []
        for pid, p in pat.items():
            t, s = cl(p.get("OS_MONTHS", "")), cl(p.get("OS_STATUS", ""))
            if not t or not s:
                continue
            try:
                tt = float(t)
            except ValueError:
                continue
            if tt <= 0:
                continue
            rows.append({"T": tt, "E": 1 if s.startswith("1") else 0,
                         **{g: int(pid in per[g]) for g in GENES.values()}})
        df = pd.DataFrame(rows)
        run.gate(f"G2_{study}_survival", "> 400 usable", len(df), len(df) > 400,
                 "OS time > 0 and status present")
        cohort_meta[study] = {"n": int(len(df)), "events": int(df["E"].sum()),
                              "mutated": {g: int(df[g].sum()) for g in GENES.values()},
                              "api_records": nrec}
        run.log("cohort", study=study, **cohort_meta[study])

        res, ps, keys = {}, [], []
        for g in list(CHROM) + CONTROLS:
            if df[g].sum() < 10:
                res[g] = {"skipped": "fewer than 10 mutated patients",
                          "n_mutated": int(df[g].sum())}
                continue
            cph = CoxPHFitter().fit(df[["T", "E", g]], "T", "E")
            km = KaplanMeierFitter()
            med = {}
            for v in (1, 0):
                m = df[g] == v
                km.fit(df["T"][m], df["E"][m])
                med[v] = float(km.median_survival_time_)
            lr = logrank_test(df["T"][df[g] == 1], df["T"][df[g] == 0],
                              df["E"][df[g] == 1], df["E"][df[g] == 0])
            res[g] = {"n_mutated": int(df[g].sum()), "n_total": int(len(df)),
                      "hr": float(np.exp(cph.params_[g])),
                      "ci95": [float(np.exp(cph.confidence_intervals_.iloc[0, 0])),
                               float(np.exp(cph.confidence_intervals_.iloc[0, 1]))],
                      "p_cox": float(cph.summary.loc[g, "p"]),
                      "p_logrank": float(lr.p_value),
                      "median_months_mut": med[1], "median_months_wt": med[0]}
            ps.append(res[g]["p_cox"]); keys.append(g)
        for g, q in zip(keys, bh(ps)):
            res[g]["q_bh"] = float(q)
            res[g]["significant_q05"] = bool(q < 0.05)
        allres[study] = res

    # positive-control gate on the larger cohort
    big = "nsclc_ctdx_msk_2022"
    ctrl_ok = any(allres[big].get(c, {}).get("significant_q05") for c in CONTROLS)
    run.gate("G3_positive_controls", "TP53 or STK11 significant in the larger cohort",
             {c: {"q": round(allres[big].get(c, {}).get("q_bh", 1), 4),
                  "hr": round(allres[big].get(c, {}).get("hr", 0), 3)} for c in CONTROLS},
             ctrl_ok, "established NSCLC biology")

    # replication across cohorts
    repl = {}
    for g in CHROM:
        a, b = allres["luad_mskcc_2020"].get(g, {}), allres[big].get(g, {})
        both = bool(a.get("significant_q05") and b.get("significant_q05"))
        same_dir = (a.get("hr", 1) > 1) == (b.get("hr", 1) > 1) if "hr" in a and "hr" in b else None
        repl[g] = {"luad_mskcc_2020": {k: a.get(k) for k in ("hr", "p_cox", "q_bh", "n_mutated")},
                   "nsclc_ctdx_msk_2022": {k: b.get(k) for k in ("hr", "p_cox", "q_bh", "n_mutated")},
                   "significant_in_both": both, "same_direction": same_dir}

    results = {"cohorts": cohort_meta, "per_cohort": allres, "replication": repl,
               "_decision": {
                   "chromatin_genes_replicating": [g for g in CHROM if repl[g]["significant_in_both"]],
                   "interpretation": "chromatin-regulator mutation(s) replicate as prognostic"
                   if any(repl[g]["significant_in_both"] for g in CHROM)
                   else "no chromatin-regulator mutation replicates as prognostic across both cohorts",
                   "prespecified_in_config_sha256": cfg_hash},
               "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir}}
    run.write("results.json", results)
    run.log("decision", **results["_decision"])
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
