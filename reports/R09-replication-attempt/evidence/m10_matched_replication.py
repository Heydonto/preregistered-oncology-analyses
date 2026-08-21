#!/usr/bin/env python3
"""R10 — setting-matched replication of R06, in the cohorts a systematic screen selected.

Sequence of events, all recorded:
  R09  TCGA cannot adjudicate R06's KEAP1/SMARCA4 finding. 385 events, and both positive
       controls still failed, in the full cohort and in the pre-declared late-stage subgroup.
       TCGA lung is resected disease: median OS 50.3 months against 28.8 in R06's cohort.
  R10 first attempt  luad_mskcc_2023_met_organotropism was rejected by a pre-registered setting
       gate: median OS 93.8 months, so despite the name it is not an advanced-disease cohort.
       That run is retained.
  screen_replication_cohorts.py  screened all 41 lung studies in cBioPortal against five
       requirements decided before any hazard ratio was computed. Two cohorts qualified.

This report analyses those two:
  bm_nsclc_mskcc_2023  NSCLC brain metastasis, 209 patients, 114 events, median OS 29.7 months
  lung_msk_pdx         PDX-derived series,      171 patients, 117 events, median OS 24.7 months

Model is R06's, unchanged: per-gene unadjusted Cox on OS, BH across the same five genes.
The pooled analysis stratifies by cohort, which keeps every comparison within-cohort (so the
two very different selections are never contrasted with each other) while recovering power.

Power computed BEFORE running, so the control rule cannot be rationalised afterwards. At 231
pooled events: STK11 (R06 HR 1.63) gives z~2.6, so it should fire if the effect is real; TP53
(R06 HR 1.25) gives z~1.7 and may not. The control rule is therefore a disjunction, and STK11
is the load-bearing control. If neither fires, this is a power failure of the available public
data, not evidence against R06.

Both cohorts are MSK, as R06's was. Patients are disjoint and gated as such, but this is not
institutional independence, and the report says so rather than implying more than it has.
"""
import json, os, subprocess, sys
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_kit import Run, bh

COHORTS = ["bm_nsclc_mskcc_2023", "lung_msk_pdx"]
R06_STUDY = "nsclc_ctdx_msk_2022"
GENES = {9817: "KEAP1", 6597: "SMARCA4", 8085: "KMT2D", 6794: "STK11", 7157: "TP53"}
CHROM = ["KEAP1", "SMARCA4", "KMT2D"]
CONTROLS = ["TP53", "STK11"]
SILENT = {"Silent", "Intron", "3'UTR", "5'UTR", "5'Flank", "IGR", "RNA"}
SEED = 20260807
R06 = {"KEAP1": {"hr": 2.07, "ci": [1.63, 2.62], "q": 0.0001},
       "SMARCA4": {"hr": 1.73, "ci": [1.27, 2.35], "q": 0.0009},
       "KMT2D": {"hr": None, "ci": None, "q": None},
       "STK11": {"hr": 1.63, "ci": None, "q": 0.0004},
       "TP53": {"hr": 1.25, "ci": None, "q": 0.0074}}
EXPECT = {"pooled_events": 231, "cohorts": 2,
          "prerun_power": {"STK11_z": 2.6, "TP53_z": 1.7, "KEAP1_z": 3.8}}


def curl(url, body=None):
    cmd = ["curl", "-s", "--max-time", "240", "--retry", "3", url]
    if body is not None:
        cmd += ["-X", "POST", "-H", "Content-Type: application/json", "-d", body]
    return json.loads(subprocess.run(cmd, capture_output=True, timeout=300).stdout)


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def load(study, r06_pat):
    cl = curl(f"https://www.cbioportal.org/api/studies/{study}/clinical-data"
              f"?clinicalDataType=PATIENT&projection=SUMMARY&pageSize=200000")
    om = {x["patientId"]: x["value"] for x in cl if x["clinicalAttributeId"] == "OS_MONTHS"}
    st = {x["patientId"]: x["value"] for x in cl if x["clinicalAttributeId"] == "OS_STATUS"}
    usable = [p for p in om if st.get(p) and (num(om[p]) or 0) > 0]
    disj = [p for p in usable if p not in r06_pat]
    d = curl(f"https://www.cbioportal.org/api/molecular-profiles/{study}_mutations/"
             f"mutations/fetch?projection=DETAILED",
             json.dumps({"sampleListId": f"{study}_all", "entrezGeneIds": list(GENES)}))
    per = {g: set() for g in GENES.values()}
    for x in d:
        g = GENES.get(x.get("entrezGeneId"))
        if g and x.get("mutationType") not in SILENT and x.get("patientId"):
            per[g].add(x["patientId"])
    df = pd.DataFrame([{"case": p, "cohort": study, "T": num(om[p]),
                        "E": int(st[p].startswith("1")),
                        **{g: int(p in per[g]) for g in GENES.values()}} for p in disj])
    return df, len(usable) - len(disj), len(d)


def per_gene(frame, run, label, strata=None):
    res, ps, keys = {}, [], []
    for g in list(CHROM) + CONTROLS:
        if frame[g].sum() < 10:
            res[g] = {"skipped": "fewer than 10 mutated", "n_mutated": int(frame[g].sum())}
            continue
        cols = ["T", "E", g] + ([strata] if strata else [])
        cph = CoxPHFitter().fit(frame[cols], "T", "E",
                                strata=[strata] if strata else None)
        m = {}
        for v in (1, 0):
            sel = frame[g] == v
            k = KaplanMeierFitter().fit(frame["T"][sel], frame["E"][sel])
            m[v] = float(k.median_survival_time_)
        lr = logrank_test(frame["T"][frame[g] == 1], frame["T"][frame[g] == 0],
                          frame["E"][frame[g] == 1], frame["E"][frame[g] == 0])
        res[g] = {"n_mutated": int(frame[g].sum()), "n_total": int(len(frame)),
                  "events": int(frame["E"].sum()),
                  "hr": float(np.exp(cph.params_[g])),
                  "ci95": [float(np.exp(cph.confidence_intervals_.loc[g].iloc[0])),
                           float(np.exp(cph.confidence_intervals_.loc[g].iloc[1]))],
                  "p_cox": float(cph.summary.loc[g, "p"]),
                  "p_logrank": float(lr.p_value),
                  "median_months_mut": m[1], "median_months_wt": m[0]}
        ps.append(res[g]["p_cox"]); keys.append(g)
    for g, q in zip(keys, bh(ps)):
        res[g]["q_bh"] = float(q)
        res[g]["significant_q05"] = bool(q < 0.05)
    for g in keys:
        run.log(label, gene=g, hr=round(res[g]["hr"], 3),
                ci=[round(x, 2) for x in res[g]["ci95"]],
                q=round(res[g]["q_bh"], 5), n_mut=res[g]["n_mutated"])
    return res


def main():
    run = Run("R10matchedrepl")
    cfg = {
        "report": "R10",
        "question": "does R06's KEAP1/SMARCA4 association replicate in a setting-matched, "
                    "patient-disjoint advanced-disease cohort?",
        "history": {
            "R09": "TCGA cannot adjudicate - 385 events, both controls failed in the full "
                   "cohort and the late-stage subgroup; median OS 50.3mo (resected disease)",
            "R10_first_attempt": "luad_mskcc_2023_met_organotropism rejected by a pre-registered "
                                 "setting gate (median OS 93.8mo); that run is retained",
            "screen": "all 41 cBioPortal lung studies screened on 5 requirements fixed before "
                      "any hazard ratio was computed; 2 qualified"},
        "cohorts": COHORTS,
        "model": "R06's exactly - per-gene UNADJUSTED Cox on OS, BH across the same 5 genes",
        "pooled_model": "Cox stratified by cohort, so all comparisons stay within-cohort while "
                        "recovering power; the two selections are never contrasted",
        "positive_controls": CONTROLS,
        "prerun_power_note": "computed BEFORE running: at 231 pooled events STK11 (HR 1.63) "
                             "gives z~2.6 and should fire; TP53 (HR 1.25) gives z~1.7 and may "
                             "not. STK11 is the load-bearing control",
        "control_rule": "TP53 or STK11 at q<0.05 in the POOLED analysis (where power is "
                        "greatest). If neither fires this is a power failure of public data, "
                        "not evidence against R06",
        "replication_rule": "a gene replicates only if q<0.05 AND same direction as R06",
        "overlap_rule": f"patients present in {R06_STUDY} removed before any outcome is read",
        "institutional_caveat": "both cohorts are MSK, as R06's was. Patients are disjoint; "
                               "this is NOT institutional independence",
        "r06_reference_values": R06, "seed": SEED, "expectations": EXPECT,
    }
    cfg_hash = run.start(cfg, [])

    r06_pat = {x["patientId"] for x in
               curl(f"https://www.cbioportal.org/api/studies/{R06_STUDY}/patients?pageSize=100000")}
    frames, meta = [], {}
    for c in COHORTS:
        df, removed, nrec = load(c, r06_pat)
        km = KaplanMeierFitter().fit(df["T"], df["E"])
        med = float(km.median_survival_time_)
        meta[c] = {"n": len(df), "events": int(df["E"].sum()), "median_os_months": med,
                   "overlap_removed": removed, "mutation_records": nrec,
                   "mutated": {g: int(df[g].sum()) for g in GENES.values()}}
        run.log("cohort", study=c, **{k: v for k, v in meta[c].items() if k != "mutated"})
        run.gate(f"G_setting_{c}", "median OS nearer R06's 28.8mo than TCGA's 50.3mo",
                 round(med, 1), abs(med - 28.8) < abs(med - 50.3),
                 "the gate that rejected the first R10 candidate at 93.8mo")
        frames.append(df)

    pool = pd.concat(frames, ignore_index=True)
    run.gate("G0_pooled_events", EXPECT["pooled_events"], int(pool["E"].sum()),
             int(pool["E"].sum()) == EXPECT["pooled_events"],
             "R06's failed replication cohort had 45")
    run.gate("G1_cohorts", EXPECT["cohorts"], pool["cohort"].nunique(),
             pool["cohort"].nunique() == EXPECT["cohorts"])
    run.gate("G2_no_overlap_remains", 0, int(pool["case"].isin(r06_pat).sum()),
             int(pool["case"].isin(r06_pat).sum()) == 0,
             "no R06 patient survives into the pooled frame")

    per_cohort = {c: per_gene(f.reset_index(drop=True), run, f"cohort_{c}")
                  for c, f in zip(COHORTS, frames)}
    pooled = per_gene(pool, run, "pooled", strata="cohort")

    ctrl_ok = any(pooled.get(c, {}).get("significant_q05") for c in CONTROLS)
    run.gate("G3_positive_controls", "TP53 or STK11 at q<0.05 in the pooled analysis",
             {c: {"hr": round(pooled.get(c, {}).get("hr", 0), 3),
                  "q": round(pooled.get(c, {}).get("q_bh", 1), 5)} for c in CONTROLS},
             ctrl_ok, "STK11 is load-bearing per the pre-run power note")

    verdict = {}
    for g in CHROM:
        p = pooled.get(g, {})
        same = (p["hr"] > 1) == (R06[g]["hr"] > 1) if R06[g]["hr"] and "hr" in p else None
        verdict[g] = {
            "r06_hr": R06[g]["hr"], "r06_ci": R06[g]["ci"],
            "pooled_hr": p.get("hr"), "pooled_ci95": p.get("ci95"),
            "pooled_q": p.get("q_bh"), "same_direction": same,
            "replicates": bool(p.get("significant_q05") and same),
            "hr_inside_r06_ci": (bool(R06[g]["ci"][0] <= p["hr"] <= R06[g]["ci"][1])
                                 if R06[g].get("ci") and p.get("hr") else None),
            "per_cohort_hr": {c: per_cohort[c].get(g, {}).get("hr") for c in COHORTS}}

    rep = [g for g in CHROM if verdict[g]["replicates"]]
    results = {
        "cohort_meta": meta,
        "pooled": {"n": len(pool), "events": int(pool["E"].sum())},
        "positive_controls_pooled": {c: {"hr": pooled.get(c, {}).get("hr"),
                                         "q": pooled.get(c, {}).get("q_bh")}
                                     for c in CONTROLS},
        "controls_recovered": ctrl_ok,
        "per_cohort": per_cohort,
        "pooled_per_gene": pooled,
        "replication_verdict": verdict,
        "_decision": {
            "rule": "replicates only if q<0.05 and same direction as R06",
            "replicating": rep,
            "interpretation": (f"replicates in setting-matched, patient-disjoint cohorts: "
                               f"{', '.join(rep)}" if rep else
                               "no replication; read against the control result and the "
                               "pre-run power note before treating this as a refutation"),
            "caveat": "both cohorts are MSK; disjoint patients but not institutionally "
                      "independent",
            "prespecified_in_config_sha256": cfg_hash},
        "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir},
    }
    run.write("results.json", results)
    run.write("survival_frame_r10.csv", pool.to_csv(index=False))
    run.log("decision", **{k: v for k, v in results["_decision"].items() if k != "rule"})
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
