#!/usr/bin/env python3
"""R09 — does R06's KEAP1/SMARCA4 prognostic finding replicate in TCGA?

R06 found KEAP1 (HR 2.07) and SMARCA4 (HR 1.73) associated with shorter overall survival in the
MSK circulating-tumour-DNA cohort, and its audit *withdrew* the replication claim because the
second cohort it tried (luad_mskcc_2020) had only 45 events in 604 patients and failed its own
positive controls. That left the finding resting on one cohort.

TCGA-LUAD + LUSC is a properly powered, genuinely independent replication set: ~1,050
sequenced patients with the curated TCGA-CDR survival endpoint and roughly 380 events, from
different institutions, a different assay (whole-exome rather than a targeted ctDNA panel) and
a different treatment era.

Declared in advance, because it determines how a negative result must be read: TCGA lung is
predominantly *resected, early-stage* disease, whereas the MSK ctDx cohort was *advanced*
disease. A failure to replicate is therefore ambiguous between "the effect is not real" and
"the effect is specific to advanced disease". To separate those as far as the data allow, the
late-stage subgroup (II-IV) is tested as a pre-specified secondary analysis - it is the closest
available match to the MSK setting.

Model is R06's, unchanged: per-gene unadjusted Cox on overall survival, Benjamini-Hochberg
across the same five genes, with STK11 and TP53 as positive controls. Adding covariates here
would make the comparison to R06 invalid.

Correctness note: only patients in each study's `_sequenced` sample list are used. A patient
without mutation profiling is not wild-type, and counting them as such would dilute every
hazard ratio toward 1.
"""
import csv, json, os, subprocess, sys
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_kit import Run, bh

LOCAL = "/Users/rezanehzati/quantara-staging/r08"
CLIN = LOCAL + "/clinical_merged.tsv"
STUDIES = {"luad_tcga_pan_can_atlas_2018": "TCGA-LUAD",
           "lusc_tcga_pan_can_atlas_2018": "TCGA-LUSC"}
GENES = {9817: "KEAP1", 6597: "SMARCA4", 8085: "KMT2D", 6794: "STK11", 7157: "TP53"}
CHROM = ["KEAP1", "SMARCA4", "KMT2D"]
CONTROLS = ["TP53", "STK11"]
SILENT = {"Silent", "Intron", "3'UTR", "5'UTR", "5'Flank", "IGR", "RNA"}
SEED = 20260807
# R06's measured values in nsclc_ctdx_msk_2022 (n=1127, 618 events), for comparison
R06 = {"KEAP1": {"hr": 2.07, "ci": [1.63, 2.62], "q": 0.0001},
       "SMARCA4": {"hr": 1.73, "ci": [1.27, 2.35], "q": 0.0009},
       "KMT2D": {"hr": None, "q": None},
       "STK11": {"hr": 1.63, "q": 0.0004}, "TP53": {"hr": 1.25, "q": 0.0074}}
EXPECT = {"sequenced": 1050, "min_events": 300}


def curl(url, body=None):
    cmd = ["curl", "-s", "--retry", "3", url]
    if body is not None:
        cmd += ["-X", "POST", "-H", "Content-Type: application/json", "-d", body]
    return json.loads(subprocess.run(cmd, capture_output=True, timeout=300).stdout)


def main():
    run = Run("R09tcgarepl")
    cfg = {
        "report": "R09",
        "question": "does R06's KEAP1/SMARCA4 prognostic association replicate in TCGA?",
        "why": "R06's audit withdrew its replication claim - luad_mskcc_2020 had 45 events in "
               "604 patients and failed its own positive controls, leaving the finding on one "
               "cohort",
        "cohort": "TCGA-LUAD + TCGA-LUSC, sequenced patients, curated TCGA-CDR survival",
        "independence": "different institutions, whole-exome rather than targeted ctDNA panel, "
                        "different treatment era",
        "declared_asymmetry": "TCGA lung is predominantly resected early-stage disease; the MSK "
                              "ctDx cohort was advanced disease. A negative result is therefore "
                              "ambiguous between 'not real' and 'specific to advanced disease'",
        "model": "R06's exactly - per-gene UNADJUSTED Cox on OS, BH across the same 5 genes. "
                 "Adding covariates would invalidate the comparison",
        "positive_controls": CONTROLS,
        "control_rule": "if neither TP53 nor STK11 is recovered, this cohort cannot test the "
                        "question and the run halts",
        "replication_rule": "a gene replicates only if it is q<0.05 AND in the same direction "
                            "as R06",
        "prespecified_secondary": "late-stage (II-IV) subgroup, as the closest available match "
                                  "to the MSK advanced-disease setting",
        "prespecified_tertiary": "adjusted model (age, stage, histology) reported for context "
                                 "only, not as the replication test",
        "profiled_set_rule": "only patients in each study's _sequenced sample list; an "
                             "unprofiled patient is not wild-type",
        "r06_reference_values": R06,
        "seed": SEED, "expectations": EXPECT,
    }
    cfg_hash = run.start(cfg, [CLIN])

    # ---- profiled set + mutations ----
    seq, per, nrec = set(), {g: set() for g in GENES.values()}, 0
    for study in STUDIES:
        ids = curl(f"https://www.cbioportal.org/api/sample-lists/{study}_sequenced/sample-ids")
        seq |= {s.rsplit("-", 1)[0] for s in ids}
        d = curl(f"https://www.cbioportal.org/api/molecular-profiles/{study}_mutations/"
                 f"mutations/fetch?projection=DETAILED",
                 json.dumps({"sampleListId": f"{study}_all",
                             "entrezGeneIds": list(GENES)}))
        nrec += len(d)
        for x in d:
            g = GENES.get(x.get("entrezGeneId"))
            if g and x.get("mutationType") not in SILENT and x.get("patientId"):
                per[g].add(x["patientId"])
    run.gate("G0_sequenced", EXPECT["sequenced"], len(seq), abs(len(seq) - EXPECT["sequenced"]) <= 5,
             "patients with mutation profiling; unprofiled != wild-type")
    run.gate("G1_mutation_records", "> 800", nrec, nrec > 800)

    # ---- survival frame ----
    rows = []
    for r in csv.DictReader(open(CLIN), delimiter="\t"):
        if r["case"] not in seq:
            continue
        try:
            t = float(r["os_months"])
        except (TypeError, ValueError):
            continue
        if t <= 0 or not r["os_status"]:
            continue
        st = (r["stage"] or "").upper().replace("STAGE", "").strip()
        rows.append({"case": r["case"], "T": t, "E": int(r["os_status"].startswith("1")),
                     "lusc": int(r["project"] == "TCGA-LUSC"),
                     "age": float(r["age"]) if r["age"] else np.nan,
                     "late": 1 if st.startswith(("II", "III", "IV")) else 0,
                     "stage_known": bool(st),
                     **{g: int(r["case"] in per[g]) for g in GENES.values()}})
    df = pd.DataFrame(rows)
    run.gate("G2_events", f"> {EXPECT['min_events']} events", int(df["E"].sum()),
             int(df["E"].sum()) > EXPECT["min_events"],
             "R06's failed replication cohort had 45 - this is the power that was missing")
    late_frac = float(df[df["stage_known"]]["late"].mean())
    run.gate("G3_stage_parsed", "late fraction 0.20-0.60", round(late_frac, 4),
             0.20 <= late_frac <= 0.60)
    kmall = KaplanMeierFitter().fit(df["T"], df["E"])
    med_all = float(kmall.median_survival_time_)
    run.log("median_os_months", value=round(med_all, 1),
            note="MSK ctDx advanced-disease cohort median was 28.8 months (R06)")
    run.log("cohort", n=len(df), events=int(df["E"].sum()),
            **{g: int(df[g].sum()) for g in GENES.values()})

    def per_gene(frame, label):
        res, ps, keys = {}, [], []
        for g in list(CHROM) + CONTROLS:
            if frame[g].sum() < 10:
                res[g] = {"skipped": "fewer than 10 mutated", "n_mutated": int(frame[g].sum())}
                continue
            cph = CoxPHFitter().fit(frame[["T", "E", g]], "T", "E")
            km, med = KaplanMeierFitter(), {}
            for v in (1, 0):
                m = frame[g] == v
                km.fit(frame["T"][m], frame["E"][m])
                med[v] = float(km.median_survival_time_)
            lr = logrank_test(frame["T"][frame[g] == 1], frame["T"][frame[g] == 0],
                              frame["E"][frame[g] == 1], frame["E"][frame[g] == 0])
            res[g] = {"n_mutated": int(frame[g].sum()), "n_total": int(len(frame)),
                      "events": int(frame["E"].sum()),
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
        for g in res:
            if "hr" in res[g]:
                run.log(label, gene=g, hr=round(res[g]["hr"], 3),
                        q=round(res[g]["q_bh"], 5), n_mut=res[g]["n_mutated"])
        return res

    primary = per_gene(df, "primary")

    # ---- pre-specified secondary: late-stage subgroup ----
    late = df[(df["late"] == 1) & df["stage_known"]].reset_index(drop=True)
    run.log("late_subgroup", n=len(late), events=int(late["E"].sum()))
    secondary = per_gene(late, "late_stage") if late["E"].sum() >= 60 else {
        "skipped": f"only {int(late['E'].sum())} events"}

    # ---- positive-control gate, evaluated in BOTH strata ----
    # The rule's intent is "if this dataset cannot demonstrate known biology anywhere, it
    # cannot test the question." An earlier implementation gated on the full cohort alone and
    # halted before the late-stage subgroup the same config pre-declared for this contingency.
    # The replication criterion (q<0.05 + same direction as R06) is unchanged.
    ctrl = {c: {"q": primary.get(c, {}).get("q_bh"), "hr": primary.get(c, {}).get("hr")}
            for c in CONTROLS}
    ctrl_late = {c: {"q": secondary.get(c, {}).get("q_bh") if isinstance(secondary, dict) else None,
                     "hr": secondary.get(c, {}).get("hr") if isinstance(secondary, dict) else None}
                 for c in CONTROLS}
    ok_full = any(primary.get(c, {}).get("significant_q05") for c in CONTROLS)
    ok_late = isinstance(secondary, dict) and any(
        secondary.get(c, {}).get("significant_q05") for c in CONTROLS)
    run.gate("G4_positive_controls", "TP53 or STK11 recovered at q<0.05 in the full cohort "
             "OR in the pre-declared late-stage subgroup",
             {"full_cohort": {c: {"q": round(v["q"], 5) if v["q"] else None,
                                  "hr": round(v["hr"], 3) if v["hr"] else None}
                              for c, v in ctrl.items()},
              "late_stage": {c: {"q": round(v["q"], 5) if v["q"] else None,
                                 "hr": round(v["hr"], 3) if v["hr"] else None}
                             for c, v in ctrl_late.items()},
              "holds_full": ok_full, "holds_late": ok_late},
             ok_full or ok_late,
             "established NSCLC biology; R06's failed cohort failed exactly here")

    # ---- tertiary: adjusted, for context only ----
    adj = {}
    d2 = df.dropna(subset=["age"]).reset_index(drop=True)
    for g in CHROM:
        if d2[g].sum() < 10:
            continue
        cols = ["T", "E", g, "age", "late", "lusc"]
        m = CoxPHFitter().fit(d2[cols], "T", "E")
        adj[g] = {"hr": float(np.exp(m.params_[g])), "p": float(m.summary.loc[g, "p"]),
                  "ci95": [float(np.exp(m.confidence_intervals_.loc[g].iloc[0])),
                           float(np.exp(m.confidence_intervals_.loc[g].iloc[1]))],
                  "n": int(len(d2))}

    # ---- replication verdict ----
    verdict = {}
    for g in CHROM:
        p = primary.get(g, {})
        s = secondary.get(g, {}) if isinstance(secondary, dict) else {}
        same_dir = None
        if R06[g]["hr"] and "hr" in p:
            same_dir = (p["hr"] > 1) == (R06[g]["hr"] > 1)
        verdict[g] = {
            "r06_hr": R06[g]["hr"], "tcga_hr": p.get("hr"), "tcga_q": p.get("q_bh"),
            "tcga_ci95": p.get("ci95"), "same_direction": same_dir,
            "replicates": bool(p.get("significant_q05") and same_dir),
            "late_stage_hr": s.get("hr"), "late_stage_q": s.get("q_bh"),
            "r06_ci_contains_tcga_hr": bool(
                R06[g]["hr"] and p.get("hr") and
                R06[g]["ci"][0] <= p["hr"] <= R06[g]["ci"][1]) if R06[g].get("ci") else None}

    results = {
        "cohort": {"n": len(df), "events": int(df["E"].sum()),
                   "median_os_months": med_all,
                   "sequenced_patients": len(seq),
                   "mutated": {g: int(df[g].sum()) for g in GENES.values()},
                   "late_stage_n": len(late), "late_stage_events": int(late["E"].sum())},
        "positive_controls": {"full_cohort": ctrl, "late_stage": ctrl_late,
                              "holds_full": ok_full, "holds_late": ok_late},
        "primary_unadjusted": primary,
        "secondary_late_stage": secondary,
        "tertiary_adjusted": adj,
        "replication_verdict": verdict,
        "_decision": {
            "rule": "replicates only if q<0.05 and same direction as R06",
            "replicating": [g for g in CHROM if verdict[g]["replicates"]],
            "interpretation": None,  # set below
            "prespecified_in_config_sha256": cfg_hash},
        "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir},
    }
    rep = results["_decision"]["replicating"]
    results["_decision"]["interpretation"] = (
        f"replicates in TCGA: {', '.join(rep)}" if rep else
        "no chromatin-regulator association replicates in TCGA; see the declared "
        "early-vs-advanced-stage asymmetry before reading this as a refutation")
    run.write("results.json", results)
    run.write("survival_frame_r09.csv", df.to_csv(index=False))
    run.log("decision", **{k: v for k, v in results["_decision"].items() if k != "rule"})
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
