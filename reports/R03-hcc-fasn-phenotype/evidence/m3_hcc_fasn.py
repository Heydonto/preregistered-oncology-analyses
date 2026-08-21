#!/usr/bin/env python3
"""R03 — FASN and the lipogenic phenotype in HCC: testing the manuscript's three claims.

The withdrawn HCC manuscript claimed FASN associates with (i) tumour stage,
(ii) vascular invasion and (iii) survival. A prior analysis of TCGA-LIHC found none of
these (stage p=0.22, vascular invasion p=0.42). This report tests the same three claims
independently in the Chinese Liver Cancer Atlas, which has microvascular invasion graded
M0/M1/M2 on 237 of 238 RNA-sequenced tumours.

Positive control: AFP gene expression versus the independently measured serum AFP
concentration. These are different assays of the same biology and must correlate. If they
do not, the expression matrix or the clinical join is wrong and nothing else here is
trustworthy.
"""
import csv, hashlib, json, os, platform, subprocess, sys, time
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu, kruskal

ROOT = os.path.dirname(os.path.abspath(__file__))
D = "/Users/rezanehzati/quantara-staging/staged/labels/hcc/hcc_clca_2024/"
MASTER_SEED = 20260807
LIPO = ["ACACA", "SCD", "SREBF1"]
EXPECT = {"n_expr": 238, "n_mvi": 237, "n_bclc": 238, "n_os": 183,
          "lipo_rho_min": 0.40, "control_rho_min": 0.30}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Run:
    def __init__(self, tag):
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                             capture_output=True).stdout.decode().strip() or "nogit"
        self.dir = os.path.join(ROOT, "runs",
                               time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()) + f"-{tag}-{sha}")
        os.makedirs(self.dir, exist_ok=True)
        self.gates, self.logf = [], open(os.path.join(self.dir, "log.jsonl"), "a")

    def log(self, e, **kw):
        r = {"t": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "event": e, **kw}
        self.logf.write(json.dumps(r, default=str) + "\n"); self.logf.flush()
        print(f"[{r['t']}] {e}: " + " ".join(f"{k}={v}" for k, v in kw.items()), flush=True)

    def gate(self, n, exp, obs, ok, note=""):
        g = {"gate": n, "expected": exp, "observed": obs, "result": "PASS" if ok else "FAIL",
             "note": note, "t": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        self.gates.append(g); self.log("gate", gate=n, expected=exp, observed=obs, result=g["result"])
        json.dump(self.gates, open(os.path.join(self.dir, "gates.json"), "w"), indent=2, default=str)
        if not ok:
            self.log("HALT", reason=f"gate {n} failed"); self.finalize(); sys.exit(2)

    def write(self, n, o):
        p = os.path.join(self.dir, n)
        with open(p, "w") as fh:
            json.dump(o, fh, indent=2, default=str) if n.endswith(".json") else fh.write(o)
        return p

    def finalize(self):
        L = []
        for dp, _, fs in os.walk(self.dir):
            for f in sorted(fs):
                if f != "MANIFEST.sha256":
                    fp = os.path.join(dp, f)
                    L.append(f"{sha256(fp)}  {os.path.relpath(fp, self.dir)}")
        open(os.path.join(self.dir, "MANIFEST.sha256"), "w").write("\n".join(sorted(L)) + "\n")


def numeric(v):
    """Serum markers carry censored values such as '>1210'. Strip the operator: for a
    rank-based test the ceiling value preserves ordering."""
    v = v.strip().replace(">", "").replace("<", "")
    try:
        return float(v)
    except Exception:
        return np.nan


def main():
    run = Run("R03hcc")
    cfg = {
        "report": "R03",
        "question": "Does FASN expression associate with stage, vascular invasion or "
                    "survival status in HCC? (the withdrawn manuscript's three claims)",
        "cohort": "Chinese Liver Cancer Atlas (cBioPortal hcc_clca_2024), RNA-seq TPM",
        "positive_control": "AFP gene expression vs independently measured serum AFP",
        "control_rule": "if the AFP control fails (rho <= 0.30) the join or matrix is wrong "
                        "and no other result here may be reported",
        "primary": "FASN vs microvascular invasion (M0/M1/M2), Kruskal-Wallis",
        "secondary": ["FASN vs BCLC stage", "FASN vs Edmondson grade",
                      "FASN vs OS status", "FASN vs recurrence status",
                      "FASN lipogenic co-expression (replication)"],
        "multiplicity": "Benjamini-Hochberg across the five FASN association tests",
        "prior": "TCGA-LIHC found FASN null for stage (p=0.22) and vascular invasion (p=0.42)",
        "master_seed": MASTER_SEED, "expectations": EXPECT,
    }
    cfg_hash = sha256(run.write("config.yaml", json.dumps(cfg, indent=2)))
    run.log("config_written", sha256=cfg_hash)
    run.write("env.txt", subprocess.run([sys.executable, "-m", "pip", "freeze"],
                                        capture_output=True).stdout.decode())

    files = {f: {"bytes": os.path.getsize(D + f), "sha256": sha256(D + f)}
             for f in ("data_mrna_seq_tpm.txt", "data_clinical_sample.txt",
                       "data_clinical_patient.txt")}
    run.write("inputs.json", files)
    run.gate("G0_inputs", "3 files present", len(files), len(files) == 3)

    # ---- expression -----------------------------------------------------------
    want = {"FASN", "AFP"} | set(LIPO)
    expr, samples = {}, None
    with open(D + "data_mrna_seq_tpm.txt") as fh:
        samples = fh.readline().rstrip("\n").split("\t")[2:]
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0] in want:
                expr[p[0]] = np.array([numeric(x) for x in p[2:]])
    run.gate("G1_expr_samples", EXPECT["n_expr"], len(samples), len(samples) == EXPECT["n_expr"])
    run.gate("G2_genes", sorted(want), sorted(expr), set(expr) == want)

    samp = {r["SAMPLE_ID"]: r for r in csv.DictReader(
        (l for l in open(D + "data_clinical_sample.txt") if not l.startswith("#")), delimiter="\t")}
    pat = {r["PATIENT_ID"]: r for r in csv.DictReader(
        (l for l in open(D + "data_clinical_patient.txt") if not l.startswith("#")), delimiter="\t")}
    run.gate("G3_join", EXPECT["n_expr"], sum(1 for s in samples if s in samp),
             sum(1 for s in samples if s in samp) == EXPECT["n_expr"],
             "every expression sample present in the clinical table")

    def sfield(f):
        return np.array([samp[s].get(f, "").strip() if s in samp else "" for s in samples], object)

    def pfield(f):
        return np.array([pat.get(samp[s]["PATIENT_ID"], {}).get(f, "").strip()
                         if s in samp else "" for s in samples], object)

    fasn = expr["FASN"]

    # ---- POSITIVE CONTROL -----------------------------------------------------
    serum = np.array([numeric(x) for x in sfield("AFP")])
    m = ~np.isnan(expr["AFP"]) & ~np.isnan(serum)
    rho_c, p_c = spearmanr(expr["AFP"][m], serum[m])
    run.gate("G4_positive_control", f"rho > {EXPECT['control_rho_min']}",
             {"rho": round(float(rho_c), 4), "p": f"{p_c:.3e}", "n": int(m.sum())},
             float(rho_c) > EXPECT["control_rho_min"],
             "AFP gene expression vs serum AFP: two assays of the same biology")

    # ---- lipogenic co-expression (replication) --------------------------------
    lipo = {}
    for g in LIPO:
        mm = ~np.isnan(fasn) & ~np.isnan(expr[g])
        r, p = spearmanr(fasn[mm], expr[g][mm])
        lipo[g] = {"rho": float(r), "p": float(p), "n": int(mm.sum())}
    run.gate("G5_lipogenic_replication", f"all rho > {EXPECT['lipo_rho_min']}",
             {g: round(v["rho"], 3) for g, v in lipo.items()},
             all(v["rho"] > EXPECT["lipo_rho_min"] for v in lipo.values()),
             "FASN co-expression with ACACA/SCD/SREBF1 should replicate")

    # ---- the three manuscript claims -----------------------------------------
    tests = {}

    mvi = sfield("MVI")
    ok = np.array([v in ("M0", "M1", "M2") for v in mvi]) & ~np.isnan(fasn)
    groups = [fasn[ok & (mvi == g)] for g in ("M0", "M1", "M2")]
    H, p = kruskal(*groups)
    tests["mvi_M0_M1_M2"] = {"test": "Kruskal-Wallis", "H": float(H), "p": float(p),
                             "n": int(ok.sum()),
                             "medians": {g: float(np.median(x)) for g, x in
                                         zip(("M0", "M1", "M2"), groups)},
                             "group_n": {g: int(len(x)) for g, x in
                                         zip(("M0", "M1", "M2"), groups)}}
    # M0 vs any invasion, the dichotomy the manuscript implies
    a, b = fasn[ok & (mvi == "M0")], fasn[ok & ((mvi == "M1") | (mvi == "M2"))]
    U, p2 = mannwhitneyu(a, b)
    tests["mvi_absent_vs_present"] = {"test": "Mann-Whitney", "p": float(p2),
                                      "n_absent": int(len(a)), "n_present": int(len(b)),
                                      "median_absent": float(np.median(a)),
                                      "median_present": float(np.median(b))}

    bclc = sfield("BCLC")
    ok = np.array([v in ("0", "A", "B", "C") for v in bclc]) & ~np.isnan(fasn)
    gs = [fasn[ok & (bclc == g)] for g in ("0", "A", "B", "C")]
    gs = [g for g in gs if len(g) >= 5]
    H, p = kruskal(*gs)
    tests["bclc_stage"] = {"test": "Kruskal-Wallis", "H": float(H), "p": float(p),
                           "groups_used": len(gs), "n": int(sum(len(g) for g in gs))}

    ed = sfield("EDMONDSON")
    ok = np.array([v.startswith("Level") for v in ed]) & ~np.isnan(fasn)
    gs = [fasn[ok & (ed == g)] for g in sorted(set(ed[ok]))]
    gs = [g for g in gs if len(g) >= 5]
    H, p = kruskal(*gs)
    tests["edmondson_grade"] = {"test": "Kruskal-Wallis", "H": float(H), "p": float(p),
                                "groups_used": len(gs)}

    os_ = pfield("OS_STATUS")
    ok = np.array([v in ("0:LIVING", "1:DECEASED") for v in os_]) & ~np.isnan(fasn)
    a, b = fasn[ok & (os_ == "0:LIVING")], fasn[ok & (os_ == "1:DECEASED")]
    U, p = mannwhitneyu(a, b)
    tests["os_status"] = {"test": "Mann-Whitney", "p": float(p), "n_living": int(len(a)),
                          "n_deceased": int(len(b)), "median_living": float(np.median(a)),
                          "median_deceased": float(np.median(b))}

    rfs = pfield("RFS_STATUS")
    ok = np.array([v.startswith(("0:", "1:")) for v in rfs]) & ~np.isnan(fasn)
    a, b = fasn[ok & (rfs == "0:DiseaseFree")], fasn[ok & (rfs == "1:Recurred/Progressed")]
    U, p = mannwhitneyu(a, b)
    tests["recurrence_status"] = {"test": "Mann-Whitney", "p": float(p),
                                  "n_diseasefree": int(len(a)), "n_recurred": int(len(b))}

    # ---- Benjamini-Hochberg across the five FASN association tests ------------
    keys = ["mvi_M0_M1_M2", "bclc_stage", "edmondson_grade", "os_status", "recurrence_status"]
    ps = np.array([tests[k]["p"] for k in keys])
    order = np.argsort(ps)
    q = np.empty_like(ps)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        i = len(ps) - rank
        prev = min(prev, ps[idx] * len(ps) / i)
        q[idx] = prev
    for k, qq in zip(keys, q):
        tests[k]["q_bh"] = float(qq)
        tests[k]["significant_q05"] = bool(qq < 0.05)

    any_sig = any(tests[k]["significant_q05"] for k in keys)
    results = {"positive_control_afp": {"rho": float(rho_c), "p": float(p_c),
                                        "n": int(m.sum()), "recovered": True},
               "lipogenic_coexpression": lipo,
               "fasn_association_tests": tests,
               "_decision": {"any_association_significant_after_BH": any_sig,
                             "interpretation": ("FASN associates with at least one "
                                                "clinicopathological feature")
                             if any_sig else
                             ("FASN shows no association with stage, vascular invasion, "
                              "grade, survival status or recurrence after correction"),
                             "prespecified_in_config_sha256": cfg_hash},
               "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir,
                               "n_expression_samples": len(samples), "inputs": files}}
    np.savetxt(os.path.join(run.dir, "fasn_by_sample.csv"),
               np.column_stack([fasn, expr["AFP"], serum]), delimiter=",",
               header="FASN_tpm,AFP_tpm,serum_AFP", comments="")
    run.write("results.json", results)
    run.log("decision", **results["_decision"])
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
