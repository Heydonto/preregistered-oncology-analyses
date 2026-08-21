#!/usr/bin/env python3
"""R24 -- does epigenetic modulation reverse the drug-resistance programme?

Executes PREREG_PLAN_R24.yaml, SHA-256 bffbaad89936c522..., hashed after reading only the TPM
column headers and before any expression value was compared between groups.

SMARCA4/BRG1 is the ATPase of the SWI/SNF chromatin remodelling complex, so shSMARCA4 in an
osimertinib-resistant line is epigenetic modulation of a resistant state -- the manuscript's
thesis, in a real experiment.
"""
import collections, gzip, hashlib, json, os, sys
import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, "/Users/rezanehzati/Projects/quantara/paper-draft/02-revised-protocol-open-datasets")
from audit_kit import Run

TPM = "/Users/rezanehzati/quantara-staging/resist/GSE202859_TPMs_allsamples.txt.gz"
PLAN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PREREG_PLAN_R24.yaml")
PLAN_SHA = "bffbaad89936c5222c6e7fec678a88aba5c613358c8a8502ab39a7a83baa4087"
SEED, N_PERM, FC = 20260820, 2000, 1.0

PAIRS = {"PC9":    (["PC9.1","PC9.2","PC9.3"], ["PC9.OR.1","PC9.OR.2","PC9.OR.3"]),
         "H1975":  (["H1975.1","H1975.2","H1975.3"], ["H1975.OR.1","H1975.OR.2","H1975.OR.3"]),
         "HCC827": (["HCC827.1","HCC827.2","HCC827.3"], ["HCC827.OR.1","HCC827.OR.2","HCC827.OR.3"])}
KD = {"PC9OR": (["PC9OR.SCR.UNT.1","PC9OR.SCR.UNT.2","PC9OR.SCR.OSI.1","PC9OR.SCR.OSI.2"],
                ["PC9OR.B1A.UNT.1","PC9OR.B1A.UNT.2","PC9OR.B1A.OSI.1","PC9OR.B1A.OSI.2",
                 "PC9OR.B1B.UNT.1","PC9OR.B1B.UNT.2","PC9OR.B1B.OSI.1","PC9OR.B1B.OSI.2"]),
      "YU005C": ([f"YU005C.SCR.{t}.{i}" for t in ("UNT","OSI") for i in (1,2,3)],
                 [f"YU005C.{b}.{t}.{i}" for b in ("B1A","B1B") for t in ("UNT","OSI") for i in (1,2,3)])}


def main():
    run = Run("R24smarca4")
    got = hashlib.sha256(open(PLAN, "rb").read()).hexdigest()
    cfg = run.start({"executes": "PREREG_PLAN_R24.yaml", "declared_sha256": PLAN_SHA,
                     "observed_sha256": got, "seed": SEED, "n_perm": N_PERM, "fc_cut": FC,
                     "isolation": "hash-before-outcome only; author and executor are the same",
                     "status": "HELD FOR IP - not for publication pending patent"}, [PLAN])
    run.gate("G0_plan_unmodified", PLAN_SHA, got, got == PLAN_SHA)

    fh = gzip.open(TPM, "rt", errors="ignore")
    cols = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(cols) if c}
    need = [s for v in PAIRS.values() for g in v for s in g] + \
           [s for v in KD.values() for g in v for s in g]
    missing = [s for s in need if s not in idx]
    run.gate("G1_libraries_present", "every named library is in the TPM file", missing, not missing)

    sym, rows = [], []
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < len(cols):
            continue
        try:
            rows.append([float(p[idx[s]]) for s in need])
        except ValueError:
            continue
        sym.append(p[2])
    M = np.array(rows, np.float64)
    pos = {s: i for i, s in enumerate(need)}
    keep = (M > 1).sum(1) >= 3
    M, sym = M[keep], [s for s, k in zip(sym, keep) if k]
    # collapse duplicate symbols by max mean expression
    best = {}
    for i, s in enumerate(sym):
        m = M[i].mean()
        if s not in best or m > best[s][1]:
            best[s] = (i, m)
    sel = sorted(v[0] for v in best.values())
    M, sym = M[sel], [sym[i] for i in sel]
    L = np.log2(M + 1.0)
    run.gate("G2_genes", "10,000-30,000 genes retained", len(sym), 10000 <= len(sym) <= 30000)

    def lfc(a, b):
        return L[:, [pos[x] for x in b]].mean(1) - L[:, [pos[x] for x in a]].mean(1)

    # ---- PC1: did the knockdown work? ----
    si = sym.index("SMARCA4") if "SMARCA4" in sym else None
    kd_pc9 = lfc(*KD["PC9OR"])
    kd_yu = lfc(*KD["YU005C"])
    run.gate("PC1_knockdown_worked", "SMARCA4 log2FC <= -0.5 in PC9OR shSMARCA4 vs scramble",
             None if si is None else round(float(kd_pc9[si]), 4),
             si is not None and kd_pc9[si] <= -0.5,
             "plan: if the knockdown did not work, halt and interpret nothing downstream")
    run.log("PC1_detail", smarca4_lfc_PC9OR=round(float(kd_pc9[si]), 4),
            smarca4_lfc_YU005C=round(float(kd_yu[si]), 4))

    # ---- PC2 + H1: resistance programmes ----
    res_lfc = {k: lfc(*v) for k, v in PAIRS.items()}
    ndeg = {k: int((np.abs(v) > FC).sum()) for k, v in res_lfc.items()}
    run.gate("PC2_resistance_is_real", "at least 200 genes |log2FC|>1 in each line", ndeg,
             all(v >= 200 for v in ndeg.values()))
    lines = list(PAIRS)
    h1 = {}
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            r = float(spearmanr(res_lfc[lines[i]], res_lfc[lines[j]]).statistic)
            h1[f"{lines[i]}_vs_{lines[j]}"] = r
    run.log("H1_consistency", **{k: round(v, 4) for k, v in h1.items()})

    # ---- H2 PRIMARY: does knockdown reverse the resistance programme? ----
    up = np.where(res_lfc["PC9"] > FC)[0]
    dn = np.where(res_lfc["PC9"] < -FC)[0]
    run.gate("G4_set_provenance",
             "resistance sets defined from the parental pair only, never from knockdown samples",
             {"n_up": len(up), "n_down": len(dn), "defined_from": "PC9 vs PC9.OR"}, True)
    run.gate("G5_perm_floor", "permutation floor below 0.05",
             round(1 / (N_PERM + 1), 5), 1 / (N_PERM + 1) < 0.05)

    rng = np.random.default_rng(SEED)

    def test(kd, s):
        obs = float(kd[s].mean())
        null = np.array([kd[rng.choice(len(kd), len(s), replace=False)].mean()
                         for _ in range(N_PERM)])
        return {"n": len(s), "mean_kd_lfc": obs, "null_mean": float(null.mean()),
                "null_lo": float(np.percentile(null, 2.5)),
                "null_hi": float(np.percentile(null, 97.5)),
                "p_two_sided": float((1 + int((np.abs(null - null.mean())
                                               >= abs(obs - null.mean())).sum())) / (N_PERM + 1))}

    res = {}
    for model, kd in (("PC9OR", kd_pc9), ("YU005C", kd_yu)):
        res[model] = {"resistance_up": test(kd, up), "resistance_down": test(kd, dn),
                      "all_genes_mean": float(kd.mean())}
        for s in ("resistance_up", "resistance_down"):
            d = res[model][s]
            run.log("H2" if model == "PC9OR" else "H3", model=model, set=s, n=d["n"],
                    mean_kd_lfc=round(d["mean_kd_lfc"], 4),
                    null_ci=[round(d["null_lo"], 4), round(d["null_hi"], 4)],
                    p=round(d["p_two_sided"], 5))

    run.gate("G3_nulls_centred", "permutation nulls centred on the all-gene mean (within 0.02)",
             round(max(abs(res[m][s]["null_mean"] - res[m]["all_genes_mean"])
                       for m in res for s in ("resistance_up", "resistance_down")), 4),
             all(abs(res[m][s]["null_mean"] - res[m]["all_genes_mean"]) < 0.02
                 for m in res for s in ("resistance_up", "resistance_down")))

    p = res["PC9OR"]
    up_rev = p["resistance_up"]["mean_kd_lfc"] < p["resistance_up"]["null_lo"]
    dn_rev = p["resistance_down"]["mean_kd_lfc"] > p["resistance_down"]["null_hi"]
    up_rei = p["resistance_up"]["mean_kd_lfc"] > p["resistance_up"]["null_hi"]
    dn_rei = p["resistance_down"]["mean_kd_lfc"] < p["resistance_down"]["null_lo"]
    if up_rev and dn_rev:
        verdict = "THESIS_SUPPORTED_knockdown_reverses_the_resistance_programme"
    elif up_rei and dn_rei:
        verdict = "THESIS_CONTRADICTED_knockdown_reinforces_it"
    elif up_rev or dn_rev:
        verdict = "PARTIAL_one_direction_only"
    else:
        verdict = "NULL_no_reversal"
    run.log("VERDICT", primary=verdict)

    out = {"status": "HELD FOR IP - not for publication pending patent",
           "plan_sha256": got,
           "VERDICT": {"H2_primary": verdict,
                       "rules_pre_declared_in": "PREREG_PLAN_R24.yaml " + PLAN_SHA[:16]},
           "genes_retained": len(sym),
           "PC1_smarca4_knockdown": {"PC9OR_log2FC": float(kd_pc9[si]),
                                     "YU005C_log2FC": float(kd_yu[si])},
           "PC2_deg_counts": ndeg, "H1_cross_line_consistency": h1,
           "H2_PC9OR": res["PC9OR"], "H3_YU005C_replication": res["YU005C"],
           "resistance_sets": {"n_up": int(len(up)), "n_down": int(len(dn)),
                               "defined_from": "PC9 parental vs PC9.OR only"},
           "_provenance": {"config_sha256": cfg, "run_dir": run.dir}}
    np.save(os.path.join(run.dir, "kd_lfc_PC9OR.npy"), kd_pc9)
    np.save(os.path.join(run.dir, "res_lfc_PC9.npy"), res_lfc["PC9"])
    open(os.path.join(run.dir, "genes.txt"), "w").write("\n".join(sym))
    run.write("results.json", out)
    run.finalize()
    print("\nVERDICT:", verdict)
    print("  SMARCA4 knockdown log2FC: PC9OR %.3f  YU005C %.3f" % (kd_pc9[si], kd_yu[si]))
    for s in ("resistance_up", "resistance_down"):
        d = p[s]
        print(f"  {s:16s} n={d['n']:5d} mean kd log2FC {d['mean_kd_lfc']:+.4f} "
              f"null 95% [{d['null_lo']:+.4f},{d['null_hi']:+.4f}] p={d['p_two_sided']:.5f}")
    print("RUN DIR:", run.dir)


if __name__ == "__main__":
    main()
