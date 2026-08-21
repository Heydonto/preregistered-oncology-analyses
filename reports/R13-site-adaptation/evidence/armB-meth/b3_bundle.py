#!/usr/bin/env python3
"""Assemble the Arm B methylation-track evidence bundle (R-series layout) and mirror it.

usage: python3 b3_bundle.py <stage1_dir> <stage3_dir> <stage4_dir>
"""
import json, os, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from audit_kit import sha256                                            # noqa: E402

EV = os.path.join(HERE, "evidence")
WORK = os.path.join(HERE, "work")
GCS = "gs://heydonto-quantara-lungcdx/nsclc-rwpr-study/armB-meth/"
GCS_PROV = "gs://heydonto-quantara-lungcdx/nsclc-rwpr-study/_provenance/unsealing_events.jsonl"


def main(d1, d3, d4):
    os.makedirs(EV, exist_ok=True)
    # the plan, as executed and hashed before any label was read
    shutil.copy(f"{d4}/config.yaml", f"{EV}/config.yaml")
    for src, dst in ((f"{d4}/results.json", "results.json"),
                     (f"{d4}/audit_armb_meth.json", "audit_armb_meth.json"),
                     (f"{d4}/permutation_null.json", "permutation_null.json"),
                     (f"{d4}/oof_epic_scores.csv", "oof_epic_scores.csv"),
                     (f"{d4}/per_sample_qc.csv", "per_sample_qc.csv"),
                     (f"{d4}/env.txt", "env.txt"),
                     (f"{d1}/stage1_results.json", "stage1_tcga_results.json"),
                     (f"{d3}/stage3_parse.json", "stage3_deposit_parse.json"),
                     (f"{d3}/epic_metadata_raw.json", "epic_deposited_metadata.json"),
                     (f"{WORK}/unsealing_event_2.json", "unsealing_event_2.json"),
                     (f"{WORK}/oof_keap1_repro_full.csv", "oof_keap1_repro_tcga_fullarray.csv"),
                     (f"{WORK}/oof_keap1_overlap.csv", "oof_keap1_tcga_epicoverlap.csv"),
                     (f"{WORK}/oof_sex_repro_full.csv", "oof_sex_control_tcga.csv"),
                     (f"{WORK}/oof_os_cox_tcga.csv", "oof_os_cox_tcga.csv"),
                     (f"{WORK}/manifest_overlap.json", "manifest_overlap_counts.json")):
        if os.path.exists(src):
            shutil.copy(src, f"{EV}/{dst}")

    # merged gates across all three stages, in execution order
    gates = []
    for d, stage in ((d1, "stage1_tcga"), (d3, "stage3_unseal_parse"), (d4, "stage4_eval")):
        if os.path.exists(f"{d}/gates.json"):
            for g in json.load(open(f"{d}/gates.json")):
                gates.append({**g, "stage": stage})
    json.dump(gates, open(f"{EV}/gates.json", "w"), indent=2, default=str)

    # merged log
    with open(f"{EV}/log.jsonl", "w") as out:
        for d in (d1, d3, d4):
            if os.path.exists(f"{d}/log.jsonl"):
                out.write(open(f"{d}/log.jsonl").read())

    # merged input hashes
    inp = {}
    for d in (d1, d3, d4):
        if os.path.exists(f"{d}/inputs.json"):
            inp.update(json.load(open(f"{d}/inputs.json")))
    json.dump(inp, open(f"{EV}/inputs.json", "w"), indent=2)

    # snapshot of the append-only unsealing log as it stands
    open(f"{EV}/unsealing_events_snapshot.jsonl", "w").write(
        subprocess.run(["gsutil", "cat", GCS_PROV], capture_output=True, text=True).stdout)

    # the frozen deployable models (coefficients only — no patient data)
    art = json.load(open(f"{WORK}/stage1_artifacts.json"))
    json.dump({"note": "coefficients of the models frozen before unsealing; probe ids and "
                       "weights only, no patient-level values",
               "models": {nm: {"probes": (m.get("selected_probes") or m.get("probes")),
                               "coef": m["coef"], "mean": m["mean"], "scale": m["scale"],
                               "intercept": m.get("intercept"), "C": m.get("C"),
                               "n_nonzero": m["n_nonzero"]}
                          for nm, m in art["models"].items()}},
              open(f"{EV}/frozen_models.json", "w"))
    json.dump({"keap1_active_probes":
               json.load(open(f"{d4}/results.json"))["model_probe_coverage_on_epic"][
                   "active_probe_identity_check"]},
              open(f"{EV}/keap1_active_probes.json", "w"), indent=2)

    for f in ("b1_armb_meth.py", "b2_eval.py", "audit_armb_meth.py", "audit_kit.py",
              "b3_bundle.py"):
        shutil.copy(os.path.join(HERE, f), f"{EV}/{f}")

    L = [f"{sha256(os.path.join(EV, f))}  {f}" for f in sorted(os.listdir(EV))
         if f != "MANIFEST.sha256"]
    open(f"{EV}/MANIFEST.sha256", "w").write("\n".join(L) + "\n")
    print(f"{len(L)} files in {EV}")
    for line in L:
        print(" ", line[:12], line.split("  ", 1)[1])

    r = subprocess.run(["gsutil", "-m", "rsync", "-r", EV, GCS + "evidence/"],
                       capture_output=True, text=True)
    print("mirror rc", r.returncode, r.stderr.strip()[-400:])
    for f in ("b1_armb_meth.py", "b2_eval.py", "audit_armb_meth.py"):
        subprocess.run(["gsutil", "cp", os.path.join(HERE, f), GCS + "code/" + f],
                       capture_output=True)
    print(subprocess.run(["gsutil", "ls", "-r", GCS], capture_output=True, text=True).stdout)


if __name__ == "__main__":
    main(*sys.argv[1:4])
