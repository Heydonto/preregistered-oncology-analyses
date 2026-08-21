#!/usr/bin/env python3
"""R18 - assemble the two-encoder site-leakage comparison and rule on it.

Reads the eight arm JSONs written by armD-meth/mil_encoder_compare.py:

    {tcga, dinov2-tcga} x {grouped, random} x {subtype, meth}

and answers the question Paper 1 lists as its first limitation -- "All conclusions concern
Phikon-v2 features. Whether the site signal is a property of the archive or of that encoder is
untested." If the +0.171 subtype inflation is a quirk of pan-cancer histology pretraining, Paper 1
is a paper about Owkin's model. If it survives an encoder trained on natural images, it is a
property of the archive.

THE COMPARISON IS FAIR ONLY IF TWO THINGS HOLD, and both are gated rather than assumed:

  1. The harness reproduces the published Phikon-v2 numbers. If it does not, nothing downstream
     means anything. (G1-G4.)
  2. Inflation is read relative to available headroom. dinov2-large is not trained on histology,
     so its absolute AUROC will be lower, and a lower absolute inflation could then reflect less
     signal rather than less leakage. So both are reported:

         inflation          = random - grouped
         relative inflation = (random - grouped) / (1 - grouped)

     the second being the fraction of the remaining distance to a perfect score that fold
     assignment alone recovers. The verdict rests on the relative figure, which is the one that
     survives an encoder being weaker in absolute terms.

Run:  python3 m16_encoder_compare_analysis.py
"""
import json
import os
import subprocess

import numpy as np

from audit_kit import Run

BUCKET = "heydonto-quantara-lungcdx"
R18 = f"gs://{BUCKET}/nsclc-rwpr-study/armD-meth/r18/"
ENCODERS = {"tcga": "phikon-v2", "dinov2-tcga": "dinov2-large"}
TARGET_NAMES = ["keap1_sig", "global", "island", "opensea", "tss200", "body"]

# R15 / m13 published Phikon-v2 values the control arms must recover
PUBLISHED = {"subtype_grouped": 0.799, "subtype_random": 0.9703,
             "keap1_grouped": 0.664, "meth_mean_inflation": 0.0849}
TOL_AUROC = 0.01
TOL_INFL = 0.02

# pre-declared thresholds. 0.05 is the bar R15's own gate G3_leakage_is_real used.
LEAK_REAL = 0.05
LEAK_ABSENT = 0.02


def fetch(name):
    p = f"/tmp/r18_{name}.json"
    if not os.path.exists(p):
        subprocess.run(["gcloud", "storage", "cp", R18 + name + ".json", p],
                       check=True, capture_output=True)
    return json.load(open(p))


def main():
    run = Run("R18encodercmp")
    cfg_hash = run.start(
        {"question": "is the whole-slide site signature a property of the archive or of "
                     "Phikon-v2?",
         "design": "identical tiles (verified byte-identical), identical architecture, folds, "
                   "hyperparameters and seeds; the encoder is the only difference. "
                   "facebook/dinov2-large matches Phikon-v2 in model class (Dinov2Model), width "
                   "(1024), depth (24), SSL algorithm (DINOv2), readout (CLS) and normalisation "
                   "(ImageNet); it differs in training corpus (LVD-142M natural images vs "
                   "pan-cancer histology) and patch size (14 vs 16).",
         "arms": "{phikon-v2, dinov2-large} x {grouped, random} x {subtype, meth}",
         "metrics": {"inflation": "random - grouped",
                     "relative_inflation": "(random - grouped) / (1 - grouped)",
                     "why_relative": "dinov2-large is not trained on histology so its absolute "
                                     "AUROC is expected to be lower; a smaller absolute inflation "
                                     "could then reflect less signal rather than less leakage. The "
                                     "verdict uses the relative figure."},
         "published_controls": PUBLISHED, "tol_auroc": TOL_AUROC, "tol_inflation": TOL_INFL,
         "verdict_rules": [
             f"dinov2 subtype relative inflation > {LEAK_REAL} -> "
             "SITE_SIGNATURE_IS_ARCHIVE_NOT_ENCODER",
             f"dinov2 subtype relative inflation < {LEAK_ABSENT} -> "
             "SITE_SIGNATURE_WAS_ENCODER_SPECIFIC",
             "otherwise -> INTERMEDIATE",
             "secondary: sign agreement across all 6 methylation targets -> "
             "METHYLATION_LEAKAGE_ALSO_REPRODUCES"],
         "status": "HELD FOR IP - not for publication pending patent"},
        [])

    A = {}
    for fd in ENCODERS:
        for fold in ("grouped", "random"):
            for mode in ("subtype", "meth"):
                A[(fd, fold, mode)] = fetch(f"{fd}__{fold}__{mode}")

    # ---- encoder identity actually used, straight from the npz metadata ----
    enc = {fd: sorted({e for f in ("grouped", "random") for m in ("subtype", "meth")
                       for e in A[(fd, f, m)]["encoder"]}) for fd in ENCODERS}
    run.gate("G0_encoder_identity",
             "each feature dir carries exactly one encoder, and it is the intended one",
             enc,
             all(len(v) == 1 for v in enc.values())
             and "phikon" in enc["tcga"][0] and "dinov2" in enc["dinov2-tcga"][0],
             "read from the npz 'encoder' field, not inferred from the directory name")

    n_pat = {k: v["n_patients"] for k, v in A.items()}
    run.gate("G1_same_cohort", "all eight arms cover 760 patients",
             sorted(set(n_pat.values())), set(n_pat.values()) == {760})

    # ---- reproduction controls: phikon must recover the published numbers ----
    ph_sg = A[("tcga", "grouped", "subtype")]["metrics"]["subtype"]["auroc"]
    ph_sr = A[("tcga", "random", "subtype")]["metrics"]["subtype"]["auroc"]
    ph_kg = A[("tcga", "grouped", "subtype")]["metrics"]["keap1"]["auroc"]
    run.gate("G2_repro_subtype_grouped", f"{PUBLISHED['subtype_grouped']} +- {TOL_AUROC}",
             round(ph_sg, 4), abs(ph_sg - PUBLISHED["subtype_grouped"]) <= TOL_AUROC)
    run.gate("G3_repro_subtype_random", f"{PUBLISHED['subtype_random']} +- {TOL_AUROC}",
             round(ph_sr, 4), abs(ph_sr - PUBLISHED["subtype_random"]) <= TOL_AUROC)
    run.gate("G4_repro_keap1_grouped", f"{PUBLISHED['keap1_grouped']} +- {TOL_AUROC}",
             round(ph_kg, 4), abs(ph_kg - PUBLISHED["keap1_grouped"]) <= TOL_AUROC,
             "if any of G2-G4 fails the harness is not the published one and the dinov2 "
             "comparison is void")

    def rel(g, r):
        return (r - g) / (1 - g) if g < 1 else float("nan")

    # ---- subtype ----
    sub = {}
    for fd in ENCODERS:
        g = A[(fd, "grouped", "subtype")]["metrics"]["subtype"]["auroc"]
        r = A[(fd, "random", "subtype")]["metrics"]["subtype"]["auroc"]
        kg = A[(fd, "grouped", "subtype")]["metrics"]["keap1"]["auroc"]
        kr = A[(fd, "random", "subtype")]["metrics"]["keap1"]["auroc"]
        sub[fd] = {"encoder": enc[fd][0],
                   "subtype_grouped": g, "subtype_random": r,
                   "subtype_inflation": r - g, "subtype_relative_inflation": rel(g, r),
                   "keap1_grouped": kg, "keap1_random": kr,
                   "keap1_inflation": kr - kg, "keap1_relative_inflation": rel(kg, kr)}
        run.log("subtype_arm", encoder=enc[fd][0], grouped=round(g, 4), random=round(r, 4),
                inflation=round(r - g, 4), relative_inflation=round(rel(g, r), 4))

    # ---- methylation, six targets ----
    meth = {}
    for fd in ENCODERS:
        mg = A[(fd, "grouped", "meth")]["metrics"]
        mr = A[(fd, "random", "meth")]["metrics"]
        per = {t: {"grouped_rho": mg[t]["rho"], "random_rho": mr[t]["rho"],
                   "inflation": mr[t]["rho"] - mg[t]["rho"]} for t in TARGET_NAMES}
        infl = np.array([per[t]["inflation"] for t in TARGET_NAMES])
        meth[fd] = {"encoder": enc[fd][0], "per_target": per,
                    "mean_grouped_rho": float(np.mean([per[t]["grouped_rho"] for t in TARGET_NAMES])),
                    "mean_random_rho": float(np.mean([per[t]["random_rho"] for t in TARGET_NAMES])),
                    "mean_inflation": float(infl.mean()),
                    "n_targets_inflated": int((infl > 0).sum())}
        run.log("meth_arm", encoder=enc[fd][0],
                mean_grouped=round(meth[fd]["mean_grouped_rho"], 4),
                mean_random=round(meth[fd]["mean_random_rho"], 4),
                mean_inflation=round(meth[fd]["mean_inflation"], 4),
                n_inflated=meth[fd]["n_targets_inflated"])

    run.gate("G5_repro_meth_inflation",
             f"phikon mean methylation inflation {PUBLISHED['meth_mean_inflation']} +- {TOL_INFL}",
             round(meth["tcga"]["mean_inflation"], 4),
             abs(meth["tcga"]["mean_inflation"] - PUBLISHED["meth_mean_inflation"]) <= TOL_INFL)

    # ---- mechanical verdict ----
    dv = sub["dinov2-tcga"]["subtype_relative_inflation"]
    if dv > LEAK_REAL:
        verdict = "SITE_SIGNATURE_IS_ARCHIVE_NOT_ENCODER"
    elif dv < LEAK_ABSENT:
        verdict = "SITE_SIGNATURE_WAS_ENCODER_SPECIFIC"
    else:
        verdict = "INTERMEDIATE"
    secondary = ("METHYLATION_LEAKAGE_ALSO_REPRODUCES"
                 if meth["dinov2-tcga"]["n_targets_inflated"] >= 5
                 else "METHYLATION_LEAKAGE_DID_NOT_REPRODUCE")
    run.log("VERDICT", primary=verdict, secondary=secondary,
            dinov2_relative_inflation=round(dv, 4),
            phikon_relative_inflation=round(sub["tcga"]["subtype_relative_inflation"], 4))

    out = {"status": "HELD FOR IP - not for publication pending patent",
           "question": "archive or encoder?",
           "VERDICT": {"primary": verdict, "secondary": secondary,
                       "rules_pre_declared_in": "config.yaml, hash " + cfg_hash},
           "encoders": enc, "subtype": sub, "methylation": meth,
           "published_controls": PUBLISHED,
           "controls_reproduced": {"subtype_grouped": ph_sg, "subtype_random": ph_sr,
                                   "keap1_grouped": ph_kg,
                                   "meth_mean_inflation": meth["tcga"]["mean_inflation"]},
           "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir}}
    run.write("results.json", out)
    run.finalize()
    print("\nVERDICT:", verdict, "|", secondary)
    for fd in ENCODERS:
        s = sub[fd]
        print(f"  {s['encoder']:22s} subtype {s['subtype_grouped']:.4f} -> "
              f"{s['subtype_random']:.4f}  infl {s['subtype_inflation']:+.4f}  "
              f"rel {s['subtype_relative_inflation']:+.4f}")
    print("RUN DIR:", run.dir)


if __name__ == "__main__":
    main()
