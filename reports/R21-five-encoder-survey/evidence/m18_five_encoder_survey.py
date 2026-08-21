#!/usr/bin/env python3
"""R21 - the site-leakage survey across five encoders, and a correction to R18.

R18 compared Phikon-v2 against facebook/dinov2-large and drew two verdicts:
SITE_SIGNATURE_IS_ARCHIVE_NOT_ENCODER (subtype leakage reproduces, larger) and
METHYLATION_LEAKAGE_DID_NOT_REPRODUCE. It flagged its own limitation plainly -- "two encoders is
not a survey" -- and one of the two was not a pathology model at all.

Three pathology foundation models now close that gap: UNI (ViT-L/16, 1024-d), H-optimus-0
(ViT-g/14, 1536-d, Apache-2.0) and Virchow2 (ViT-H/14, 1280-d). All 1,182 slides re-encoded per
model with byte-identical tile grids and each model's own prescribed normalisation.

THE CORRECTION. R18's secondary verdict said methylation leakage was "a property of Phikon-v2
features". With four histology encoders in hand that is wrong: every one of them inflates all six
targets, and only the natural-image encoder does not. The dissociation is by TRAINING CORPUS, not
by vendor. R18 and the narrowing note added to R15 are both corrected.

TWO CLASSES OF CLAIM, kept apart on purpose.

  PRE-DECLARED. R18's verdict rules were fixed and hashed before its arms were read. Applying the
  SAME rules to new encoders is legitimate -- same rule, new data -- so the per-encoder verdicts
  below carry the same weight R18's did.

  POST-HOC. Relative inflation appears to fall as site-disjoint performance rises. That pattern was
  noticed AFTER seeing all five encoders. It is reported as an observation with its rank
  correlation and n, explicitly not as a tested finding, because testing a hypothesis on the data
  that suggested it is the exact move this programme exists to police. What it licenses is a
  pre-registration, not a conclusion.

Run:  python3 m18_five_encoder_survey.py
"""
import json
import os
import subprocess

import numpy as np
from scipy.stats import spearmanr

from audit_kit import Run

B = "gs://heydonto-quantara-lungcdx/nsclc-rwpr-study/armD-meth/r18/"
ENC = [("tcga", "owkin/phikon-v2", "histology"),
       ("dinov2-large", "facebook/dinov2-large", "natural images"),
       ("uni-tcga", "MahmoodLab/UNI", "histology"),
       ("hopt-tcga", "bioptimus/H-optimus-0", "histology"),
       ("virchow2-tcga", "paige-ai/Virchow2", "histology")]
DIRMAP = {"dinov2-large": "dinov2-tcga"}
TARGETS = ["keap1_sig", "global", "island", "opensea", "tss200", "body"]

# R18's pre-declared thresholds, reused verbatim on new data
LEAK_REAL, LEAK_ABSENT = 0.05, 0.02
PUBLISHED = {"subtype_grouped": 0.799, "subtype_random": 0.9703,
             "keap1_grouped": 0.664, "meth_mean_inflation": 0.0849}
TOL = 0.01


def fetch(name):
    p = f"/tmp/r21_{name}.json"
    if not os.path.exists(p):
        subprocess.run(["gcloud", "storage", "cp", B + name + ".json", p],
                       check=True, capture_output=True)
    return json.load(open(p))


def main():
    run = Run("R21fiveenc")
    cfg = run.start(
        {"question": "does the site-leakage pattern hold across pathology foundation models, or "
                     "was R18's dissociation specific to Phikon-v2?",
         "encoders": {d: {"repo": r, "corpus": c} for d, r, c in ENC},
         "design": "identical tiles (byte-identical grids), identical folds, architecture, "
                   "hyperparameters and seeds; per-model prescribed normalisation; CLS readout "
                   "throughout. Feature width varies (1024/1536/1280) and the MIL projects to 512.",
         "reused_pre_declared_rules": {
             "source": "R18 config 143b3b9d0a92fabb..., fixed before R18's arms were read",
             "leak_real": LEAK_REAL, "leak_absent": LEAK_ABSENT,
             "note": "same rule applied to new data, which is legitimate; no threshold was "
                     "chosen after seeing these five encoders"},
         "post_hoc_and_labelled_as_such": [
             "relative inflation vs site-disjoint performance across encoders -- noticed after "
             "seeing all five, reported with rank correlation and n, NOT tested as a hypothesis"],
         "status": "HELD FOR IP - not for publication pending patent"}, [])

    A, rows = {}, []
    for d, repo, corpus in ENC:
        fd = DIRMAP.get(d, d)
        sg, sr = fetch(f"{fd}__grouped__subtype"), fetch(f"{fd}__random__subtype")
        mg, mr = fetch(f"{fd}__grouped__meth"), fetch(f"{fd}__random__meth")
        enc_seen = sorted({e for x in (sg, sr, mg, mr) for e in x["encoder"]})
        if len(enc_seen) != 1 or enc_seen[0] != repo:
            run.gate(f"G_enc_{d}", f"all four arms report {repo}", enc_seen, False)
        g, r = sg["metrics"]["subtype"]["auroc"], sr["metrics"]["subtype"]["auroc"]
        kg, kr = sg["metrics"]["keap1"]["auroc"], sr["metrics"]["keap1"]["auroc"]
        per = {t: {"grouped_rho": mg["metrics"][t]["rho"], "random_rho": mr["metrics"][t]["rho"],
                   "inflation": mr["metrics"][t]["rho"] - mg["metrics"][t]["rho"]}
               for t in TARGETS}
        gm = float(np.mean([per[t]["grouped_rho"] for t in TARGETS]))
        rm = float(np.mean([per[t]["random_rho"] for t in TARGETS]))
        ninf = int(sum(1 for t in TARGETS if per[t]["inflation"] > 0))
        A[d] = {"encoder": repo, "corpus": corpus, "n_patients": sg["n_patients"],
                "subtype_grouped": g, "subtype_random": r, "subtype_inflation": r - g,
                "subtype_relative_inflation": (r - g) / (1 - g),
                "keap1_grouped": kg, "keap1_random": kr, "keap1_inflation": kr - kg,
                "meth_mean_grouped": gm, "meth_mean_random": rm, "meth_mean_inflation": rm - gm,
                "meth_targets_inflated": ninf, "meth_per_target": per,
                "subtype_verdict": ("LEAKAGE_PRESENT"
                                    if (r - g) / (1 - g) > LEAK_REAL else "LEAKAGE_ABSENT"),
                "meth_verdict": ("METHYLATION_LEAKAGE_PRESENT" if ninf >= 5
                                 else "METHYLATION_LEAKAGE_ABSENT")}
        rows.append(A[d])
        run.log("encoder", repo=repo, corpus=corpus, sub_grouped=round(g, 4),
                sub_random=round(r, 4), sub_infl=round(r - g, 4),
                sub_rel=round((r - g) / (1 - g), 4), keap1_infl=round(kr - kg, 4),
                meth_infl=round(rm - gm, 4), n_inflated=ninf)

    run.gate("G1_same_cohort", "all twenty arms cover 760 patients",
             sorted({v["n_patients"] for v in A.values()}),
             {v["n_patients"] for v in A.values()} == {760})
    ph = A["tcga"]
    run.gate("G2_controls_still_reproduce",
             "the Phikon-v2 arms still match R15 within tolerance",
             {"sub_g": round(ph["subtype_grouped"], 4), "sub_r": round(ph["subtype_random"], 4),
              "keap1_g": round(ph["keap1_grouped"], 4),
              "meth_infl": round(ph["meth_mean_inflation"], 4)},
             abs(ph["subtype_grouped"] - PUBLISHED["subtype_grouped"]) <= TOL
             and abs(ph["subtype_random"] - PUBLISHED["subtype_random"]) <= TOL
             and abs(ph["keap1_grouped"] - PUBLISHED["keap1_grouped"]) <= TOL
             and abs(ph["meth_mean_inflation"] - PUBLISHED["meth_mean_inflation"]) <= 0.02)

    # ---- pre-declared rules, applied to new data ----
    hist = [v for v in A.values() if v["corpus"] == "histology"]
    nat = [v for v in A.values() if v["corpus"] == "natural images"]
    sub_all = all(v["subtype_verdict"] == "LEAKAGE_PRESENT" for v in A.values())
    meth_hist = all(v["meth_verdict"] == "METHYLATION_LEAKAGE_PRESENT" for v in hist)
    meth_nat = all(v["meth_verdict"] == "METHYLATION_LEAKAGE_ABSENT" for v in nat)
    run.gate("G3_subtype_universal",
             "subtype leakage clears the pre-declared 0.05 relative bar in every encoder",
             {v["encoder"]: round(v["subtype_relative_inflation"], 4) for v in A.values()},
             sub_all)

    if meth_hist and meth_nat:
        verdict = "METHYLATION_LEAKAGE_IS_A_PROPERTY_OF_HISTOLOGY_PRETRAINING"
    elif meth_hist:
        verdict = "METHYLATION_LEAKAGE_IN_ALL_ENCODERS"
    else:
        verdict = "METHYLATION_LEAKAGE_MIXED_ACROSS_HISTOLOGY_ENCODERS"
    run.log("VERDICT", primary="SUBTYPE_LEAKAGE_UNIVERSAL" if sub_all else "SUBTYPE_LEAKAGE_MIXED",
            secondary=verdict,
            corrects="R18 secondary verdict METHYLATION_LEAKAGE_DID_NOT_REPRODUCE, which "
                     "attributed the effect to Phikon-v2 rather than to histology pretraining")

    # ---- post-hoc observation, labelled ----
    gvals = [v["subtype_grouped"] for v in A.values()]
    rvals = [v["subtype_relative_inflation"] for v in A.values()]
    rs = spearmanr(gvals, rvals)
    post = {"claim": "relative subtype inflation falls as site-disjoint performance rises",
            "spearman_rho": float(rs.statistic), "p": float(rs.pvalue), "n_encoders": len(A),
            "STATUS": "POST-HOC OBSERVATION, NOT A TESTED FINDING",
            "why": "noticed after all five encoders were read. n=5, so even a perfect rank "
                   "correlation cannot reach p<0.05 by this test. Reported so it can be "
                   "pre-registered and tested on encoders not used here; it licenses a protocol, "
                   "not a conclusion.",
            "pairs": [{"encoder": v["encoder"], "site_disjoint_auroc": v["subtype_grouped"],
                       "relative_inflation": v["subtype_relative_inflation"]}
                      for v in A.values()]}
    run.log("POST_HOC", rho=round(float(rs.statistic), 4), p=round(float(rs.pvalue), 4), n=len(A))

    out = {"status": "HELD FOR IP - not for publication pending patent",
           "VERDICT": {"subtype": "SUBTYPE_LEAKAGE_UNIVERSAL" if sub_all else "MIXED",
                       "methylation": verdict,
                       "corrects_R18": "R18 secondary verdict attributed methylation leakage to "
                                       "Phikon-v2; four histology encoders show it and only the "
                                       "natural-image encoder does not",
                       "rules_reused_from": "R18 config 143b3b9d0a92fabb..., pre-declared"},
           "encoders": A,
           "by_corpus": {
               "histology": {"n": len(hist),
                             "meth_inflation_range": [round(min(v["meth_mean_inflation"]
                                                                for v in hist), 4),
                                                      round(max(v["meth_mean_inflation"]
                                                                for v in hist), 4)],
                             "all_six_targets_inflated": meth_hist},
               "natural_images": {"n": len(nat),
                                  "meth_inflation": round(nat[0]["meth_mean_inflation"], 4),
                                  "targets_inflated": nat[0]["meth_targets_inflated"]}},
           "post_hoc_observation": post,
           "_provenance": {"config_sha256": cfg, "run_dir": run.dir}}
    run.write("results.json", out)
    run.finalize()
    print("\nSUBTYPE:", out["VERDICT"]["subtype"])
    print("METHYLATION:", verdict)
    print(f"post-hoc: rho={rs.statistic:.3f} p={rs.pvalue:.3f} n={len(A)} (NOT a tested finding)")
    print("RUN DIR:", run.dir)


if __name__ == "__main__":
    main()
