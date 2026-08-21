#!/usr/bin/env python3
"""Resolve each load-bearing claim in both manuscripts to a NAMED file and field.

WHY THIS IS THE SECOND VERSION. R19 found that Paper 1 quoted two numbers with no persisted
artefact, and that the earlier verification pass had called them verified because they matched what
I believed. The first version of this script tried to fix that by searching every evidence file for
each number. Its negative control killed it:

    fraction of FABRICATED numbers wrongly reported as grounded
      4-dp values in [0,1] (AUROC / rho / balanced accuracy)   100.0%
      3-dp values in [0,1]                                     100.0%
      hazard-ratio-like 3-dp values in [1,3]                   100.0%
      6-digit counts                                            93.7%
      7-digit counts                                            21.6%

With 1.5 million distinct numeric strings across 1,356 artefacts, the haystack is saturated: every
possible AUROC to four decimals appears somewhere. A test that cannot fail is not a test, and
reporting its "0 ungrounded" as reassurance would have been the same error R19 uncovered, committed
twice.

So verification is scoped instead. Each claim names one file and one field. The value must match
there, to the precision the manuscript states. A file that does not exist, a field that does not
exist, and a value that disagrees are three distinct failures and are reported separately.

The negative control is part of the run, not an afterthought: every claim is re-checked with its
expected value perturbed, and the suite must reject all of those. If it does not, the run fails.

Run:  python3 papers/verify_claims.py
"""
import json
import os
import sys

import pathlib
ROOT = str(pathlib.Path(__file__).resolve().parent)
RES = f"{ROOT}/reports"
ARM = f"{ROOT}/reports/_shared"

# (paper, description, file, dotted-path, expected, tolerance)
# tolerance 0 means "must match at the stated precision after rounding to its decimals"
CLAIMS = [
    # ---------------- Paper 1: site leakage, from R15 ----------------
    ("paper1", "subtype AUROC, site-grouped folds",
     f"{RES}/R15-wsi-methylation/evidence/results.json",
     "HEADLINE_site_leakage.subtype.grouped", 0.799, 0.001),
    ("paper1", "subtype AUROC, random folds",
     f"{RES}/R15-wsi-methylation/evidence/results.json",
     "HEADLINE_site_leakage.subtype.random", 0.9703, 0.001),
    ("paper1", "subtype inflation from fold assignment alone",
     f"{RES}/R15-wsi-methylation/evidence/results.json",
     "HEADLINE_site_leakage.subtype.inflation", 0.171, 0.001),
    ("paper1", "KEAP1 AUROC, site-grouped",
     f"{RES}/R15-wsi-methylation/evidence/results.json",
     "HEADLINE_site_leakage.keap1.grouped", 0.664, 0.001),
    ("paper1", "mean methylation inflation across 6 targets",
     f"{ARM}/armD-meth/meth_leakage_arm.json", "mean_inflation", 0.0849, 0.0001),

    # ---------------- Paper 1: site probes ----------------
    ("paper1", "within-TCGA site probe, balanced accuracy",
     f"{RES}/R12-federated-simulation/evidence/results.json",
     "silo_signature_probe.balanced_accuracy", 0.7062, 0.0001),
    ("paper1", "within-TCGA site probe, permutation p",
     f"{RES}/R12-federated-simulation/evidence/results.json",
     "silo_signature_probe.permutation_p", 0.0099, 0.0001),
    ("paper1", "EAGLE vs HTMCP, published 78-slide set",
     f"{RES}/R19-external-site-probe/evidence/results.json",
     "arms.ARM_ALL.balanced_accuracy", 1.0000, 0.0001),
    ("paper1", "EAGLE vs HTMCP, null mean, published set",
     f"{RES}/R19-external-site-probe/evidence/results.json",
     "arms.ARM_ALL.null_mean", 0.5074, 0.0001),
    ("paper1", "EAGLE vs HTMCP, H&E only",
     f"{RES}/R19-external-site-probe/evidence/results.json",
     "arms.ARM_HE.balanced_accuracy", 1.0000, 0.0001),
    ("paper1", "EAGLE vs HTMCP, null mean, H&E only",
     f"{RES}/R19-external-site-probe/evidence/results.json",
     "arms.ARM_HE.null_mean", 0.5058, 0.0001),
    ("paper1", "three-way size-balanced probe, published set",
     f"{RES}/R19-external-site-probe/evidence/results.json",
     "arms.ARM3_ALL.balanced_accuracy", 0.9816, 0.0001),
    ("paper1", "three-way size-balanced probe, H&E only",
     f"{RES}/R19-external-site-probe/evidence/results.json",
     "arms.ARM3_HE.balanced_accuracy", 0.9814, 0.0001),
    ("paper1", "non-H&E slides in the published probe set",
     f"{RES}/R19-external-site-probe/evidence/results.json",
     "stain_contamination.n_non_he_in_published_set", 3, 0),

    # ---------------- Paper 1: methylation from morphology, R15 ----------------
    # NOTE: the purity control lives in its own artefact, not in results.json. A first version of
    # this manifest guessed results.json and reported NO_FIELD -- which is the manifest being
    # wrong, not the evidence being absent. Recorded because a NO_FIELD must be chased to ground
    # before it is called a missing artefact.
    ("paper1", "partial rho, KEAP1 signature, subtype only",
     f"{RES}/R15-wsi-methylation/evidence/purity_control.json",
     "targets.keap1_sig.partial_subtype_only", 0.252, 0.001),
    ("paper1", "partial rho, KEAP1 signature, subtype + purity",
     f"{RES}/R15-wsi-methylation/evidence/purity_control.json",
     "targets.keap1_sig.partial_subtype_plus_purity", 0.221, 0.001),
    ("paper1", "largest loss from purity adjustment (KEAP1 signature)",
     f"{RES}/R15-wsi-methylation/evidence/purity_control.json",
     "targets.keap1_sig.delta", -0.031, 0.001),
    ("paper1", "patients with ABSOLUTE purity available",
     f"{RES}/R15-wsi-methylation/evidence/purity_control.json", "n_with_purity", 741, 0),
    ("paper1", "genome-wide median rho, subtype-supervised (R15's negative)",
     f"{RES}/R15-wsi-methylation/evidence/results.json", "Q2_genomewide.median_rho", 0.006, 0.001),
    ("paper1", "CpGs above null, subtype-supervised",
     f"{RES}/R15-wsi-methylation/evidence/results.json", "Q2_genomewide.n_predictable", 110212, 0),

    # ---------------- Paper 2: Vanguri cohort, R14 ----------------
    # PC1 is 'PD-L1 TPS raw score -> response' -- a SINGLE-variable positive control. Paper 2
    # originally attributed its 0.747 to the multivariable clinical baseline; corrected 2026-08-17.
    ("paper2", "PC1, PD-L1 TPS alone -> response (positive control)",
     f"{RES}/R14-autonomous-generation/evidence/results.json",
     "positive_controls.PC1.auroc", 0.747, 0.001),
    ("paper2", "PC1 permutation p",
     f"{RES}/R14-autonomous-generation/evidence/results.json",
     "positive_controls.PC1.perm_p", 0.0001, 0.00002),
    ("paper2", "PC1 n",
     f"{RES}/R14-autonomous-generation/evidence/results.json", "positive_controls.PC1.n", 246, 0),
    ("paper2", "multivariable baseline M_base, AUROC of mean OOF",
     f"{RES}/R14-autonomous-generation/evidence/results.json",
     "models.E1_M_base.auroc_of_mean_oof", 0.745, 0.001),
    ("paper2", "multivariable baseline M_base, mean of repeats",
     f"{RES}/R14-autonomous-generation/evidence/results.json",
     "models.E1_M_base.auroc_mean_of_repeats", 0.742, 0.001),
    ("paper2", "primary multimodal delta (H1)",
     f"{RES}/R14-autonomous-generation/evidence/results.json", "hypotheses.H1.delta", -0.036, 0.001),
    ("paper2", "primary comparison n (all modalities present)",
     f"{RES}/R14-autonomous-generation/evidence/results.json", "hypotheses.H1.n", 73, 0),
    ("paper2", "primary one-sided p",
     f"{RES}/R14-autonomous-generation/evidence/results.json",
     "hypotheses.H1.p_onesided", 0.638, 0.001),

    # ---------------- Paper 2: ACRIN 6668, R16 ----------------
    ("paper2", "primary hazard ratio per doubling of post-RT peak SUV",
     f"{RES}/R16-acrin-arm-a/evidence/results.json", "PRIMARY.hr_per_doubling", 1.202, 0.001),
    ("paper2", "primary permutation p",
     f"{RES}/R16-acrin-arm-a/evidence/results.json", "PRIMARY.permutation_p", 0.0005, 0.0001),
    ("paper2", "analysed patients",
     f"{RES}/R16-acrin-arm-a/evidence/results.json", "cohort.n_analysed", 166, 0),
    ("paper2", "events",
     f"{RES}/R16-acrin-arm-a/evidence/results.json", "cohort.events", 125, 0),
    ("paper2", "PC1 hazard ratio (progression by day 365)",
     f"{RES}/R16-acrin-arm-a/evidence/results.json", "PC1_detail.hr", 2.463, 0.001),
    ("paper2", "PC2 hazard ratio (performance status), the value a bug once reported as primary",
     f"{RES}/R16-acrin-arm-a/evidence/results.json", "PC2_detail.hr", 1.382, 0.001),
    ("paper2", "PC3 exposed n (uninformative by pre-declared rule)",
     f"{RES}/R16-acrin-arm-a/evidence/results.json", "PC3_detail.exposed_n", 26, 0),
    ("paper2", "exposure values floored at 0.5",
     f"{RES}/R16-acrin-arm-a/evidence/results.json", "cohort.floored_at_0.5", 106, 0),

    # ---------------- Paper 1: R18, the encoder question ----------------
    ("paper1", "dinov2-large subtype AUROC, site-grouped",
     f"{RES}/R18-encoder-robustness/evidence/results.json",
     "subtype.dinov2-tcga.subtype_grouped", 0.7356, 0.001),
    ("paper1", "dinov2-large subtype AUROC, random folds",
     f"{RES}/R18-encoder-robustness/evidence/results.json",
     "subtype.dinov2-tcga.subtype_random", 0.9390, 0.001),
    ("paper1", "dinov2-large subtype inflation (larger than Phikon-v2's)",
     f"{RES}/R18-encoder-robustness/evidence/results.json",
     "subtype.dinov2-tcga.subtype_inflation", 0.2034, 0.001),
    ("paper1", "dinov2-large relative subtype inflation",
     f"{RES}/R18-encoder-robustness/evidence/results.json",
     "subtype.dinov2-tcga.subtype_relative_inflation", 0.769, 0.001),
    ("paper1", "Phikon-v2 relative subtype inflation",
     f"{RES}/R18-encoder-robustness/evidence/results.json",
     "subtype.tcga.subtype_relative_inflation", 0.852, 0.001),
    ("paper1", "dinov2-large KEAP1 inflation (asymmetry reproduces)",
     f"{RES}/R18-encoder-robustness/evidence/results.json",
     "subtype.dinov2-tcga.keap1_inflation", 0.0202, 0.001),
    ("paper1", "dinov2-large mean methylation inflation (leakage did NOT reproduce)",
     f"{RES}/R18-encoder-robustness/evidence/results.json",
     "methylation.dinov2-tcga.mean_inflation", 0.0004, 0.0002),
    ("paper1", "Phikon-v2 mean methylation inflation",
     f"{RES}/R18-encoder-robustness/evidence/results.json",
     "methylation.tcga.mean_inflation", 0.0849, 0.0001),

    # ---------------- Paper 1: R20, the external probe across encoders ----------------
    ("paper1", "Phikon-v2 external probe, all H&E in the archive",
     f"{RES}/R20-external-probe-two-encoders/evidence/results.json",
     "arms.phikon-v2.SET_HE.balanced_accuracy", 1.0000, 0.0001),
    ("paper1", "Phikon-v2 external probe, null mean on that set",
     f"{RES}/R20-external-probe-two-encoders/evidence/results.json",
     "arms.phikon-v2.SET_HE.null_mean", 0.5023, 0.0001),
    ("paper1", "dinov2-large external probe, all H&E (Paper 1's qualification)",
     f"{RES}/R20-external-probe-two-encoders/evidence/results.json",
     "arms.dinov2-large.SET_HE.balanced_accuracy", 0.8768, 0.0001),
    ("paper1", "dinov2-large external probe, null mean",
     f"{RES}/R20-external-probe-two-encoders/evidence/results.json",
     "arms.dinov2-large.SET_HE.null_mean", 0.4990, 0.0001),
    ("paper1", "HTMCP slides R19's input omitted",
     f"{RES}/R20-external-probe-two-encoders/evidence/results.json",
     "r19_subset_defect.htmcp_omitted_by_r19", 51, 0),

    # ---------------- Paper 1: R21, the five-encoder survey ----------------
    ("paper1", "UNI subtype AUROC, site-grouped",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "encoders.uni-tcga.subtype_grouped", 0.9328, 0.001),
    ("paper1", "UNI subtype inflation",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "encoders.uni-tcga.subtype_inflation", 0.0413, 0.001),
    ("paper1", "UNI mean methylation inflation (histology encoder DOES leak)",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "encoders.uni-tcga.meth_mean_inflation", 0.0606, 0.001),
    ("paper1", "H-optimus-0 subtype AUROC, site-grouped",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "encoders.hopt-tcga.subtype_grouped", 0.9211, 0.001),
    ("paper1", "H-optimus-0 mean methylation inflation (smallest histology value)",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "encoders.hopt-tcga.meth_mean_inflation", 0.0173, 0.001),
    ("paper1", "Virchow2 subtype AUROC, site-grouped (best site-disjoint)",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "encoders.virchow2-tcga.subtype_grouped", 0.9497, 0.001),
    ("paper1", "Virchow2 relative subtype inflation (smallest of five)",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "encoders.virchow2-tcga.subtype_relative_inflation", 0.444, 0.001),
    ("paper1", "Virchow2 mean methylation inflation",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "encoders.virchow2-tcga.meth_mean_inflation", 0.0773, 0.001),
    ("paper1", "post-hoc capability-vs-inflation rank correlation (NOT a claim)",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "post_hoc_observation.spearman_rho", -0.80, 0.01),
    ("paper1", "post-hoc p, which is why it is not claimed",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "post_hoc_observation.p", 0.104, 0.001),
    ("paper1", "histology encoders inflating all six methylation targets",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "by_corpus.histology.n", 4, 0),

    # ---------------- Paper 1: R22, the external probe over 25 CV splits ----------------
    # R20 quoted single-split values; these are the distributions Paper 1 now cites.
    ("paper1", "dinov2-large external, 25-split mean (replaces R20's single 0.8768)",
     f"{RES}/R22-external-probe-five-encoders/evidence/results.json",
     "encoders.facebook/dinov2-large.set_he_seed_distribution.mean", 0.9176, 0.001),
    ("paper1", "dinov2-large external, 25-split sd (why one split was not enough)",
     f"{RES}/R22-external-probe-five-encoders/evidence/results.json",
     "encoders.facebook/dinov2-large.set_he_seed_distribution.sd", 0.0247, 0.001),
    ("paper1", "Phikon-v2 external is genuinely saturated, sd exactly zero",
     f"{RES}/R22-external-probe-five-encoders/evidence/results.json",
     "encoders.owkin/phikon-v2.set_he_seed_distribution.sd", 0.0, 0.0001),
    ("paper1", "UNI external, 25-split mean",
     f"{RES}/R22-external-probe-five-encoders/evidence/results.json",
     "encoders.MahmoodLab/UNI.set_he_seed_distribution.mean", 0.9937, 0.001),
    ("paper1", "H-optimus-0 external, 25-split mean",
     f"{RES}/R22-external-probe-five-encoders/evidence/results.json",
     "encoders.bioptimus/H-optimus-0.set_he_seed_distribution.mean", 0.9961, 0.001),
    ("paper1", "Virchow2 external, 25-split mean",
     f"{RES}/R22-external-probe-five-encoders/evidence/results.json",
     "encoders.paige-ai/Virchow2.set_he_seed_distribution.mean", 0.9601, 0.001),
    ("paper1", "R20's single-split value reproduced exactly before being superseded",
     f"{RES}/R22-external-probe-five-encoders/evidence/results.json",
     "r20_reproduction.facebook/dinov2-large", 0.8768, 0.0001),

    # ---------------- Paper 2: the two runs where the discipline bound ----------------
    ("paper2", "H-optimus-0 inflation that a magnitude bar would have failed",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "encoders.hopt-tcga.meth_mean_inflation", 0.0173, 0.001),
    ("paper2", "post-hoc rank correlation Paper 2 declines to claim",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "post_hoc_observation.spearman_rho", -0.80, 0.01),
    ("paper2", "its p, which is the reason it is not claimed",
     f"{RES}/R21-five-encoder-survey/evidence/results.json",
     "post_hoc_observation.p", 0.104, 0.001),
    ("paper2", "the seed-sensitive external value the gate halted on",
     f"{RES}/R22-external-probe-five-encoders/evidence/results.json",
     "encoders.facebook/dinov2-large.set_he_seed_distribution.mean", 0.9176, 0.001),
    ("paper2", "its across-split sd",
     f"{RES}/R22-external-probe-five-encoders/evidence/results.json",
     "encoders.facebook/dinov2-large.set_he_seed_distribution.sd", 0.0247, 0.001),

    # ---------------- Paper 2: R23/R24, the two limitation cases ----------------
    ("paper2", "R23 pre-registered primary, the confounded one",
     f"{RES}/R23-methylation-drug-resistance/evidence/results.json", "VERDICT.rho", 0.651, 0.001),
    ("paper2", "R23 permutation p",
     f"{RES}/R23-methylation-drug-resistance/evidence/results.json",
     "VERDICT.perm_p", 0.0005, 0.0001),
    ("paper2", "R23 lines analysed",
     f"{RES}/R23-methylation-drug-resistance/evidence/results.json", "VERDICT.n_lines", 187, 0),
    ("paper2", "R23 post-hoc partial, after removing general drug sensitivity",
     f"{RES}/R23-methylation-drug-resistance/evidence/results.json",
     "H3_POSTHOC_general_sensitivity_control.partial_rho", -0.160, 0.001),
    ("paper2", "R23 control drugs used for that adjustment",
     f"{RES}/R23-methylation-drug-resistance/evidence/results.json",
     "H3_POSTHOC_general_sensitivity_control.n_control_drugs", 227, 0),
    ("paper2", "R24 SMARCA4 knockdown, the gated positive control",
     f"{RES}/R24-smarca4-resistance-reversal/evidence/results.json",
     "PC1_smarca4_knockdown.PC9OR_log2FC", -2.37, 0.01),
    ("paper2", "R24 resistance-UP genes fell below null as the thesis predicts",
     f"{RES}/R24-smarca4-resistance-reversal/evidence/results.json",
     "H2_PC9OR.resistance_up.mean_kd_lfc", -0.0043, 0.0005),
    ("paper2", "R24 resistance-DOWN genes also fell, which is why the verdict is PARTIAL",
     f"{RES}/R24-smarca4-resistance-reversal/evidence/results.json",
     "H2_PC9OR.resistance_down.mean_kd_lfc", -0.0795, 0.0005),
    ("paper2", "R24 replication inverts the sign in YU005C",
     f"{RES}/R24-smarca4-resistance-reversal/evidence/results.json",
     "H3_YU005C_replication.resistance_up.mean_kd_lfc", 0.0284, 0.0005),

    # ---------------- Paper 2: R17, the supervision correction ----------------
    ("paper2", "methylation-supervised CpGs above null",
     f"{RES}/R17-percpg-methylation-supervised/evidence/results.json",
     "arms.meth.n_above_null", 83369, 0),
    ("paper2", "subtype-supervised CpGs above null",
     f"{RES}/R17-percpg-methylation-supervised/evidence/results.json",
     "arms.sub.n_above_null", 90137, 0),
    ("paper2", "observed-global-mean baseline CpGs above null",
     f"{RES}/R17-percpg-methylation-supervised/evidence/results.json",
     "arms.globalmean.n_above_null", 298320, 0),
]


def dig(obj, path):
    """Follow a dotted path; keys may themselves contain dots (e.g. floored_at_0.5)."""
    cur = obj
    parts = path.split(".")
    i = 0
    while i < len(parts):
        if not isinstance(cur, dict):
            return None, f"path stops at a non-object before '{parts[i]}'"
        for j in range(len(parts), i, -1):          # greedy: try longest key first
            k = ".".join(parts[i:j])
            if k in cur:
                cur = cur[k]
                i = j
                break
        else:
            return None, f"field '{parts[i]}' not present"
    return cur, None


def check(path, field, expected, tol):
    if not os.path.exists(path):
        return "NO_FILE", f"{os.path.relpath(path, ROOT)} does not exist"
    try:
        obj = json.load(open(path))
    except Exception as e:
        return "BAD_FILE", f"unreadable: {e}"
    got, err = dig(obj, field)
    if err:
        return "NO_FIELD", err
    if isinstance(got, bool) or not isinstance(got, (int, float)):
        return "NOT_NUMERIC", f"field holds {type(got).__name__}: {str(got)[:40]}"
    if tol == 0:
        ok = float(got) == float(expected)
    else:
        ok = abs(float(got) - float(expected)) <= tol
    return ("OK", f"{got}") if ok else ("MISMATCH", f"file says {got}, manuscript says {expected}")


def main():
    print(f"scoped verification of {len(CLAIMS)} load-bearing claims\n")
    fails, bypaper = [], {}
    for paper, desc, path, field, exp, tol in CLAIMS:
        status, detail = check(path, field, exp, tol)
        bypaper.setdefault(paper, []).append(status)
        mark = "ok  " if status == "OK" else "FAIL"
        if status != "OK":
            fails.append((paper, desc, status, detail))
        print(f"  [{mark}] {paper}  {desc}")
        print(f"         {os.path.relpath(path, ROOT)} :: {field}")
        print(f"         expect {exp}  ->  {status}: {detail}")

    print("\n--- negative control: every claim re-checked with a perturbed expectation ---")
    leaks = []
    for paper, desc, path, field, exp, tol in CLAIMS:
        bad = (float(exp) + max(tol * 10, 0.05)) if float(exp) != 0 else 0.5
        if isinstance(exp, int) and tol == 0:
            bad = int(exp) + 7
        st, _ = check(path, field, bad, tol)
        if st == "OK":
            leaks.append((paper, desc))
    print(f"  perturbed claims wrongly accepted: {len(leaks)} of {len(CLAIMS)}")
    for p, d in leaks:
        print(f"    LEAK {p} {d}")
    power = 100.0 * (1 - len(leaks) / len(CLAIMS))
    print(f"  detection power: {power:.1f}%")

    print("\n--- summary ---")
    for p, sts in sorted(bypaper.items()):
        print(f"  {p}: {sts.count('OK')}/{len(sts)} resolved to a named file and field")
    if fails:
        print(f"\n  {len(fails)} FAILING claims:")
        for p, d, s, det in fails:
            print(f"    {p}  {d}\n       {s}: {det}")
    ok = not fails and not leaks
    print(f"\n  RESULT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
