# Claim manifest

Every load-bearing claim in both manuscripts, with the file and field it resolves to. Generated
from `verify_claims.py` so it cannot drift from the checker.

Run `python3 verify_claims.py` to verify all of them, including the negative control.

**Tolerance `0`** means the value must match exactly. Otherwise the absolute difference must fall
within the stated tolerance.


## paper1 — 51 claims

| # | Claim | Resolves to | Field | Expected | Tol |
|---|---|---|---|---|---|
| 1 | subtype AUROC, site-grouped folds | `reports/R15-wsi-methylation/evidence/results.json` | `HEADLINE_site_leakage.subtype.grouped` | 0.799 | 0.001 |
| 2 | subtype AUROC, random folds | `reports/R15-wsi-methylation/evidence/results.json` | `HEADLINE_site_leakage.subtype.random` | 0.9703 | 0.001 |
| 3 | subtype inflation from fold assignment alone | `reports/R15-wsi-methylation/evidence/results.json` | `HEADLINE_site_leakage.subtype.inflation` | 0.171 | 0.001 |
| 4 | KEAP1 AUROC, site-grouped | `reports/R15-wsi-methylation/evidence/results.json` | `HEADLINE_site_leakage.keap1.grouped` | 0.664 | 0.001 |
| 5 | mean methylation inflation across 6 targets | `reports/_shared/armD-meth/meth_leakage_arm.json` | `mean_inflation` | 0.0849 | 0.0001 |
| 6 | within-TCGA site probe, balanced accuracy | `reports/R12-federated-simulation/evidence/results.json` | `silo_signature_probe.balanced_accuracy` | 0.7062 | 0.0001 |
| 7 | within-TCGA site probe, permutation p | `reports/R12-federated-simulation/evidence/results.json` | `silo_signature_probe.permutation_p` | 0.0099 | 0.0001 |
| 8 | EAGLE vs HTMCP, published 78-slide set | `reports/R19-external-site-probe/evidence/results.json` | `arms.ARM_ALL.balanced_accuracy` | 1.0 | 0.0001 |
| 9 | EAGLE vs HTMCP, null mean, published set | `reports/R19-external-site-probe/evidence/results.json` | `arms.ARM_ALL.null_mean` | 0.5074 | 0.0001 |
| 10 | EAGLE vs HTMCP, H&E only | `reports/R19-external-site-probe/evidence/results.json` | `arms.ARM_HE.balanced_accuracy` | 1.0 | 0.0001 |
| 11 | EAGLE vs HTMCP, null mean, H&E only | `reports/R19-external-site-probe/evidence/results.json` | `arms.ARM_HE.null_mean` | 0.5058 | 0.0001 |
| 12 | three-way size-balanced probe, published set | `reports/R19-external-site-probe/evidence/results.json` | `arms.ARM3_ALL.balanced_accuracy` | 0.9816 | 0.0001 |
| 13 | three-way size-balanced probe, H&E only | `reports/R19-external-site-probe/evidence/results.json` | `arms.ARM3_HE.balanced_accuracy` | 0.9814 | 0.0001 |
| 14 | non-H&E slides in the published probe set | `reports/R19-external-site-probe/evidence/results.json` | `stain_contamination.n_non_he_in_published_set` | 3 | 0 |
| 15 | partial rho, KEAP1 signature, subtype only | `reports/R15-wsi-methylation/evidence/purity_control.json` | `targets.keap1_sig.partial_subtype_only` | 0.252 | 0.001 |
| 16 | partial rho, KEAP1 signature, subtype + purity | `reports/R15-wsi-methylation/evidence/purity_control.json` | `targets.keap1_sig.partial_subtype_plus_purity` | 0.221 | 0.001 |
| 17 | largest loss from purity adjustment (KEAP1 signature) | `reports/R15-wsi-methylation/evidence/purity_control.json` | `targets.keap1_sig.delta` | -0.031 | 0.001 |
| 18 | patients with ABSOLUTE purity available | `reports/R15-wsi-methylation/evidence/purity_control.json` | `n_with_purity` | 741 | 0 |
| 19 | genome-wide median rho, subtype-supervised (R15's negative) | `reports/R15-wsi-methylation/evidence/results.json` | `Q2_genomewide.median_rho` | 0.006 | 0.001 |
| 20 | CpGs above null, subtype-supervised | `reports/R15-wsi-methylation/evidence/results.json` | `Q2_genomewide.n_predictable` | 110212 | 0 |
| 21 | dinov2-large subtype AUROC, site-grouped | `reports/R18-encoder-robustness/evidence/results.json` | `subtype.dinov2-tcga.subtype_grouped` | 0.7356 | 0.001 |
| 22 | dinov2-large subtype AUROC, random folds | `reports/R18-encoder-robustness/evidence/results.json` | `subtype.dinov2-tcga.subtype_random` | 0.939 | 0.001 |
| 23 | dinov2-large subtype inflation (larger than Phikon-v2's) | `reports/R18-encoder-robustness/evidence/results.json` | `subtype.dinov2-tcga.subtype_inflation` | 0.2034 | 0.001 |
| 24 | dinov2-large relative subtype inflation | `reports/R18-encoder-robustness/evidence/results.json` | `subtype.dinov2-tcga.subtype_relative_inflation` | 0.769 | 0.001 |
| 25 | Phikon-v2 relative subtype inflation | `reports/R18-encoder-robustness/evidence/results.json` | `subtype.tcga.subtype_relative_inflation` | 0.852 | 0.001 |
| 26 | dinov2-large KEAP1 inflation (asymmetry reproduces) | `reports/R18-encoder-robustness/evidence/results.json` | `subtype.dinov2-tcga.keap1_inflation` | 0.0202 | 0.001 |
| 27 | dinov2-large mean methylation inflation (leakage did NOT reproduce) | `reports/R18-encoder-robustness/evidence/results.json` | `methylation.dinov2-tcga.mean_inflation` | 0.0004 | 0.0002 |
| 28 | Phikon-v2 mean methylation inflation | `reports/R18-encoder-robustness/evidence/results.json` | `methylation.tcga.mean_inflation` | 0.0849 | 0.0001 |
| 29 | Phikon-v2 external probe, all H&E in the archive | `reports/R20-external-probe-two-encoders/evidence/results.json` | `arms.phikon-v2.SET_HE.balanced_accuracy` | 1.0 | 0.0001 |
| 30 | Phikon-v2 external probe, null mean on that set | `reports/R20-external-probe-two-encoders/evidence/results.json` | `arms.phikon-v2.SET_HE.null_mean` | 0.5023 | 0.0001 |
| 31 | dinov2-large external probe, all H&E (Paper 1's qualification) | `reports/R20-external-probe-two-encoders/evidence/results.json` | `arms.dinov2-large.SET_HE.balanced_accuracy` | 0.8768 | 0.0001 |
| 32 | dinov2-large external probe, null mean | `reports/R20-external-probe-two-encoders/evidence/results.json` | `arms.dinov2-large.SET_HE.null_mean` | 0.499 | 0.0001 |
| 33 | HTMCP slides R19's input omitted | `reports/R20-external-probe-two-encoders/evidence/results.json` | `r19_subset_defect.htmcp_omitted_by_r19` | 51 | 0 |
| 34 | UNI subtype AUROC, site-grouped | `reports/R21-five-encoder-survey/evidence/results.json` | `encoders.uni-tcga.subtype_grouped` | 0.9328 | 0.001 |
| 35 | UNI subtype inflation | `reports/R21-five-encoder-survey/evidence/results.json` | `encoders.uni-tcga.subtype_inflation` | 0.0413 | 0.001 |
| 36 | UNI mean methylation inflation (histology encoder DOES leak) | `reports/R21-five-encoder-survey/evidence/results.json` | `encoders.uni-tcga.meth_mean_inflation` | 0.0606 | 0.001 |
| 37 | H-optimus-0 subtype AUROC, site-grouped | `reports/R21-five-encoder-survey/evidence/results.json` | `encoders.hopt-tcga.subtype_grouped` | 0.9211 | 0.001 |
| 38 | H-optimus-0 mean methylation inflation (smallest histology value) | `reports/R21-five-encoder-survey/evidence/results.json` | `encoders.hopt-tcga.meth_mean_inflation` | 0.0173 | 0.001 |
| 39 | Virchow2 subtype AUROC, site-grouped (best site-disjoint) | `reports/R21-five-encoder-survey/evidence/results.json` | `encoders.virchow2-tcga.subtype_grouped` | 0.9497 | 0.001 |
| 40 | Virchow2 relative subtype inflation (smallest of five) | `reports/R21-five-encoder-survey/evidence/results.json` | `encoders.virchow2-tcga.subtype_relative_inflation` | 0.444 | 0.001 |
| 41 | Virchow2 mean methylation inflation | `reports/R21-five-encoder-survey/evidence/results.json` | `encoders.virchow2-tcga.meth_mean_inflation` | 0.0773 | 0.001 |
| 42 | post-hoc capability-vs-inflation rank correlation (NOT a claim) | `reports/R21-five-encoder-survey/evidence/results.json` | `post_hoc_observation.spearman_rho` | -0.8 | 0.01 |
| 43 | post-hoc p, which is why it is not claimed | `reports/R21-five-encoder-survey/evidence/results.json` | `post_hoc_observation.p` | 0.104 | 0.001 |
| 44 | histology encoders inflating all six methylation targets | `reports/R21-five-encoder-survey/evidence/results.json` | `by_corpus.histology.n` | 4 | 0 |
| 45 | dinov2-large external, 25-split mean (replaces R20's single 0.8768) | `reports/R22-external-probe-five-encoders/evidence/results.json` | `encoders.facebook/dinov2-large.set_he_seed_distribution.mean` | 0.9176 | 0.001 |
| 46 | dinov2-large external, 25-split sd (why one split was not enough) | `reports/R22-external-probe-five-encoders/evidence/results.json` | `encoders.facebook/dinov2-large.set_he_seed_distribution.sd` | 0.0247 | 0.001 |
| 47 | Phikon-v2 external is genuinely saturated, sd exactly zero | `reports/R22-external-probe-five-encoders/evidence/results.json` | `encoders.owkin/phikon-v2.set_he_seed_distribution.sd` | 0.0 | 0.0001 |
| 48 | UNI external, 25-split mean | `reports/R22-external-probe-five-encoders/evidence/results.json` | `encoders.MahmoodLab/UNI.set_he_seed_distribution.mean` | 0.9937 | 0.001 |
| 49 | H-optimus-0 external, 25-split mean | `reports/R22-external-probe-five-encoders/evidence/results.json` | `encoders.bioptimus/H-optimus-0.set_he_seed_distribution.mean` | 0.9961 | 0.001 |
| 50 | Virchow2 external, 25-split mean | `reports/R22-external-probe-five-encoders/evidence/results.json` | `encoders.paige-ai/Virchow2.set_he_seed_distribution.mean` | 0.9601 | 0.001 |
| 51 | R20's single-split value reproduced exactly before being superseded | `reports/R22-external-probe-five-encoders/evidence/results.json` | `r20_reproduction.facebook/dinov2-large` | 0.8768 | 0.0001 |

## paper2 — 33 claims

| # | Claim | Resolves to | Field | Expected | Tol |
|---|---|---|---|---|---|
| 1 | PC1, PD-L1 TPS alone -> response (positive control) | `reports/R14-autonomous-generation/evidence/results.json` | `positive_controls.PC1.auroc` | 0.747 | 0.001 |
| 2 | PC1 permutation p | `reports/R14-autonomous-generation/evidence/results.json` | `positive_controls.PC1.perm_p` | 0.0001 | 2e-05 |
| 3 | PC1 n | `reports/R14-autonomous-generation/evidence/results.json` | `positive_controls.PC1.n` | 246 | 0 |
| 4 | multivariable baseline M_base, AUROC of mean OOF | `reports/R14-autonomous-generation/evidence/results.json` | `models.E1_M_base.auroc_of_mean_oof` | 0.745 | 0.001 |
| 5 | multivariable baseline M_base, mean of repeats | `reports/R14-autonomous-generation/evidence/results.json` | `models.E1_M_base.auroc_mean_of_repeats` | 0.742 | 0.001 |
| 6 | primary multimodal delta (H1) | `reports/R14-autonomous-generation/evidence/results.json` | `hypotheses.H1.delta` | -0.036 | 0.001 |
| 7 | primary comparison n (all modalities present) | `reports/R14-autonomous-generation/evidence/results.json` | `hypotheses.H1.n` | 73 | 0 |
| 8 | primary one-sided p | `reports/R14-autonomous-generation/evidence/results.json` | `hypotheses.H1.p_onesided` | 0.638 | 0.001 |
| 9 | primary hazard ratio per doubling of post-RT peak SUV | `reports/R16-acrin-arm-a/evidence/results.json` | `PRIMARY.hr_per_doubling` | 1.202 | 0.001 |
| 10 | primary permutation p | `reports/R16-acrin-arm-a/evidence/results.json` | `PRIMARY.permutation_p` | 0.0005 | 0.0001 |
| 11 | analysed patients | `reports/R16-acrin-arm-a/evidence/results.json` | `cohort.n_analysed` | 166 | 0 |
| 12 | events | `reports/R16-acrin-arm-a/evidence/results.json` | `cohort.events` | 125 | 0 |
| 13 | PC1 hazard ratio (progression by day 365) | `reports/R16-acrin-arm-a/evidence/results.json` | `PC1_detail.hr` | 2.463 | 0.001 |
| 14 | PC2 hazard ratio (performance status), the value a bug once reported as primary | `reports/R16-acrin-arm-a/evidence/results.json` | `PC2_detail.hr` | 1.382 | 0.001 |
| 15 | PC3 exposed n (uninformative by pre-declared rule) | `reports/R16-acrin-arm-a/evidence/results.json` | `PC3_detail.exposed_n` | 26 | 0 |
| 16 | exposure values floored at 0.5 | `reports/R16-acrin-arm-a/evidence/results.json` | `cohort.floored_at_0.5` | 106 | 0 |
| 17 | H-optimus-0 inflation that a magnitude bar would have failed | `reports/R21-five-encoder-survey/evidence/results.json` | `encoders.hopt-tcga.meth_mean_inflation` | 0.0173 | 0.001 |
| 18 | post-hoc rank correlation Paper 2 declines to claim | `reports/R21-five-encoder-survey/evidence/results.json` | `post_hoc_observation.spearman_rho` | -0.8 | 0.01 |
| 19 | its p, which is the reason it is not claimed | `reports/R21-five-encoder-survey/evidence/results.json` | `post_hoc_observation.p` | 0.104 | 0.001 |
| 20 | the seed-sensitive external value the gate halted on | `reports/R22-external-probe-five-encoders/evidence/results.json` | `encoders.facebook/dinov2-large.set_he_seed_distribution.mean` | 0.9176 | 0.001 |
| 21 | its across-split sd | `reports/R22-external-probe-five-encoders/evidence/results.json` | `encoders.facebook/dinov2-large.set_he_seed_distribution.sd` | 0.0247 | 0.001 |
| 22 | R23 pre-registered primary, the confounded one | `reports/R23-methylation-drug-resistance/evidence/results.json` | `VERDICT.rho` | 0.651 | 0.001 |
| 23 | R23 permutation p | `reports/R23-methylation-drug-resistance/evidence/results.json` | `VERDICT.perm_p` | 0.0005 | 0.0001 |
| 24 | R23 lines analysed | `reports/R23-methylation-drug-resistance/evidence/results.json` | `VERDICT.n_lines` | 187 | 0 |
| 25 | R23 post-hoc partial, after removing general drug sensitivity | `reports/R23-methylation-drug-resistance/evidence/results.json` | `H3_POSTHOC_general_sensitivity_control.partial_rho` | -0.16 | 0.001 |
| 26 | R23 control drugs used for that adjustment | `reports/R23-methylation-drug-resistance/evidence/results.json` | `H3_POSTHOC_general_sensitivity_control.n_control_drugs` | 227 | 0 |
| 27 | R24 SMARCA4 knockdown, the gated positive control | `reports/R24-smarca4-resistance-reversal/evidence/results.json` | `PC1_smarca4_knockdown.PC9OR_log2FC` | -2.37 | 0.01 |
| 28 | R24 resistance-UP genes fell below null as the thesis predicts | `reports/R24-smarca4-resistance-reversal/evidence/results.json` | `H2_PC9OR.resistance_up.mean_kd_lfc` | -0.0043 | 0.0005 |
| 29 | R24 resistance-DOWN genes also fell, which is why the verdict is PARTIAL | `reports/R24-smarca4-resistance-reversal/evidence/results.json` | `H2_PC9OR.resistance_down.mean_kd_lfc` | -0.0795 | 0.0005 |
| 30 | R24 replication inverts the sign in YU005C | `reports/R24-smarca4-resistance-reversal/evidence/results.json` | `H3_YU005C_replication.resistance_up.mean_kd_lfc` | 0.0284 | 0.0005 |
| 31 | methylation-supervised CpGs above null | `reports/R17-percpg-methylation-supervised/evidence/results.json` | `arms.meth.n_above_null` | 83369 | 0 |
| 32 | subtype-supervised CpGs above null | `reports/R17-percpg-methylation-supervised/evidence/results.json` | `arms.sub.n_above_null` | 90137 | 0 |
| 33 | observed-global-mean baseline CpGs above null | `reports/R17-percpg-methylation-supervised/evidence/results.json` | `arms.globalmean.n_above_null` | 298320 | 0 |
