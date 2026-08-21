#!/usr/bin/env python3
"""Regenerate the R10 figure from evidence/ only (no matrix access needed)."""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

r = json.load(open("evidence/results.json"))
a = json.load(open("evidence/audit_addendum.json"))
NAVY, GREY, RED, BLUE = "#14284b", "#5a5a5a", "#961e1e", "#c8d4e6"

fig, ax = plt.subplots(2, 2, figsize=(11.6, 8.0))

# ---- A: how strong is the phenotype ----
q1, pc = r["Q1_keap1_classifier"], r["PC_sex_classifier"]
C, D = a["C_no_smoking_probes"], a["D_seed_stability"]
labels = ["sex\n(positive control)", "KEAP1 status", "KEAP1, smoking\nprobes deleted"]
vals = [pc["auroc"], q1["auroc"], C["auroc_smoking_removed"]]
errs = [[pc["auroc"] - pc["ci95"][0]], [q1["auroc"] - q1["ci95"][0]], [0]]
hi = [[pc["ci95"][1] - pc["auroc"]], [q1["ci95"][1] - q1["auroc"]], [0]]
A = ax[0, 0]
A.bar(range(3), vals, 0.55, color=[BLUE, NAVY, NAVY], edgecolor=GREY, lw=.5)
for k in range(2):
    A.errorbar([k], [vals[k]], yerr=[errs[k], hi[k]], fmt="none", ecolor=GREY,
               capsize=4, lw=1.4)
A.axhspan(D["range"][0], D["range"][1], color=RED, alpha=.13, zorder=0)
A.text(-0.42, (D["range"][0] + D["range"][1]) / 2,
       f"seed range\n{D['range'][0]:.3f}–{D['range'][1]:.3f}",
       fontsize=7, color=RED, ha="left", va="center")
A.axhline(0.5, color=GREY, ls=":", lw=1.2)
tops = [pc["ci95"][1], q1["ci95"][1], C["auroc_smoking_removed"]]
for k, v in enumerate(vals):
    A.text(k, tops[k] + .018, f"{v:.3f}", ha="center", fontsize=9.5, color=NAVY)
A.set_xticks(range(3)); A.set_xticklabels(labels, fontsize=8.5)
A.set_ylim(0.45, 1.08); A.set_ylabel("AUROC (nested CV)")
A.set_title("A  The phenotype is strong, and it is not smoking\n"
            f"KEAP1 predictable at {q1['auroc']:.3f} from methylation alone",
            fontsize=10, color=NAVY, loc="left")

# ---- B: where in the genome ----
ctx = r["Q2_genomic_context"]["cpg_island_context"]
order = ["Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"]
B = ax[0, 1]
fr = [ctx[c]["frac_sig"] * 100 for c in order]
cols = [RED if ctx[c]["odds_ratio"] < 1 else NAVY for c in order]
B.bar(range(len(order)), fr, 0.6, color=cols, edgecolor=GREY, lw=.5)
overall = r["Q2_genomic_context"]["overall_frac_sig"] * 100
B.axhline(overall, color=GREY, ls="--", lw=1.2)
B.text(-0.35, overall + .7, f"array-wide {overall:.1f}%", fontsize=7.5, color=GREY,
       ha="left")
for k, c in enumerate(order):
    B.text(k, fr[k] + .5, f"{ctx[c]['odds_ratio']:.2f}", ha="center", fontsize=7.5,
           color=cols[k])
B.set_xticks(range(len(order)))
B.set_xticklabels(["Island", "N shore", "S shore", "N shelf", "S shelf", "Open sea"],
                  fontsize=8, rotation=20)
B.set_ylabel("% of probes differential")
B.set_ylim(0, 40)
B.set_title("B  It spares promoter CpG islands\nand concentrates in open sea "
            "(odds ratios shown)", fontsize=10, color=NAVY, loc="left")

# ---- C: NRF2 enrichment vs matched null ----
C_ = ax[1, 0]
obs = r["Q3_nrf2_enrichment"]["nrf2_mean_frac_sig"]
mn, mp = a["A_nrf2"]["matched_null_mean"], a["A_nrf2"]["matched_p"]
rng = np.random.default_rng(0)
# the archived null summary is a mean; draw the reported distribution shape around it
C_.axvline(obs, color=RED, lw=2.2)
C_.axvline(mn, color=NAVY, ls="--", lw=1.6)
C_.axvline(r["Q3_nrf2_enrichment"]["random_null_p95"], color=GREY, ls=":", lw=1.4)
C_.set_xlim(0.205, 0.405)
C_.set_ylim(0, 1)
C_.set_yticks([])
C_.annotate(f"NRF2 gene set\n{obs:.3f}", xy=(obs, .82), xytext=(obs + .012, .82),
            fontsize=9, color=RED, va="center")
C_.annotate(f"probe-count-matched\nnull mean {mn:.3f}", xy=(mn, .55),
            xytext=(mn + .009, .55), fontsize=8, color=NAVY, va="center", ha="left")
C_.annotate(f"null 95th pct\n{r['Q3_nrf2_enrichment']['random_null_p95']:.3f}",
            xy=(r["Q3_nrf2_enrichment"]["random_null_p95"], .16),
            xytext=(r["Q3_nrf2_enrichment"]["random_null_p95"] + .008, .18),
            fontsize=7.5, color=GREY, va="center", ha="left")
C_.set_xlabel("mean fraction of probes differential, per gene")
C_.set_title(f"C  The pre-declared NRF2 pathway is enriched\np = {mp:.3f} against "
             f"1,000 size-matched sets ({a['A_nrf2']['n_nrf2_genes']} genes)",
             fontsize=10, color=NAVY, loc="left")

# ---- D: overlap with the smoking signature ----
D_ = ax[1, 1]
B4 = a["B_overlap"]
bars = [("KEAP1 signature", B4["keap1_size"], NAVY),
        ("smoking signature\n(in KEAP1-WT only)", B4["smoking_size"], BLUE),
        ("observed overlap", B4["observed"], RED),
        ("overlap expected\nif independent", B4["expected_if_independent"], GREY)]
D_.barh(range(len(bars))[::-1], [b[1] for b in bars], 0.6,
        color=[b[2] for b in bars], edgecolor=GREY, lw=.5)
for k, b in enumerate(bars):
    D_.text(b[1] + 3000, list(range(len(bars)))[::-1][k], f"{b[1]:,.0f}", va="center",
            fontsize=8.5, color=NAVY)
D_.set_yticks(range(len(bars))[::-1])
D_.set_yticklabels([b[0] for b in bars], fontsize=8)
D_.set_xlim(0, 150000)
D_.set_xlabel("probes")
D_.set_title(f"D  Smoking explains {B4['frac_of_keap1']*100:.1f}% of it, "
             f"only {B4['obs_over_exp']:.2f}× chance\n"
             "R08 called this confound unresolvable", fontsize=10, color=NAVY, loc="left")

for row in ax:
    for a_ in row:
        a_.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("figures/r10_panels.pdf")
plt.savefig("figures/r10_panels.png", dpi=150)
print("figures/r10_panels.{pdf,png}")
