#!/usr/bin/env python3
"""Regenerate the R08 figure. Panel A reads the beta matrix; B and C read only evidence JSON."""
import csv, json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "evidence/"
res = json.load(open(E + "results.json"))
aud = json.load(open(E + "audit_addendum.json"))
NAVY, GREY, RED, BLUE = "#14284b", "#5a5a5a", "#961e1e", "#c8d4e6"
LOCAL = "/Users/rezanehzati/quantara-staging/r08"

fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0))

# ---- A: the AHRR smoking control, the one R07 could not power ----
pc3 = res["PC3_AHRR_smoking"]
if os.path.exists(f"{LOCAL}/beta_450k.npy"):
    probes = open(f"{LOCAL}/probes.txt").read().split()
    i = probes.index("cg05575921")
    M = np.load(f"{LOCAL}/beta_450k.npy", mmap_mode="r")
    samples = list(csv.DictReader(open(f"{LOCAL}/samples.tsv"), delimiter="\t"))
    by = {}
    for s in samples:
        if s["sample"] not in by or s["file_id"] < by[s["sample"]]["file_id"]:
            by[s["sample"]] = s
    cols = sorted(int(s["col"]) for s in by.values())
    smeta = [samples[c] for c in cols]
    row = np.asarray(M[i, cols], np.float64)
    clin = {r["case"]: r for r in csv.DictReader(open(f"{LOCAL}/clinical_merged.tsv"),
                                                delimiter="\t")}
    tum = np.array([s["sample_type"] != "Solid Tissue Normal" for s in smeta])
    sm = np.array([clin.get(s["case"], {}).get("smoke_status", "") for s in smeta], object)
    cur = row[tum & (sm == "Current Smoker")]
    nev = row[tum & np.array(["Non-Smoker" in x for x in sm])]
    cur, nev = cur[~np.isnan(cur)], nev[~np.isnan(nev)]
    bp = ax[0].boxplot([nev, cur], tick_labels=[f"Never\n(n={len(nev)})",
                                               f"Current smoker\n(n={len(cur)})"],
                       widths=.55, patch_artist=True, medianprops=dict(color=NAVY, lw=2))
    for b, c in zip(bp["boxes"], [BLUE, "#e6c8c8"]):
        b.set_facecolor(c); b.set_edgecolor(GREY)
    rng = np.random.default_rng(1)
    for k, v in enumerate([nev, cur]):
        ax[0].scatter(np.full(len(v), k + 1) + rng.normal(0, .07, len(v)), v, s=8,
                      color=GREY, alpha=.35, zorder=3)
ax[0].set_ylabel("cg05575921 beta (AHRR)")
ax[0].set_title(f"A  The control R07 could not power\nnow holds: p = {pc3['p']:.4f}",
                fontsize=10, color=NAVY, loc="left")

# ---- B: Q1 pooled vs within-histology (the correction) ----
genes = ["KEAP1", "SMARCA4", "KMT2D"]
A = aud["A_histology_confound"]
x = np.arange(len(genes)); w = 0.27
pooled = [A[g]["pooled_frac_sig"] * 100 for g in genes]
luad = [(A[g]["TCGA-LUAD_frac_sig"] or 0) * 100 for g in genes]
lusc = [(A[g]["TCGA-LUSC_frac_sig"] if A[g]["TCGA-LUSC_frac_sig"] is not None else np.nan) * 100
        for g in genes]
ax[1].bar(x - w, pooled, w, color=GREY, label="pooled (confounded)", edgecolor="none")
ax[1].bar(x, luad, w, color=NAVY, label="within LUAD", edgecolor="none")
ax[1].bar(x + w, lusc, w, color=BLUE, label="within LUSC", edgecolor=GREY, lw=.5)
for k, g in enumerate(genes):
    if A[g]["TCGA-LUSC_frac_sig"] is None:
        ax[1].text(k + w, 0.6, f"n={A[g]['TCGA-LUSC_n_mut']}\ntoo few", fontsize=6.5,
                   ha="center", color=GREY)
# mark the zero bars explicitly so they read as measured zero, not missing data
for k, g in enumerate(genes):
    if luad[k] < 0.05:
        ax[1].text(k, 0.35, "0.0", fontsize=7, ha="center", color=NAVY)
ax[1].annotate("entirely\nhistology", xy=(2 - w, pooled[2]), xytext=(1.02, 24.5),
               fontsize=8, color=RED, ha="center",
               arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
ax[1].set_ylim(0, 35)
ax[1].set_xticks(x); ax[1].set_xticklabels(genes)
ax[1].set_ylabel("% of 20,000 probes differential (q<0.05)")
ax[1].legend(fontsize=7.5, frameon=False, loc="upper center", ncol=1,
             bbox_to_anchor=(0.62, 1.0))
ax[1].set_title("B  Q1 corrected: KMT2D's phenotype\nis a LUAD/LUSC artefact",
                fontsize=10, color=NAVY, loc="left")

# ---- C: Q2/Q3 — methylation adds nothing prognostically ----
q2 = res["Q2_methylation_signature"]; q3 = res["Q3_incremental"]
d = aud["D_q2_power"]
labels = ["methylation\nalone", "mutation+stage\n+age", "  + methylation"]
vals = [q2["cindex"], q3["cindex_base"], q3["cindex_with_methylation"]]
cols_ = [BLUE, NAVY, NAVY]
ax[2].bar(range(3), vals, 0.55, color=cols_, edgecolor=GREY, lw=.5)
ax[2].errorbar([0], [q2["cindex"]],
               yerr=[[q2["cindex"] - q2["ci95"][0]], [q2["ci95"][1] - q2["cindex"]]],
               fmt="none", ecolor=GREY, capsize=4, lw=1.5)
ax[2].axhline(0.5, color=GREY, ls=":", lw=1.2)
ax[2].axhline(d["min_detectable_cindex_80pct"], color=RED, ls="--", lw=1.2)
ax[2].text(2.42, d["min_detectable_cindex_80pct"] + .004,
           f"detectable\nceiling {d['min_detectable_cindex_80pct']:.3f}",
           fontsize=7, color=RED, ha="right")
for k, v in enumerate(vals):
    ax[2].text(k, v + .008, f"{v:.3f}", ha="center", fontsize=9, color=NAVY)
ax[2].set_xticks(range(3)); ax[2].set_xticklabels(labels, fontsize=8.5)
ax[2].set_ylim(0.40, 0.68); ax[2].set_ylabel("C-index (cross-validated)")
ax[2].set_title(f"C  Q3: methylation adds nothing\nLR p = {q3['lr_p']:.2f}",
                fontsize=10, color=NAVY, loc="left")

for a in ax:
    a.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("figures/r08_panels.pdf")
plt.savefig("figures/r08_panels.png", dpi=150)
print("figures/r08_panels.{pdf,png}")
