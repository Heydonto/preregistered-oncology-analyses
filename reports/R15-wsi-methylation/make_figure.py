#!/usr/bin/env python3
"""R15 figure. Panels A-C read evidence/ only; panel D reads the methylation-MIL output."""
import json, os
import numpy as np
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "evidence/"
r = json.load(open(E + "results.json"))
NAVY, GREY, RED, BLUE, GREEN = "#14284b", "#5a5a5a", "#961e1e", "#c8d4e6", "#2f6b4f"

fig, ax = plt.subplots(2, 2, figsize=(12.4, 8.4))

# ---- A: what fold assignment alone buys ----
a = ax[0, 0]
L = r["HEADLINE_site_leakage"]
x = np.arange(2); w = 0.35
gr = [L["subtype"]["grouped"], L["keap1"]["grouped"]]
rn = [L["subtype"]["random"], L["keap1"]["random"]]
a.bar(x - w/2, gr, w, color=NAVY, label="site-grouped (honest)", edgecolor=GREY, lw=.5)
a.bar(x + w/2, rn, w, color=BLUE, label="random folds (sites shared)", edgecolor=GREY, lw=.5)
for i, (g, n) in enumerate(zip(gr, rn)):
    a.text(i - w/2 - 0.02, g + .012, f"{g:.3f}", ha="right", fontsize=9, color=NAVY)
    a.text(i + w/2, n + .012, f"{n:.3f}", ha="center", fontsize=9, color=NAVY)
    a.annotate("", xy=(i + w/2, n), xytext=(i - w/2, g),
               arrowprops=dict(arrowstyle="->", color=RED, lw=1.3))
    a.text(i, max(g, n) + .055, f"+{n-g:.3f}", ha="center", fontsize=9.5,
           color=RED, fontweight="bold")
a.axhline(0.5, color=GREY, ls=":", lw=1.2)
a.set_xticks(x); a.set_xticklabels(["LUAD vs LUSC subtype", "KEAP1 mutation"])
a.set_ylabel("AUROC"); a.set_ylim(0.45, 1.16)
a.legend(fontsize=7.5, frameon=False, loc="upper left", bbox_to_anchor=(0.02, 0.99))
a.set_title("A  Identical data and model; only the folds differ\n"
            "subtype is site-confounded, KEAP1 is not", fontsize=10, color=NAVY, loc="left")

# ---- B: KEAP1 from slide vs from array ----
b = ax[0, 1]
q1 = r["Q1_keap1_from_wsi"]
vals = [q1["grouped_auroc"], q1["r10_from_methylation"]]
b.bar([0, 1], vals, 0.5, color=[NAVY, GREEN], edgecolor=GREY, lw=.5)
lo, hi = q1["grouped_ci95"]
b.errorbar([0], [vals[0]], yerr=[[vals[0]-lo], [hi-vals[0]]], fmt="none",
           ecolor=GREY, capsize=5, lw=1.5)
for i, v in enumerate(vals):
    b.text(i, v + .015, f"{v:.3f}", ha="center", fontsize=10, color=NAVY)
b.axhline(0.5, color=GREY, ls=":", lw=1.2)
b.set_xticks([0, 1])
b.set_xticklabels(["from H&E slide\n(this work)", "from 450k array\n(R10)"], fontsize=8.5)
b.set_ylabel("AUROC"); b.set_ylim(0.45, 1.0)
b.text(0.5, 0.52, f"shortfall {q1['shortfall_vs_methylation']:.3f}", ha="center",
       fontsize=8.5, color=RED)
b.set_title("B  Morphology does not replace the assay\nfor KEAP1 status",
            fontsize=10, color=NAVY, loc="left")

# ---- C: the first, negative answer — subtype-supervised embedding ----
c = ax[1, 0]
q2 = r["Q2_genomewide"]
n_tot = q2["probes"]
vals = [q2["n_predictable"] / n_tot * 100,
        q2["subtype_only_baseline"]["n_predictable"] / n_tot * 100]
c.bar([0, 1], vals, 0.5, color=[BLUE, GREY], edgecolor=GREY, lw=.5)
for i, v in enumerate(vals):
    c.text(i, v + 1.2, f"{v:.1f}%", ha="center", fontsize=10, color=NAVY)
c.set_xticks([0, 1])
c.set_xticklabels(["subtype-supervised\nWSI embedding", "subtype label alone\n(one binary variable)"],
                  fontsize=8.5)
c.set_ylabel(f"% of {n_tot:,} CpGs above permutation null")
c.set_ylim(0, 70)
c.text(0.5, 62, "a 512-d embedding beaten by\none binary variable", ha="center",
       fontsize=8.5, color=RED)
c.set_title("C  First attempt: supervise on subtype, ask about methylation\n"
            "the negative was an artefact of that choice",
            fontsize=10, color=NAVY, loc="left")

# ---- D: methylation-supervised, partial correlations beyond subtype ----
d = ax[1, 1]
if os.path.exists(E + "mil_meth_output.npz"):
    z = np.load(E + "mil_meth_output.npz", allow_pickle=True)
    names = [str(x) for x in z["target_names"]]
    P, O = z["pred"], z["observed"]
    pids = [str(x) for x in z["pids"]]
    lab = {p["pid"]: p for p in json.load(open(E + "mil_input.json"))["patients"]}
    sub = np.array([lab[p]["subtype"] for p in pids], float)

    def res(y, x):
        X = np.c_[np.ones(len(x)), x]
        return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]

    pr = [spearmanr(res(O[:, i], sub), res(P[:, i], sub)).statistic
          for i in range(len(names))]
    order = np.argsort(pr)[::-1]
    pretty = {"keap1_sig": "KEAP1 signature", "global": "global", "island": "CpG island",
              "opensea": "open sea", "tss200": "TSS200 promoter", "body": "gene body"}
    d.barh(range(len(order))[::-1], [pr[i] for i in order], 0.6,
           color=NAVY, edgecolor=GREY, lw=.5)
    for k, i in enumerate(order):
        d.text(pr[i] + .006, list(range(len(order)))[::-1][k], f"{pr[i]:.3f}",
               va="center", fontsize=9, color=NAVY)
    d.set_yticks(range(len(order))[::-1])
    d.set_yticklabels([pretty.get(names[i], names[i]) for i in order], fontsize=8.5)
    d.axvline(0, color=GREY, lw=1)
    d.set_xlim(0, 0.30)
    d.set_xlabel("partial Spearman $\\rho$, controlling for subtype")
    d.set_title("D  Supervise on methylation instead: 6/6 targets\ncarry signal beyond subtype",
                fontsize=10, color=NAVY, loc="left")
for row in ax:
    for x_ in row:
        x_.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("figures/r15_panels.pdf")
plt.savefig("figures/r15_panels.png", dpi=150)
print("figures/r15_panels.{pdf,png}")
