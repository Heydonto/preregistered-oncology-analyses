#!/usr/bin/env python3
"""R21 figure — reads evidence/results.json only."""
import json

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

r = json.load(open("evidence/results.json"))
E = r["encoders"]
ORDER = ["tcga", "dinov2-large", "uni-tcga", "hopt-tcga", "virchow2-tcga"]
SHORT = {"tcga": "Phikon-v2", "dinov2-large": "dinov2-large", "uni-tcga": "UNI",
         "hopt-tcga": "H-optimus-0", "virchow2-tcga": "Virchow2"}
NAVY, GREY, RED, AMBER, GREEN = "#14284b", "#5a5a5a", "#961e1e", "#b07d2b", "#2f6b4f"
col = {k: (AMBER if E[k]["corpus"] == "natural images" else NAVY) for k in ORDER}
TG = ["keap1_sig", "global", "island", "opensea", "tss200", "body"]

fig, ax = plt.subplots(1, 3, figsize=(16.2, 5.1))

# ---- A: subtype, grouped vs random, all five ----
a = ax[0]
x = np.arange(len(ORDER))
a.bar(x - .19, [E[k]["subtype_grouped"] for k in ORDER], .38, color=[col[k] for k in ORDER],
      edgecolor="k", lw=.5, label="site-grouped (sites held out)")
a.bar(x + .19, [E[k]["subtype_random"] for k in ORDER], .38, color=[col[k] for k in ORDER],
      alpha=.45, edgecolor="k", lw=.5, label="random folds (sites shared)")
for i, k in enumerate(ORDER):
    a.text(i, E[k]["subtype_random"] + .012,
           f"{E[k]['subtype_inflation']:+.3f}", ha="center", fontsize=8.6, color=RED,
           fontweight="bold")
a.set_xticks(x)
a.set_xticklabels([SHORT[k] for k in ORDER], fontsize=8.4, rotation=18)
a.set_ylabel("subtype AUROC")
a.set_ylim(0.5, 1.06)
a.legend(fontsize=8, frameon=False, loc="lower left")
a.set_title("A  Subtype leakage in all five encoders\n"
            "but the gap shrinks sharply for the newer models",
            fontsize=10, color=NAVY, loc="left")

# ---- B: methylation inflation, split by training corpus ----
b = ax[1]
vals = [E[k]["meth_mean_inflation"] for k in ORDER]
b.bar(x, vals, .55, color=[col[k] for k in ORDER], edgecolor="k", lw=.5)
for i, k in enumerate(ORDER):
    b.text(i, vals[i] + .0022, f"{vals[i]:+.4f}\n{E[k]['meth_targets_inflated']}/6",
           ha="center", fontsize=8.3, color=NAVY)
b.axhline(0, color="k", lw=1)
b.set_xticks(x)
b.set_xticklabels([SHORT[k] for k in ORDER], fontsize=8.4, rotation=18)
b.set_ylabel(r"mean methylation inflation ($\Delta\rho$)")
b.set_ylim(0, max(vals) * 1.32)
b.text(.98, .95, "navy = histology-pretrained\namber = natural images",
       transform=b.transAxes, ha="right", va="top", fontsize=8.4, color=GREY)
b.set_title("B  Methylation leakage tracks the TRAINING CORPUS\n"
            "4/4 histology encoders inflate all 6 targets; the natural-image one does not",
            fontsize=9.4, color=NAVY, loc="left")

# ---- C: the post-hoc observation, labelled as such ----
c = ax[2]
ph = r["post_hoc_observation"]
gx = [E[k]["subtype_grouped"] for k in ORDER]
gy = [E[k]["subtype_relative_inflation"] for k in ORDER]
for i, k in enumerate(ORDER):
    c.scatter(gx[i], gy[i], s=120, color=col[k], edgecolor="k", zorder=3)
    c.annotate(SHORT[k], (gx[i], gy[i]), textcoords="offset points", xytext=(8, -3),
               fontsize=8.2, color=NAVY)
c.set_xlabel("site-disjoint subtype AUROC (capability)")
c.set_ylabel("relative inflation  (random $-$ grouped)/(1 $-$ grouped)")
c.set_xlim(0.70, 1.05)
c.text(.03, .06, f"Spearman $\\rho$ = {ph['spearman_rho']:.2f}, p = {ph['p']:.3f}, n = {ph['n_encoders']}\n"
       "POST-HOC OBSERVATION — NOT A TESTED FINDING\n"
       "noticed after reading all five; at n=5 no rank\ncorrelation can reach p<0.05 by this test",
       transform=c.transAxes, fontsize=8.2, color=RED, va="bottom")
c.set_title("C  A pattern we are deliberately NOT claiming\n"
            "better encoders look less site-penalised",
            fontsize=9.4, color=NAVY, loc="left")

for x_ in ax:
    x_.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("figures/r21_panels.pdf")
plt.savefig("figures/r21_panels.png", dpi=150)
print("figures/r21_panels.{pdf,png}")
