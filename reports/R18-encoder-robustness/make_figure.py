#!/usr/bin/env python3
"""R18 figure — reads evidence/results.json only."""
import json

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

r = json.load(open("evidence/results.json"))
S, M = r["subtype"], r["methylation"]
ph, dv = S["tcga"], S["dinov2-tcga"]
mph, mdv = M["tcga"], M["dinov2-tcga"]
NAVY, GREY, RED, BLUE, GREEN, AMBER = "#14284b", "#5a5a5a", "#961e1e", "#c8d4e6", "#2f6b4f", "#b07d2b"
TG = ["keap1_sig", "global", "island", "opensea", "tss200", "body"]
PRETTY = ["KEAP1 sig", "global", "island", "open sea", "TSS200", "body"]

fig, ax = plt.subplots(1, 3, figsize=(15.2, 4.9))

# ---- A: subtype, grouped vs random, both encoders ----
a = ax[0]
x = np.arange(2)
w = 0.34
for i, (nm, d, col) in enumerate((("Phikon-v2\n(histology)", ph, NAVY),
                                  ("dinov2-large\n(natural images)", dv, AMBER))):
    a.bar(x + (i - 0.5) * w, [d["subtype_grouped"], d["subtype_random"]], w,
          color=col, edgecolor="k", lw=.5, label=nm)
    a.annotate("", xy=(x[1] + (i - 0.5) * w, d["subtype_random"]),
               xytext=(x[0] + (i - 0.5) * w, d["subtype_grouped"]),
               arrowprops=dict(arrowstyle="->", color=RED, lw=1.6))
    a.text(0.5 + (i - 0.5) * w, (d["subtype_grouped"] + d["subtype_random"]) / 2 + .02,
           f"{d['subtype_inflation']:+.3f}", color=RED, fontsize=9.5, ha="center",
           fontweight="bold")
a.set_xticks(x)
a.set_xticklabels(["site-grouped folds\n(sites held out)", "random folds\n(sites shared)"],
                  fontsize=9)
a.set_ylabel("subtype AUROC")
a.set_ylim(0.5, 1.04)
a.axhline(0.5, color=GREY, ls=":", lw=1)
a.legend(fontsize=8.5, frameon=False, loc="upper left")
a.set_title("A  Subtype leakage reproduces — and is LARGER\n"
            "an encoder that has never seen histology shows the same gap",
            fontsize=10, color=NAVY, loc="left")

# ---- B: methylation per target, inflation ----
b = ax[1]
y = np.arange(6)
b.barh(y + .18, [mph["per_target"][t]["inflation"] for t in TG], .34,
       color=NAVY, edgecolor="k", lw=.5, label="Phikon-v2")
b.barh(y - .18, [mdv["per_target"][t]["inflation"] for t in TG], .34,
       color=AMBER, edgecolor="k", lw=.5, label="dinov2-large")
b.axvline(0, color="k", lw=1)
b.set_yticks(y)
b.set_yticklabels(PRETTY, fontsize=9)
b.set_xlabel(r"inflation in Spearman $\rho$  (random $-$ grouped)")
b.legend(fontsize=8.5, frameon=False, loc="lower right")
b.set_title("B  Methylation leakage does NOT reproduce\n"
            f"mean {mph['mean_inflation']:+.4f} vs {mdv['mean_inflation']:+.4f}; "
            f"{mph['n_targets_inflated']}/6 vs {mdv['n_targets_inflated']}/6 inflated",
            fontsize=10, color=NAVY, loc="left")

# ---- C: the dissociation, on relative inflation ----
c = ax[2]
labels = ["subtype", "methylation\n(mean)"]
phv = [ph["subtype_relative_inflation"],
       mph["mean_inflation"] / (1 - mph["mean_grouped_rho"])]
dvv = [dv["subtype_relative_inflation"],
       mdv["mean_inflation"] / (1 - mdv["mean_grouped_rho"])]
xx = np.arange(2)
c.bar(xx - .18, phv, .34, color=NAVY, edgecolor="k", lw=.5, label="Phikon-v2")
c.bar(xx + .18, dvv, .34, color=AMBER, edgecolor="k", lw=.5, label="dinov2-large")
for i, (p, d) in enumerate(zip(phv, dvv)):
    c.text(i - .18, p + .015, f"{p:.3f}", ha="center", fontsize=8.7, color=NAVY)
    c.text(i + .18, d + .015, f"{d:.3f}", ha="center", fontsize=8.7, color=AMBER)
c.set_xticks(xx)
c.set_xticklabels(labels, fontsize=9.5)
c.set_ylabel("relative inflation  (random $-$ grouped) / (1 $-$ grouped)")
c.axhline(0, color="k", lw=1)
c.legend(fontsize=8.5, frameon=False)
c.set_title("C  Relative inflation, so a weaker encoder is judged fairly\n"
            "the dissociation is not a headroom artefact",
            fontsize=10, color=NAVY, loc="left")

for x_ in ax:
    x_.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("figures/r18_panels.pdf")
plt.savefig("figures/r18_panels.png", dpi=150)
print("figures/r18_panels.{pdf,png}")
