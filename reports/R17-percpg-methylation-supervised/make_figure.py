#!/usr/bin/env python3
"""R17 figure — reads evidence/ only."""
import json

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "evidence/"
r = json.load(open(E + "results.json"))
NAVY, GREY, RED, BLUE, GREEN, AMBER = "#14284b", "#5a5a5a", "#961e1e", "#c8d4e6", "#2f6b4f", "#b07d2b"

rho = {k: np.load(E + f"per_cpg_rho_{v}.npy") for k, v in (
    ("meth", "meth_supervised"), ("sub", "subtype_supervised"),
    ("sublabel", "subtype_label"), ("globalmean", "globalmean"),
    ("null", "null_single_perm"))}
thr = r["null"]["single_permutation_p95"]

fig, ax = plt.subplots(2, 2, figsize=(12.4, 8.6))

# ---- A: the decisive matched-dimensionality comparison ----
a = ax[0, 0]
bins = np.linspace(-0.25, 0.45, 130)
a.hist(rho["null"], bins=bins, color=BLUE, label="permutation null", lw=.4, edgecolor=GREY)
a.hist(rho["sub"], bins=bins, histtype="step", color=GREY, lw=1.9,
       label=f"subtype-supervised 512-d ({r['arms']['sub']['n_above_null']:,})")
a.hist(rho["meth"], bins=bins, histtype="step", color=RED, lw=1.9,
       label=f"methylation-supervised 512-d ({r['arms']['meth']['n_above_null']:,})")
a.axvline(thr, color=NAVY, ls="--", lw=1.2)
a.text(thr + .008, a.get_ylim()[1] * .62, "null p95", fontsize=8, color=NAVY)
a.set_xlabel(r"per-CpG Spearman $\rho$ (out-of-fold, site-grouped)")
a.set_ylabel("CpGs")
a.legend(fontsize=7.6, frameon=False, loc="upper right")
a.set_title("A  Supervising on methylation does not help genome-wide\n"
            "matched dimensionality, folds, ridge — only the target differs",
            fontsize=10, color=NAVY, loc="left")

# ---- B: all four arms above null ----
b = ax[1, 0] if False else ax[0, 1]
order = ["meth", "sub", "sublabel", "globalmean"]
lab = ["methylation-\nsupervised\n512-d", "subtype-\nsupervised\n512-d",
       "subtype label\n1-d", "observed global\nmean methylation\n1-d"]
vals = [r["arms"][k]["n_above_null"] for k in order]
cols = [RED, GREY, GREEN, AMBER]
bar = b.bar(range(4), vals, color=cols, edgecolor=NAVY, lw=.6)
for i, v in enumerate(vals):
    b.text(i, v + 6000, f"{v:,}\n{v/399579*100:.1f}%", ha="center", fontsize=8.4, color=NAVY)
b.set_xticks(range(4))
b.set_xticklabels(lab, fontsize=7.8)
b.set_ylabel("CpGs above the permutation null (of 399,579)")
b.set_ylim(0, 355000)
b.axhline(vals[0], color=RED, ls=":", lw=1)
b.set_title("B  Both embeddings lose to a single scalar\n"
            "the 512-d vs 1-d gap is a variance effect, not a biology one",
            fontsize=10, color=NAVY, loc="left")

# ---- C: exact reproduction of m13 ----
c = ax[1, 0]
m13 = r["m13_reference"]
x = np.arange(2)
w = .36
c.bar(x - w/2, [m13["n_above_null_subtype_emb"], m13["n_above_null_subtype_label"]],
      w, color=BLUE, edgecolor=NAVY, lw=.6, label="R15 / m13 as published")
c.bar(x + w/2, [m13["reproduced_here"]["subtype_emb"], m13["reproduced_here"]["subtype_label"]],
      w, color=GREEN, edgecolor=NAVY, lw=.6, label="recomputed here")
for i, (p, q) in enumerate(zip(
        [m13["n_above_null_subtype_emb"], m13["n_above_null_subtype_label"]],
        [m13["reproduced_here"]["subtype_emb"], m13["reproduced_here"]["subtype_label"]])):
    c.text(i, max(p, q) + 7000, f"{p:,}\nvs {q:,}", ha="center", fontsize=8.4,
           color=GREEN if p == q else RED)
c.set_xticks(x)
c.set_xticklabels(["subtype-supervised\nembedding", "subtype label\nalone"], fontsize=8.4)
c.set_ylabel("CpGs above null")
c.set_ylim(0, 265000)
c.legend(fontsize=8, frameon=False, loc="upper left")
c.set_title("C  The reimplementation reproduces R15 to the unit\n"
            "so the new arm's null result is not a coding difference",
            fontsize=10, color=NAVY, loc="left")

# ---- D: genomic context ----
d = ax[1, 1]
ctx = r["genomic_context"]
keys = ["Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"]
xs = np.arange(len(keys))
d.plot(xs, [ctx["A_meth"][k]["enrichment"] for k in keys], "o-", color=RED, lw=1.8,
       label=f"methylation-supervised (spread {r['context_spread']['A_meth']})")
d.plot(xs, [ctx["A_sub"][k]["enrichment"] for k in keys], "s-", color=GREY, lw=1.8,
       label=f"subtype-supervised (spread {r['context_spread']['A_sub']})")
# R10 reported only Island and OpenSea; plot those two points and nothing between them,
# because a connecting line would imply shore/shelf values we do not have.
d.plot([0, len(keys) - 1], [r["r10_reference_context"]["Island"],
                            r["r10_reference_context"]["OpenSea"]], "^", color=GREEN, ms=11,
       label="R10, from the assay itself (2 classes reported)")
for xx, yy in ((0, r["r10_reference_context"]["Island"]),
               (len(keys) - 1, r["r10_reference_context"]["OpenSea"])):
    d.annotate(f"{yy:.2f}", (xx, yy), textcoords="offset points",
               xytext=(12 if xx == 0 else -26, -3), fontsize=8, color=GREEN)
d.axhline(1.0, color=NAVY, ls=":", lw=1.1)
d.set_xticks(xs)
d.set_xticklabels(keys, fontsize=8, rotation=20)
d.set_ylabel("enrichment among predictable CpGs")
d.legend(fontsize=7.6, frameon=False, loc="upper right")
d.set_title("D  Neither embedding recovers the assay's structure\n"
            "and A_meth's mild slope points the opposite way to R10",
            fontsize=10, color=NAVY, loc="left")

for row in ax:
    for x_ in row:
        x_.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("figures/r17_panels.pdf")
plt.savefig("figures/r17_panels.png", dpi=150)
print("figures/r17_panels.{pdf,png}")
