#!/usr/bin/env python3
"""Regenerate the R09 figure entirely from evidence/consolidated.json."""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

c = json.load(open("evidence/consolidated.json"))
NAVY, GREY, RED, BLUE = "#14284b", "#5a5a5a", "#961e1e", "#c8d4e6"
R06 = {"KEAP1": 2.07, "SMARCA4": 1.73, "STK11": 1.63, "TP53": 1.25}
R06CI = {"KEAP1": [1.63, 2.62], "SMARCA4": [1.27, 2.35]}

fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.2))

# ---- A: positive controls across all three settings ----
tp = {x["gene"]: x for x in c["tcga"]["primary"]}
po = {x["gene"]: x for x in c["matched"]["pooled"]}
groups = ["TP53", "STK11"]
series = [("R06 cohort\n(618 ev)", [R06["TP53"], R06["STK11"]], NAVY),
          ("TCGA\n(385 ev)", [tp["TP53"]["hr"], tp["STK11"]["hr"]], GREY),
          ("matched pool\n(231 ev)", [po["TP53"]["hr"], po["STK11"]["hr"]], BLUE)]
x = np.arange(2); w = 0.26
for k, (lab, vals, col) in enumerate(series):
    ax[0].bar(x + (k - 1) * w, vals, w, color=col, label=lab,
              edgecolor=GREY if col == BLUE else "none", lw=.5)
ax[0].axhline(1.0, color=RED, ls="--", lw=1.2)
ax[0].text(1.46, 1.02, "no effect", fontsize=7, color=RED, ha="right")
ax[0].set_xticks(x); ax[0].set_xticklabels(groups)
ax[0].set_ylabel("hazard ratio")
ax[0].set_ylim(0, 2.0)
ax[0].legend(fontsize=7, frameon=False, loc="upper left")
ax[0].annotate("wrong\ndirection", xy=(1 + w, po["STK11"]["hr"]), xytext=(1.30, 0.45),
               fontsize=7.5, color=RED, ha="center",
               arrowprops=dict(arrowstyle="->", color=RED, lw=1.1))
ax[0].set_title("A  The controls fail everywhere\nbut in R06's own cohort",
                fontsize=10, color=NAVY, loc="left")

# ---- B: the chromatin genes, R06 vs both attempts ----
genes = ["KEAP1", "SMARCA4", "KMT2D"]
y = np.arange(len(genes))[::-1]
for k, g in enumerate(genes):
    yy = y[k]
    if g in R06CI:
        ax[1].plot(R06CI[g], [yy + .22, yy + .22], color=NAVY, lw=1.4, zorder=2)
        ax[1].plot([R06[g]], [yy + .22], "o", color=NAVY, ms=7, zorder=3,
                   label="R06 (measured)" if k == 0 else None)
    else:
        ax[1].plot([], [])
    t = tp[g]
    ax[1].plot([t["hr"]], [yy], "s", color=GREY, ms=6,
               label="TCGA" if k == 0 else None)
    p = po[g]
    ax[1].plot(p["ci"], [yy - .22, yy - .22], color=BLUE, lw=1.4)
    ax[1].plot([p["hr"]], [yy - .22], "D", color=BLUE, ms=6, markeredgecolor=GREY,
               markeredgewidth=.5, label="matched pool" if k == 0 else None)
ax[1].axvline(1.0, color=RED, ls="--", lw=1.2)
ax[1].set_yticks(y); ax[1].set_yticklabels(genes)
ax[1].set_xscale("log")
ax[1].minorticks_off()          # log minor ticks collide with the explicit labels
ax[1].set_xlim(0.28, 3.4)
ax[1].set_xticks([0.3, 0.5, 1, 2, 3]); ax[1].set_xticklabels(["0.3", "0.5", "1", "2", "3"])
ax[1].set_ylim(-0.6, 2.6)
ax[1].set_xlabel("hazard ratio (log scale)")
ax[1].legend(fontsize=7, frameon=False, loc="lower right")
ax[1].set_title("B  No chromatin gene reproduces\nR06's effect size",
                fontsize=10, color=NAVY, loc="left")

# ---- C: the screen ----
s = c["screen"]
rows = [r for r in s["rows"] if r.get("verdict")]
buckets = {"no OS endpoint": 0, "too few events": 0, "wrong setting": 0, "eligible": 0}
for r in rows:
    v = r["verdict"]
    if v == "no OS endpoint":
        buckets["no OS endpoint"] += 1
    elif v == "ELIGIBLE":
        buckets["eligible"] += 1
    elif v.startswith("only"):
        buckets["too few events"] += 1
    elif "not advanced" in v:
        buckets["wrong setting"] += 1
assert sum(buckets.values()) == len(rows), (buckets, len(rows))
labs = list(buckets); vals = [buckets[k] for k in labs]
pos = list(range(len(labs)))[::-1]
ax[2].barh(pos, vals, 0.6, color=[GREY] * 3 + [NAVY], edgecolor="none")
for k, v in enumerate(vals):
    ax[2].text(v + 0.4, pos[k], str(v), va="center", fontsize=9, color=NAVY)
ax[2].set_yticks(pos); ax[2].set_yticklabels(labs, fontsize=8.5)
ax[2].set_xlabel(f"lung studies (of {len(rows)} screened)")
ax[2].set_xlim(0, max(vals) * 1.25)
ax[2].set_title("C  Systematic screen: 2 of 41 cohorts\neven qualified to be tested",
                fontsize=10, color=NAVY, loc="left")

for a in ax:
    a.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("figures/r09_panels.pdf")
plt.savefig("figures/r09_panels.png", dpi=150)
print("figures/r09_panels.{pdf,png}")
