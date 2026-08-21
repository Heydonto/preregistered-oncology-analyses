#!/usr/bin/env python3
"""R24 figure — reads evidence/ only."""
import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

E = "evidence/"
r = json.load(open(E + "results.json"))
dt = json.load(open(E + "dose_trend.json"))
kd = np.load(E + "kd_lfc_PC9OR.npy"); res = np.load(E + "res_lfc_PC9.npy")
NAVY,GREY,RED,AMBER,GREEN,BLUE = "#14284b","#5a5a5a","#961e1e","#b07d2b","#2f6b4f","#c8d4e6"

fig, ax = plt.subplots(1, 4, figsize=(18.4, 4.7))

# A: knockdown verified
a = ax[0]
pc1 = r["PC1_smarca4_knockdown"]
v = [pc1["PC9OR_log2FC"], pc1["YU005C_log2FC"]]
a.bar([0,1], v, .5, color=GREEN, edgecolor="k", lw=.6)
a.axhline(-0.5, color=RED, ls=":", lw=1.3)
a.text(1.32, -0.5, "pre-declared\nbar $-0.5$", fontsize=7.6, color=RED, va="center")
for i,x in enumerate(v): a.text(i, x-.17, f"{x:.2f}", ha="center", fontsize=10, color=NAVY)
a.axhline(0, color="k", lw=1); a.set_xticks([0,1]); a.set_xticklabels(["PC9-OR","YU005C"], fontsize=9)
a.set_ylabel(r"SMARCA4 log$_2$FC (shSMARCA4 vs scramble)"); a.set_ylim(-2.9,.35)
a.set_title("A  The epigenetic modulation worked\npositive control, declared before looking",
            fontsize=9.6, color=NAVY, loc="left")

# B: resistance programme is line-specific
b = ax[1]
h1 = r["H1_cross_line_consistency"]; ks = list(h1)
b.bar(range(len(ks)), [h1[k] for k in ks], .55, color=AMBER, edgecolor="k", lw=.5)
for i,k in enumerate(ks): b.text(i, h1[k]+.012, f"{h1[k]:.2f}", ha="center", fontsize=9.5, color=NAVY)
b.axhline(1.0, color=GREY, ls="--", lw=1.1)
b.text(2.3, .95, "a shared\nprogramme", fontsize=7.6, color=GREY, va="top", ha="right")
b.set_xticks(range(len(ks))); b.set_xticklabels([k.replace("_vs_","\nvs ") for k in ks], fontsize=8)
b.set_ylabel(r"Spearman $\rho$ of per-gene log$_2$FC"); b.set_ylim(0,1.08)
b.set_title("B  Resistance is largely line-specific\nnot one programme across three lines",
            fontsize=9.6, color=NAVY, loc="left")

# C: the primary test, both models
c = ax[2]
xs, lbl = [], []
for j,(mdl,key) in enumerate((("PC9-OR","H2_PC9OR"),("YU005C","H3_YU005C_replication"))):
    d = r[key]
    for i,(nm,s) in enumerate((("UP",d["resistance_up"]),("DOWN",d["resistance_down"]))):
        x = j*2.4 + i
        mid = (s["null_lo"]+s["null_hi"])/2
        c.errorbar(x, mid, yerr=[[mid-s["null_lo"]],[s["null_hi"]-mid]], fmt="_",
                   color=GREY, ms=30, lw=2.2, capsize=5)
        c.scatter(x, s["mean_kd_lfc"], s=95, color=RED, zorder=4)
        xs.append(x); lbl.append(f"{mdl}\n{nm}")
c.axhline(0, color="k", lw=.9)
c.set_xticks(xs); c.set_xticklabels(lbl, fontsize=7.8)
c.set_ylabel(r"mean knockdown log$_2$FC")
c.plot([],[], "_", color=GREY, ms=13, lw=2, label="null 95%")
c.scatter([],[], color=RED, s=55, label="observed")
c.legend(fontsize=7.6, frameon=False, loc="lower left")
c.set_title("C  PARTIAL in PC9-OR, and it flips in YU005C\nreversal needed UP below and DOWN above",
            fontsize=9.6, color=NAVY, loc="left")

# D: the monotonic dose trend
d4 = ax[3]
bins = dt["bins"]
m = [x["mean_kd_lfc"] for x in bins]; n = [x["n"] for x in bins]
d4.bar(range(len(bins)), m, .55, color=[BLUE if q>0 else RED for q in m], edgecolor="k", lw=.5)
for i,(q,nn) in enumerate(zip(m,n)):
    d4.text(i, q+.008 if q>0 else q-.024, f"{q:+.3f}\nn={nn:,}", ha="center", fontsize=7.6,
            color=NAVY, va="bottom" if q>0 else "top")
d4.axhline(dt["all_gene_mean"], color=GREY, ls="--", lw=1.2)
d4.text(4.4, dt["all_gene_mean"], "all-gene\nmean", fontsize=7.4, color=GREY, va="center")
d4.axhline(0, color="k", lw=.9)
d4.set_xticks(range(len(bins)))
d4.set_xticklabels([x["bin"].replace("3-99",">3") for x in bins], fontsize=8)
d4.set_xlabel(r"|resistance log$_2$FC| bin"); d4.set_ylabel(r"mean knockdown log$_2$FC")
d4.set_ylim(-.30,.10)
d4.set_title("D  Monotonic, and in BOTH directions\nsuppression of the programme, not reversal",
             fontsize=9.6, color=NAVY, loc="left")

for x_ in ax: x_.spines[["top","right"]].set_visible(False)
plt.tight_layout()
plt.savefig("figures/r24_panels.pdf"); plt.savefig("figures/r24_panels.png", dpi=150)
print("figures/r24_panels.{pdf,png}")
