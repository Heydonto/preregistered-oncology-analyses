#!/usr/bin/env python3
"""R23 figure — reads evidence/ only."""
import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

r = json.load(open("evidence/results.json"))
NAVY,GREY,RED,AMBER,GREEN = "#14284b","#5a5a5a","#961e1e","#b07d2b","#2f6b4f"
fig, ax = plt.subplots(1, 3, figsize=(15.6, 5.2))

# A: per-drug predictability
a = ax[0]
soc = r["H1_standard_of_care"]; chrom = r["H2_chromatin"][:12]
names = [d["drug"] for d in soc] + [d["drug"] for d in chrom]
rhos  = [d["rho"] for d in soc] + [d["rho"] for d in chrom]
sig   = [d["holm_p"] < 0.05 for d in soc] + [d["holm_p"] < 0.05 for d in chrom]
cols  = [NAVY]*len(soc) + [AMBER]*len(chrom)
y = np.arange(len(names))[::-1]
a.barh(y, rhos, color=cols, edgecolor="k", lw=.4)
for i,(v,s) in enumerate(zip(rhos,sig)):
    a.text(v+.012, y[i], ("*" if s else "n.s."), va="center", fontsize=7.5,
           color=NAVY if s else GREY)
a.set_yticks(y); a.set_yticklabels(names, fontsize=7.4)
a.set_xlabel(r"out-of-fold Spearman $\rho$ (methylation $\rightarrow$ ln IC$_{50}$)")
a.set_xlim(0, .82)
a.plot([], [], color=NAVY, lw=6, label="standard of care")
a.plot([], [], color=AMBER, lw=6, label="chromatin-targeting")
a.legend(fontsize=8, frameon=False, loc="lower right")
a.set_title("A  Methylation predicts sensitivity for both drug classes\n"
            "* = Holm-significant, 2,000 permutations", fontsize=10, color=NAVY, loc="left")

# B: the primary test and what the control does to it
b = ax[1]
raw = r["VERDICT"]["rho"]; pc = r["H3_POSTHOC_general_sensitivity_control"]
b.bar([0,1], [raw, pc["partial_rho"]], .5,
      color=[RED, GREEN], edgecolor="k", lw=.6)
b.axhline(0, color="k", lw=1)
for t,v in ((-0.20,"thesis supported below"),(0.20,"thesis contradicted above")):
    b.axhline(t, color=GREY, ls=":", lw=1.1)
    b.text(1.42, t, f"{t:+.2f}", fontsize=7.6, color=GREY, va="center")
b.text(0, raw+.03, f"{raw:+.3f}\np={r['VERDICT']['perm_p']:.4f}", ha="center", fontsize=9, color=RED)
b.text(0.72, pc["partial_rho"]-.14, f"{pc['partial_rho']:+.3f}\np={pc['perm_p']:.3f}",
       ha="center", fontsize=9, color=GREEN)
b.set_xticks([0,1])
b.set_xticklabels(["PRE-REGISTERED\n(raw indices)","POST-HOC control\n(general sensitivity\npartialled out)"],
                  fontsize=8.4)
b.set_ylabel(r"Spearman $\rho$: resistance index vs chromatin index")
b.set_ylim(-.42, .78)
b.set_title("B  The pre-registered verdict does not survive its own control\n"
            "the sign reverses once general drug sensitivity is removed",
            fontsize=10, color=NAVY, loc="left")

# C: the confound itself
c = ax[2]
a_i = np.load("evidence/resistance_index.npy"); b_i = np.load("evidence/chromatin_index.npy")
c.scatter(a_i, b_i, s=26, color=NAVY, alpha=.62, edgecolor="none")
z = np.polyfit(a_i, b_i, 1); xs = np.linspace(a_i.min(), a_i.max(), 50)
c.plot(xs, np.polyval(z, xs), color=RED, lw=1.8)
c.set_xlabel("resistance index (mean z ln IC$_{50}$, standard of care)")
c.set_ylabel("chromatin index (mean z ln IC$_{50}$)")
c.text(.03,.95, f"n = {r['VERDICT']['n_lines']} lung cell lines\n"
       r"raw $\rho$ = " + f"{raw:+.3f}" + "\nalmost entirely a general\ndrug-sensitivity axis",
       transform=c.transAxes, fontsize=8.6, va="top", color=NAVY)
c.set_title("C  What the raw association actually looks like\n"
            "lines resistant to one class are resistant to most things",
            fontsize=10, color=NAVY, loc="left")
for x_ in ax: x_.spines[["top","right"]].set_visible(False)
plt.tight_layout(); plt.savefig("figures/r23_panels.pdf"); plt.savefig("figures/r23_panels.png", dpi=150)
print("figures/r23_panels.{pdf,png}")
