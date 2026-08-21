#!/usr/bin/env python3
"""Regenerate the R07 figure from the archived evidence in ./evidence/.
Run from results/R07-methylation-tki-response/. Panels A and B read the source data;
panel C reads only the audit addendum."""
import json, gzip, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "evidence/"
res = json.load(open(E + "results.json"))
add = json.load(open(E + "audit_addendum.json"))
nul = json.load(open(E + "permutation_null.json"))
NAVY, GREY, RED = "#14284b", "#5a5a5a", "#961e1e"

D = "/Users/rezanehzati/quantara-staging/staged/labels/lung/GSE147377/"
hdr = gzip.open(D + "suppl/GSE147377_AverageBeta_Matrix.csv.gz", "rt").readline().rstrip().split(",")
bidx = [i for i, c in enumerate(hdr) if c.endswith(".AVG_Beta")]
ahrr = None
for line in gzip.open(D + "suppl/GSE147377_AverageBeta_Matrix.csv.gz", "rt"):
    q = line.split(",")
    if q[0] == "cg05575921":
        ahrr = np.array([float(q[i]) for i in bidx]); break
smoke = None
for line in gzip.open(D + "matrix/GSE147377_series_matrix.txt.gz", "rt", errors="ignore"):
    if line.startswith("!Sample_characteristics_ch1") and "smoking" in line:
        smoke = np.array([x.split(": ", 1)[-1].strip().strip('"')
                          for x in line.split("\t")[1:]], object)
    if line.startswith("!series_matrix_table_begin"):
        break

fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))

sm, nv = ahrr[smoke == "smoker"], ahrr[smoke == "Never-smoker"]
bp = ax[0].boxplot([nv, sm], tick_labels=[f"Never-smoker\n(n={len(nv)})", f"Smoker\n(n={len(sm)})"],
                   widths=.55, patch_artist=True, medianprops=dict(color=NAVY, lw=2))
for b, c in zip(bp["boxes"], ["#c8d4e6", "#e6c8c8"]):
    b.set_facecolor(c); b.set_edgecolor(GREY)
rng = np.random.default_rng(1)
for i, v in enumerate([nv, sm]):
    ax[0].scatter(np.full(len(v), i + 1) + rng.normal(0, .06, len(v)), v, s=13,
                  color=GREY, alpha=.6, zorder=3)
ax[0].set_ylabel("cg05575921 beta (AHRR)")
ax[0].set_title("A  Control A: smoking\ndirection correct, p = "
                f"{res['control_A_AHRR_smoking']['p']:.3f}", fontsize=10, color=NAVY, loc="left")

null = np.array(nul["null"])
ax[1].hist(null, bins=20, color="#c8d4e6", edgecolor=GREY, lw=.5)
ax[1].axvline(0.5, color=GREY, ls=":", lw=1.2)
ax[1].axvline(res["primary_tki_response"]["per_repeat_mean"], color=RED, lw=2)
ax[1].set_xlabel("AUROC"); ax[1].set_ylabel("permutations")
ax[1].set_title("B  Observed sits inside the null\np = "
                f"{add['like_for_like_permutation_p']:.2f} (100 permutations)",
                fontsize=10, color=NAVY, loc="left")

lo, hi = add["auroc_ci95"]; mdt = add["min_detectable_auroc_80pct_power"]
obs = res["primary_tki_response"]["auroc"]
ax[2].axvspan(0.5, mdt, color="#eeeeee", zorder=0)
ax[2].errorbar([obs], [1], xerr=[[obs - lo], [hi - obs]], fmt="o", color=NAVY,
               capsize=5, ms=8, lw=2)
ax[2].axvline(0.5, color=GREY, ls=":", lw=1.2)
ax[2].axvline(mdt, color=RED, ls="--", lw=1.5)
ax[2].text(mdt + .006, 1.28, f"minimum detectable\nAUROC = {mdt:.2f}", fontsize=8,
           color=RED, va="center")
ax[2].text(0.505, 0.72, "undetectable at n = 69", fontsize=8, color=GREY)
ax[2].set_xlim(0.38, 0.80); ax[2].set_ylim(0.6, 1.5); ax[2].set_yticks([])
ax[2].set_xlabel("AUROC")
ax[2].set_title(f"C  Inconclusive, not negative\n{obs:.3f}  95% CI [{lo:.2f}, {hi:.2f}]",
                fontsize=10, color=NAVY, loc="left")
for a in ax:
    a.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("figures/r07_panels.pdf")
plt.savefig("figures/r07_panels.png", dpi=150)
print("figures/r07_panels.{pdf,png}")
