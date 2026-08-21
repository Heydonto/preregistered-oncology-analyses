#!/usr/bin/env python3
"""R16 figure — reads evidence/ only."""
import json
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

E = "evidence/"
r = json.load(open(E + "results.json"))
df = pd.read_csv(E + "analysis_frame.csv")
zs = np.load(E + "permutation_z.npy")
NAVY, GREY, RED, BLUE, GREEN = "#14284b", "#5a5a5a", "#961e1e", "#c8d4e6", "#2f6b4f"

fig, ax = plt.subplots(2, 2, figsize=(12.4, 8.6))

# ---- A: KM by the plan's FIXED 3.5 threshold (S3) ----
a = ax[0, 0]
k = KaplanMeierFitter()
for hi, lab, col in ((0, "post-RT peak SUV < 3.5", GREEN), (1, "≥ 3.5", RED)):
    m = (df["suv"] >= 3.5).astype(int) == hi
    k.fit(df["T_months"][m], df["event"][m], label=f"{lab} (n={int(m.sum())})")
    k.plot_survival_function(ax=a, ci_show=True, color=col, lw=1.8)
s3 = r["secondary_tests"]["S3_fixed_threshold_3.5"]
a.set_xlabel("months from post-treatment PET (landmark)")
a.set_ylabel("overall survival")
a.text(.97, .84, f"HR {s3['hr']:.2f}\nHolm p = {s3['holm_p']:.4f}", transform=a.transAxes,
       ha="right", fontsize=9, color=NAVY)
a.legend(fontsize=7.5, frameon=False, loc="lower left")
a.set_title("A  Pre-specified threshold, no cutpoint search\n"
            "3.5 was fixed by the plan before any outcome was read",
            fontsize=10, color=NAVY, loc="left")

# ---- B: permutation null vs observed ----
b = ax[0, 1]
p = r["PRIMARY"]
b.hist(zs, bins=40, color=BLUE, edgecolor=GREY, lw=.4)
b.axvline(p["z"], color=RED, lw=2.2)
b.axvline(0, color=GREY, ls=":", lw=1.2)
b.set_xlabel("Cox score z under exposure permutation")
b.set_ylabel(f"permutations (B={len(zs)})")
b.text(.03, .93, f"observed z = {p['z']:.2f}\npermutation p = {p['permutation_p']:.4f}\n"
       f"null mean {zs.mean():.3f}, SD {zs.std(ddof=1):.3f}",
       transform=b.transAxes, fontsize=8.5, color=NAVY, va="top")
b.set_title("B  The p-value of record is the permutation p\nnull correctly centred and calibrated",
            fontsize=10, color=NAVY, loc="left")

# ---- C: primary and all seven secondaries ----
c = ax[1, 0]
S = r["secondary_tests"]
rows = [("PRIMARY  log2 peak SUV", p["hr_per_doubling"], p["ci95"], True)] + [
    (k2.replace("_", " "), v["hr"], None, False) for k2, v in S.items()]
y = np.arange(len(rows))[::-1]
for i, (lab, hr, ci, isp) in enumerate(rows):
    col = NAVY if isp else GREY
    c.plot([hr], [y[i]], "o" if isp else "s", color=col, ms=8 if isp else 6)
    if ci:
        c.plot(ci, [y[i], y[i]], color=col, lw=1.6)
    c.text(hr + .02, y[i], f"{hr:.3f}", va="center", fontsize=8, color=col)
c.axvline(1.0, color=RED, ls="--", lw=1.3)
c.set_yticks(y); c.set_yticklabels([r0[0] for r0 in rows], fontsize=8)
c.set_xlabel("hazard ratio"); c.set_xlim(0.95, 1.95)
c.set_title("C  Every pre-declared robustness test agrees\nall seven Holm-significant",
            fontsize=10, color=NAVY, loc="left")

# ---- D: the exposure is 59% floored ----
d = ax[1, 1]
suv = df["suv"].to_numpy()
d.hist(np.clip(suv, 0, 12), bins=40, color=BLUE, edgecolor=GREY, lw=.4)
nfl = r["cohort"]["floored_at_0.5"]
d.axvline(3.5, color=RED, ls="--", lw=1.4)
d.text(3.7, d.get_ylim()[1]*.55, "plan's fixed\nthreshold 3.5", fontsize=8, color=RED)
d.set_xlabel("post-RT peak SUV (LE.lee11), clipped at 12 for display")
d.set_ylabel("patients")
d.text(.40, .97, f"{nfl} of 180 exposure values ≤ 0.5\n(complete metabolic response),\n"
       "floored by the plan's rule", transform=d.transAxes, fontsize=8.5, color=NAVY, va="top")
d.set_title("D  The caveat the plan forced us to report\nthe effect rests on a minority with residual uptake",
            fontsize=10, color=NAVY, loc="left")

for row in ax:
    for x_ in row:
        x_.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("figures/r16_panels.pdf"); plt.savefig("figures/r16_panels.png", dpi=150)
print("figures/r16_panels.{pdf,png}")
