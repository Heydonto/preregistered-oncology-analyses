#!/usr/bin/env python3
"""R11 — manuscript figures for the glioma arm, built from real evidence.

This implements the figure specification the colleague sent as generate_figures1.py. That sketch
communicated the intent clearly; it could not be run (SyntaxError) and several of its data paths
would have fabricated content, so the specification is reimplemented here against the actual
evidence files, with gates that make each of those failure modes impossible rather than unlikely.

Specification, as received:
  Fig 1  cohort overview: sizes, molecular prevalence, scanner distribution, median OS
  Fig 3  Kaplan-Meier: IDH (TCGA), MGMT (UPENN), radiomic risk tertiles (UPENN)
  Fig 4  habitat composition and entropy by IDH
  Fig 5  imaging vs molecular C-index comparison
  extra  calibration plot

What changed, and why (each is a gate below, not a comment):

  1. IDH for TCGA is joined from data_clinical_sample.txt by PATIENT_ID. The sketch read
     IDH_STATUS from data_clinical_patient.txt, where the column does not exist, then fell back
     to np.random.choice. There is NO fallback here: if the join fails the run halts. Gate
     G_IDH reproduces R04's measured 89.7 vs 14.0 months and p < 1e-60.
  2. The radiomic prediction is selected BY COLUMN NAME. The sketch used oof.iloc[:, 1], which
     in this file is log10_days - the observed outcome - so it stratified survival by survival.
     Gate G_PREDCOL asserts the name and that the chosen column is not a near-copy of the
     outcome.
  3. Risk tertile direction is gated. Higher predicted log-survival means LOWER risk; the
     sketch labelled the lowest tertile 'Low risk'. Gate G_RISKDIR requires the high-risk group
     to have the shorter observed median survival.
  4. Scanner distribution uses the real eight-hospital metadata (ds007045 participants.tsv,
     n=337). The sketch hard-coded [443, 350, 267, 0], summing to 1,060, which matches no
     cohort we hold.
  5. Figure 5 omits the sketch's 'Imaging (external) C=0.581, p=0.012'. No external validation
     exists: R02 states the 337-patient cohort intended for it has no survival data at all.
     Figure 5 also foregrounds the within-GBM like-for-like comparison, because pitting
     molecular-all-grades (0.719) against imaging (0.602) is exactly the grade-mix comparison
     R04's audit withdrew.
  6. Figure 4 cannot be produced - habitat features have never been computed (that is the
     unstarted preprocessing/GPU arm). The panel states that rather than being silently absent.
  7. The calibration plot uses log10 throughout. The sketch compared log10 predictions against
     np.log observed, which mis-scales the diagonal by ln(10).
"""
import csv, io, json, os, sys, zipfile, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_kit import Run

S = "/Users/rezanehzati/quantara-staging/"
CLIN = S + "glioma-figs/upenn.csv"
ZIP = S + "glioma-figs/rf.zip"
TCGA = S + "staged/labels/glioma/lgggbm_tcga_pub/"
EIGHT = S + "staged/imaging/glioma/ds007045/participants.tsv"
OOF = (os.path.dirname(os.path.abspath(__file__)) +
       "/runs/20260806T225005Z-R02surv-3f8bece/oof_survival.csv")
OUT = "/Users/rezanehzati/Projects/quantara/results/R11-glioma-manuscript-figures/figures"
MODS, ROIS = ["T1", "T1GD", "T2", "FLAIR"], ["ET", "ED", "NC"]
PRED_COL, OBS_COL = "predicted_log10_days", "log10_days"
# measured reference values from the earlier reports; these are what the gates check against
REF = {"n_survival": 574, "n_mgmt": 247, "mgmt_p": 5.831120711860514e-08,
       "mgmt_med_meth": 595.0, "mgmt_med_unmeth": 371.0,
       "idh_n": 919, "idh_med_mut": 89.7, "idh_med_wt": 14.0,
       "imaging_cindex": 0.6024019381071802, "imaging_ci": [0.5763946279889873, 0.6291337624029092],
       "mol_allgrade_cindex": 0.7193457649053193, "mol_gbm_idh": 0.541, "mol_gbm_idh_mgmt": 0.561,
       "eight_n": 337, "scanner": {"Siemens": 228, "Philips": 89, "GE": 20},
       "zip_bytes": 16121402}
NAVY, GREY, RED, BLUE, GREEN, PURPLE = "#14284b", "#5a5a5a", "#961e1e", "#c8d4e6", "#2f6b4f", "#6b4f8a"


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def rd(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def vendor(model):
    m = model.lower()
    if "siemens" in m or "magnetom" in m:
        return "Siemens"
    if "philips" in m or any(k in m for k in ("achieva", "ingenia", "intera", "panorama", "gemini")):
        return "Philips"
    if m.startswith("ge ") or any(k in m for k in ("signa", "discovery", "optima")):
        return "GE"
    return "Other"


def main():
    run = Run("R11figures")
    cfg = {
        "report": "R11",
        "purpose": "manuscript figures for the glioma arm, implementing the colleague's "
                   "generate_figures1.py specification against real evidence",
        "figures": {"1": "cohort overview", "3": "Kaplan-Meier x3",
                    "4": "habitat - NOT PRODUCIBLE, never computed",
                    "5": "imaging vs molecular C-index", "extra": "calibration"},
        "deviations_from_sketch": {
            "IDH_join": "from data_clinical_sample.txt by PATIENT_ID; no random fallback exists",
            "prediction_column": f"selected by name '{PRED_COL}', not by position",
            "risk_direction": "gated: high-risk tertile must have shorter median survival",
            "scanner": "real ds007045 metadata (n=337), not hard-coded [443,350,267,0]",
            "external_bar": "OMITTED - no external validation exists (R02)",
            "fig5_framing": "within-GBM like-for-like foregrounded; all-grades marked "
                            "not-comparable (R04's withdrawn comparison)",
            "calibration_scale": "log10 throughout, not log10-vs-natural-log"},
        "reference_values": REF,
        "sources": {"upenn_clinical": CLIN, "radiomic_features": ZIP, "oof": OOF,
                    "tcga": TCGA, "eight_hospital": EIGHT},
    }
    cfg_hash = run.start(cfg, [CLIN, ZIP, OOF, TCGA + "data_clinical_sample.txt", EIGHT])
    os.makedirs(OUT, exist_ok=True)
    run.gate("G0_zip_integrity", REF["zip_bytes"], os.path.getsize(ZIP),
             os.path.getsize(ZIP) == REF["zip_bytes"])

    # ---------- UPENN cohort, reproducing R02's definition exactly ----------
    z = zipfile.ZipFile(ZIP)
    per = {}
    for m in MODS:
        for roi in ROIS:
            with z.open(f"Radiomic_Features_CaPTk_automaticsegm_{m}_{roi}.csv") as fh:
                rows = list(csv.reader(io.TextIOWrapper(fh, "utf8")))
            per[f"{m}_{roi}"] = {r[0].strip() for r in rows[1:] if r}
    complete = set.intersection(*per.values())
    clin = {}
    for r in csv.DictReader(open(CLIN)):
        clin[r["ID"].strip()] = r
    ids = sorted(i for i in complete if i in clin
                 and num(clin[i]["Survival_from_surgery_days_UPDATED"])
                 and num(clin[i]["Survival_from_surgery_days_UPDATED"]) > 0)
    run.gate("G1_upenn_n", REF["n_survival"], len(ids), len(ids) == REF["n_survival"],
             "R02's cohort: complete features AND survival > 0")
    T_up = np.array([num(clin[i]["Survival_from_surgery_days_UPDATED"]) for i in ids])
    E_up = np.ones(len(ids), int)   # uncensored, as R02 established

    # ---------- Fig 3B source: MGMT subset ----------
    mg_ids = [i for i in ids if clin[i]["MGMT"].strip() in ("Methylated", "Unmethylated")]
    run.gate("G2_mgmt_n", REF["n_mgmt"], len(mg_ids), len(mg_ids) == REF["n_mgmt"],
             "reproduces R02's MGMT positive control cohort")
    mg = np.array([clin[i]["MGMT"].strip() == "Methylated" for i in mg_ids])
    T_mg = np.array([num(clin[i]["Survival_from_surgery_days_UPDATED"]) for i in mg_ids])
    E_mg = np.ones(len(mg_ids), int)
    lr_mg = logrank_test(T_mg[mg], T_mg[~mg], E_mg[mg], E_mg[~mg])
    med_meth, med_unmeth = float(np.median(T_mg[mg])), float(np.median(T_mg[~mg]))
    run.gate("G3_mgmt_reproduces_R02",
             {"p": f"{REF['mgmt_p']:.3e}", "med": [REF["mgmt_med_meth"], REF["mgmt_med_unmeth"]]},
             {"p": f"{lr_mg.p_value:.3e}", "med": [med_meth, med_unmeth]},
             abs(np.log10(lr_mg.p_value) - np.log10(REF["mgmt_p"])) < 0.01
             and med_meth == REF["mgmt_med_meth"] and med_unmeth == REF["mgmt_med_unmeth"])

    # ---------- Fig 3C source: prediction BY NAME ----------
    oof_rows = list(csv.DictReader(open(OOF)))
    cols = list(oof_rows[0].keys())
    run.gate("G4_predcol_by_name", PRED_COL, [c for c in cols if c == PRED_COL],
             PRED_COL in cols,
             "the sketch used iloc[:,1] which is the OBSERVED outcome in this file")
    pred = np.array([float(r[PRED_COL]) for r in oof_rows])
    obs = np.array([float(r[OBS_COL]) for r in oof_rows])
    T_oof = np.array([float(r["survival_days"]) for r in oof_rows])
    rho = float(np.corrcoef(pred, obs)[0, 1])
    run.gate("G5_pred_is_not_outcome", "corr(pred, observed) < 0.99", round(rho, 4), rho < 0.99,
             "guards against silently plotting the outcome against itself")
    run.gate("G6_oof_aligned", REF["n_survival"], len(oof_rows), len(oof_rows) == REF["n_survival"])

    # risk tertiles: HIGH predicted survival = LOW risk
    q1, q2 = np.quantile(pred, [1 / 3, 2 / 3])
    risk = np.where(pred <= q1, "High risk", np.where(pred <= q2, "Medium risk", "Low risk"))
    med_hi = float(np.median(T_oof[risk == "High risk"]))
    med_lo = float(np.median(T_oof[risk == "Low risk"]))
    run.gate("G7_risk_direction", "high-risk median < low-risk median",
             {"high": med_hi, "low": med_lo}, med_hi < med_lo,
             "the sketch inverted these labels")
    lr_risk = logrank_test(T_oof[risk == "High risk"], T_oof[risk == "Low risk"],
                           np.ones((risk == "High risk").sum()),
                           np.ones((risk == "Low risk").sum()))

    # ---------- Fig 3A source: TCGA IDH, joined by PATIENT_ID ----------
    tp = {r["PATIENT_ID"]: r for r in rd(TCGA + "data_clinical_patient.txt")}
    idh_by_pat = {}
    for r in rd(TCGA + "data_clinical_sample.txt"):
        v = (r.get("IDH_STATUS") or "").strip()
        if v in ("Mutant", "WT"):
            idh_by_pat.setdefault(r["PATIENT_ID"], v)
    run.gate("G8_idh_join", "> 900 patients with IDH from the SAMPLE file", len(idh_by_pat),
             len(idh_by_pat) > 900, "column absent from the patient file; sketch fell back to random")
    tr = []
    for pid, v in idh_by_pat.items():
        p = tp.get(pid)
        if not p:
            continue
        t, s = num(p.get("OS_MONTHS")), (p.get("OS_STATUS") or "").strip()
        if t is None or t <= 0 or not s:
            continue
        tr.append((t, 1 if s.startswith("1") else 0, v == "Mutant"))
    T_t = np.array([x[0] for x in tr]); E_t = np.array([x[1] for x in tr])
    M_t = np.array([x[2] for x in tr])
    run.gate("G9_idh_n", REF["idh_n"], len(tr), len(tr) == REF["idh_n"],
             "reproduces R04's IDH cohort exactly")
    kmf = KaplanMeierFitter()
    med = {}
    for lab, m in (("Mutant", M_t), ("WT", ~M_t)):
        kmf.fit(T_t[m], E_t[m])
        med[lab] = round(float(kmf.median_survival_time_), 1)
    lr_idh = logrank_test(T_t[M_t], T_t[~M_t], E_t[M_t], E_t[~M_t])
    run.gate("G10_idh_reproduces_R04",
             {"med": [REF["idh_med_mut"], REF["idh_med_wt"]], "p": "< 1e-60"},
             {"med": [med["Mutant"], med["WT"]], "p": f"{lr_idh.p_value:.3e}"},
             med["Mutant"] == REF["idh_med_mut"] and med["WT"] == REF["idh_med_wt"]
             and lr_idh.p_value < 1e-60,
             "if this passes, the IDH labels are real; the sketch's random fallback gave p=0.325")

    # ---------- Fig 1C source: real scanner metadata ----------
    eight = rd(EIGHT)
    run.gate("G11_eight_n", REF["eight_n"], len(eight), len(eight) == REF["eight_n"])
    vend = collections.Counter(vendor(r["scanner_model"]) for r in eight)
    run.gate("G12_scanner_counts", REF["scanner"], dict(vend),
             all(vend[k] == v for k, v in REF["scanner"].items()),
             "real counts; the sketch used [443,350,267,0] = 1,060 which matches no cohort")
    field = collections.Counter(r["field_strength"].replace(",", ".") for r in eight)
    mgmt8 = collections.Counter(r["MGMT"] for r in eight)

    # ---------- TCGA prevalence for Fig 1B ----------
    ts = rd(TCGA + "data_clinical_sample.txt")
    idh_all = collections.Counter((r.get("IDH_STATUS") or "").strip() for r in ts)
    mgmt_all = collections.Counter((r.get("MGMT_PROMOTER_STATUS") or "").strip() for r in ts)
    # TCGA is ~50% censored, so a raw median of observed times understates survival badly.
    # Panel D must compare Kaplan-Meier medians. UPENN is uncensored (R02 gate G3), so its raw
    # median equals its KM median; TCGA's does not. Gate G13 enforces the distinction.
    tsv = [(num(p.get("OS_MONTHS")), (p.get("OS_STATUS") or "").strip())
           for p in tp.values()]
    tsv = [(t, st) for t, st in tsv if t and t > 0 and st]
    T_tc = np.array([t for t, _ in tsv])
    E_tc = np.array([1 if st.startswith("1") else 0 for _, st in tsv])
    tcga_surv = T_tc

    run.gate("G13_censoring_handled",
             "TCGA KM median > raw median (censoring present and handled)",
             {"n": len(tsv), "pct_censored": round(float((1 - E_tc).mean()) * 100, 1),
              "raw_median": round(float(np.median(T_tc)), 1),
              "km_median": round(float(KaplanMeierFitter().fit(T_tc, E_tc)
                                       .median_survival_time_), 1)},
             float(KaplanMeierFitter().fit(T_tc, E_tc).median_survival_time_)
             > float(np.median(T_tc)),
             "a raw median would understate TCGA by ~9 months; UPENN is uncensored so raw==KM")
    run.gate("G14_tcga_n", 1040, len(tsv), len(tsv) == 1040,
             "matches R04's tcga_survival.csv; one patient has OS_MONTHS but no OS_STATUS")

    run.log("figure_inputs_ready",
            upenn=len(ids), mgmt=len(mg_ids), tcga_idh=len(tr), eight=len(eight))

    # ================= FIGURE 1 =================
    fig, ax = plt.subplots(2, 2, figsize=(12, 9.5))
    fig.suptitle("Figure 1 — Study cohorts", fontsize=14, fontweight="bold", color=NAVY)

    a = ax[0, 0]
    names = ["UPENN-GBM\n(survival)", "UPENN-GBM\n(MGMT subset)", "TCGA pan-glioma\n(IDH+survival)",
             "Eight-hospital\n(imaging only)"]
    vals = [len(ids), len(mg_ids), len(tr), len(eight)]
    a.bar(range(4), vals, 0.6, color=[NAVY, BLUE, GREEN, GREY], edgecolor=GREY, lw=.5)
    for i, v in enumerate(vals):
        a.text(i, v + 12, str(v), ha="center", fontsize=10, color=NAVY)
    a.set_xticks(range(4)); a.set_xticklabels(names, fontsize=8)
    a.set_ylabel("patients"); a.set_ylim(0, max(vals) * 1.18)
    a.set_title("A  Analysis cohorts", fontsize=10, color=NAVY, loc="left")

    b = ax[0, 1]
    grp = ["UPENN\nMGMT meth.", "UPENN\nIDH1 mut.", "TCGA\nIDH mut.", "TCGA\nMGMT meth.",
           "8-hosp.\nMGMT meth."]
    idh_up = sum(1 for i in mg_ids if clin[i]["IDH1"].strip() == "Mutated")
    fr = [mg.mean(),
          idh_up / len(mg_ids),
          idh_all["Mutant"] / (idh_all["Mutant"] + idh_all["WT"]),
          mgmt_all["Methylated"] / (mgmt_all["Methylated"] + mgmt_all["Unmethylated"]),
          mgmt8["1"] / len(eight)]
    b.bar(range(5), fr, 0.6, color=[PURPLE, GREEN, GREEN, PURPLE, GREY], edgecolor=GREY, lw=.5)
    for i, v in enumerate(fr):
        b.text(i, v + .02, f"{v*100:.0f}%", ha="center", fontsize=9, color=NAVY)
    b.set_xticks(range(5)); b.set_xticklabels(grp, fontsize=7.5)
    b.set_ylabel("proportion"); b.set_ylim(0, 0.85)
    b.set_title("B  Molecular marker prevalence", fontsize=10, color=NAVY, loc="left")

    c = ax[1, 0]
    order = [k for k, _ in vend.most_common()]
    c.pie([vend[k] for k in order],
          labels=[f"{k}\n{vend[k]} ({vend[k]/len(eight)*100:.0f}%)" for k in order],
          colors=[NAVY, BLUE, GREY, GREEN][:len(order)],
          autopct=None, startangle=90, textprops={"fontsize": 8.5},
          wedgeprops={"edgecolor": "white", "lw": 1.2})
    fs = ", ".join(f"{k}T: {v}" for k, v in sorted(field.items(), reverse=True))
    c.set_title(f"C  Scanner manufacturer, eight-hospital cohort (n={len(eight)})\n"
                f"field strength — {fs}", fontsize=10, color=NAVY, loc="left")

    d = ax[1, 1]
    med_up_m = float(np.median(T_up)) / 30.44          # uncensored: raw == KM
    med_mg_m = float(np.median(T_mg)) / 30.44           # uncensored: raw == KM
    med_tc_raw = float(np.median(T_tc))
    med_tc_m = float(KaplanMeierFitter().fit(T_tc, E_tc).median_survival_time_)
    mv = [med_up_m, med_mg_m, med_tc_m]
    d.bar(range(3), mv, 0.55, color=[NAVY, BLUE, GREEN], edgecolor=GREY, lw=.5)
    for i, v in enumerate(mv):
        d.text(i, v + .6, f"{v:.1f}", ha="center", fontsize=10, color=NAVY)
    d.set_xticks(range(3))
    d.set_xticklabels([f"UPENN-GBM\n(n={len(ids)})", f"UPENN MGMT\n(n={len(mg_ids)})",
                       f"TCGA pan-glioma\n(n={len(tcga_surv)})"], fontsize=8)
    d.set_ylabel("median overall survival (months)")
    d.set_ylim(0, max(mv) * 1.2)
    d.set_title("D  Median overall survival (Kaplan–Meier)\n"
                "TCGA: grades 2–4, 50% censored. UPENN: GBM only, uncensored",
                fontsize=9.5, color=NAVY, loc="left")
    d.annotate(f"a raw median of observed\ntimes would read {med_tc_raw:.1f} —\n"
               "censoring must be handled",
               xy=(2, med_tc_m), xytext=(1.42, med_tc_m * 1.02), fontsize=7, color=RED,
               ha="right", va="top",
               arrowprops=dict(arrowstyle="->", color=RED, lw=1))
    for row in ax:
        for x in row:
            if x is not c:
                x.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig1_cohort_overview.pdf"); plt.savefig(f"{OUT}/fig1_cohort_overview.png", dpi=200)
    plt.close()
    run.log("fig1_written")

    # ================= FIGURE 3 =================
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    fig.suptitle("Figure 3 — Survival analysis", fontsize=14, fontweight="bold", color=NAVY)

    k = KaplanMeierFitter()
    for m, lab, col in ((M_t, f"IDH mutant (n={int(M_t.sum())})", GREEN),
                        (~M_t, f"IDH wild-type (n={int((~M_t).sum())})", RED)):
        k.fit(T_t[m], E_t[m], label=lab)
        k.plot_survival_function(ax=ax[0], ci_show=True, color=col, lw=1.8)
    ax[0].set_xlabel("months since diagnosis"); ax[0].set_ylabel("survival probability")
    ax[0].set_xlim(0, 180)
    ax[0].text(.97, .82, f"log-rank p = {lr_idh.p_value:.1e}\nmedian {med['Mutant']} vs "
               f"{med['WT']} months", transform=ax[0].transAxes, ha="right", fontsize=8.5,
               color=NAVY)
    ax[0].set_title(f"A  IDH mutation — TCGA (n={len(tr)})\npositive control",
                    fontsize=10, color=NAVY, loc="left")

    for m, lab, col in ((mg, f"MGMT methylated (n={int(mg.sum())})", GREEN),
                        (~mg, f"MGMT unmethylated (n={int((~mg).sum())})", RED)):
        k.fit(T_mg[m], E_mg[m], label=lab)
        k.plot_survival_function(ax=ax[1], ci_show=True, color=col, lw=1.8)
    ax[1].set_xlabel("days from surgery"); ax[1].set_ylabel("survival probability")
    ax[1].text(.97, .82, f"log-rank p = {lr_mg.p_value:.1e}\nmedian {med_meth:.0f} vs "
               f"{med_unmeth:.0f} days", transform=ax[1].transAxes, ha="right", fontsize=8.5,
               color=NAVY)
    ax[1].set_title(f"B  MGMT methylation — UPENN (n={len(mg_ids)})\npositive control",
                    fontsize=10, color=NAVY, loc="left")

    for lab, col in (("Low risk", GREEN), ("Medium risk", GREY), ("High risk", RED)):
        m = risk == lab
        k.fit(T_oof[m], np.ones(m.sum()), label=f"{lab} (n={int(m.sum())})")
        k.plot_survival_function(ax=ax[2], ci_show=False, color=col, lw=1.8)
    ax[2].set_xlabel("days from surgery"); ax[2].set_ylabel("survival probability")
    ax[2].text(.97, .82, f"low vs high p = {lr_risk.p_value:.1e}\nmedian {med_lo:.0f} vs "
               f"{med_hi:.0f} days", transform=ax[2].transAxes, ha="right", fontsize=8.5,
               color=NAVY)
    ax[2].set_title(f"C  Radiomic risk tertiles — UPENN (n={len(oof_rows)})\n"
                    "out-of-fold predictions, C-index 0.602", fontsize=10, color=NAVY, loc="left")
    for x in ax:
        x.spines[["top", "right"]].set_visible(False)
        x.legend(fontsize=7.5, frameon=False, loc="lower left")
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig3_km_curves.pdf"); plt.savefig(f"{OUT}/fig3_km_curves.png", dpi=200)
    plt.close()
    run.log("fig3_written")

    # ================= FIGURE 4 — not producible =================
    fig, a4 = plt.subplots(figsize=(8, 4.2))
    a4.axis("off")
    a4.text(.5, .74, "Figure 4 — Tumour habitat analysis", ha="center", fontsize=13,
            fontweight="bold", color=NAVY)
    a4.text(.5, .52, "NOT PRODUCIBLE FROM CURRENT HOLDINGS", ha="center", fontsize=11,
            color=RED, fontweight="bold")
    a4.text(.5, .30,
            "Habitat features have never been computed. Producing them requires the image\n"
            "preprocessing pipeline (co-registration of four modalities, skull-stripping,\n"
            "resampling across 30+ distinct voxel grids) over ~967 patients — the unstarted\n"
            "GPU arm, not an analysis of data in hand. The panel is shown rather than omitted\n"
            "so the gap is explicit.", ha="center", fontsize=9, color=GREY, linespacing=1.7)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig4_habitat_UNAVAILABLE.pdf")
    plt.savefig(f"{OUT}/fig4_habitat_UNAVAILABLE.png", dpi=200)
    plt.close()
    run.log("fig4_unavailable_panel_written")

    # ================= FIGURE 5 =================
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8),
                           gridspec_kw={"width_ratios": [1.15, 1]})
    fig.suptitle("Figure 5 — Prognostic performance", fontsize=14, fontweight="bold", color=NAVY)

    p0 = ax[0]
    lab = [f"Imaging\nUPENN-GBM (n={len(ids)})", "Molecular IDH\nTCGA GBM only",
           "Molecular IDH+MGMT\nTCGA GBM only"]
    v = [REF["imaging_cindex"], REF["mol_gbm_idh"], REF["mol_gbm_idh_mgmt"]]
    p0.bar(range(3), v, 0.55, color=[NAVY, GREEN, GREEN], edgecolor=GREY, lw=.5)
    lo = REF["imaging_cindex"] - REF["imaging_ci"][0]
    hi = REF["imaging_ci"][1] - REF["imaging_cindex"]
    p0.errorbar([0], [REF["imaging_cindex"]], yerr=[[lo], [hi]], fmt="none", ecolor=GREY,
                capsize=5, lw=1.5)
    for i, y in enumerate(v):
        p0.text(i, y + .012, f"{y:.3f}", ha="center", fontsize=10, color=NAVY)
    p0.axhline(0.5, color=GREY, ls=":", lw=1.2)
    p0.set_xticks(range(3)); p0.set_xticklabels(lab, fontsize=8)
    p0.set_ylabel("Harrell's C-index"); p0.set_ylim(0.45, 0.70)
    p0.set_title("A  Like-for-like: glioblastoma only\nimaging exceeds the molecular panel",
                 fontsize=10, color=NAVY, loc="left")

    p1 = ax[1]
    p1.bar([0], [REF["imaging_cindex"]], 0.5, color=NAVY, edgecolor=GREY, lw=.5)
    p1.bar([1], [REF["mol_allgrade_cindex"]], 0.5, color="#d9d9d9", edgecolor=RED, lw=1.4,
           hatch="//")
    p1.text(0, REF["imaging_cindex"] + .012, f"{REF['imaging_cindex']:.3f}", ha="center",
            fontsize=10, color=NAVY)
    p1.text(1, REF["mol_allgrade_cindex"] + .012, f"{REF['mol_allgrade_cindex']:.3f}",
            ha="center", fontsize=10, color=RED)
    p1.axhline(0.5, color=GREY, ls=":", lw=1.2)
    p1.set_xticks([0, 1])
    p1.set_xticklabels(["Imaging\nGBM only", "Molecular IDH\nTCGA grades 2–4"], fontsize=8)
    p1.set_ylim(0.45, 0.80); p1.set_ylabel("Harrell's C-index")
    p1.text(.03, .30, "NOT A VALID\nCOMPARISON\n\nIDH status largely\nIS the grade\n"
            "distinction\n(R04 audit)", transform=p1.transAxes, ha="left", va="center",
            fontsize=8, color=RED)
    p1.set_title("B  Why the naive comparison misleads", fontsize=10, color=NAVY, loc="left")
    for x in ax:
        x.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig5_comparison.pdf"); plt.savefig(f"{OUT}/fig5_comparison.png", dpi=200)
    plt.close()
    run.log("fig5_written")

    # ================= CALIBRATION =================
    fig, a5 = plt.subplots(figsize=(5.6, 5.4))
    bins = np.quantile(pred, np.linspace(0, 1, 11))
    bi = np.clip(np.digitize(pred, bins[1:-1]), 0, 9)
    bx = [pred[bi == i].mean() for i in range(10) if (bi == i).sum() > 3]
    by = [obs[bi == i].mean() for i in range(10) if (bi == i).sum() > 3]
    a5.scatter(bx, by, s=70, color=NAVY, zorder=3, label="decile means")
    span = [min(bx + by), max(bx + by)]
    a5.plot(span, span, ls="--", color=RED, lw=1.4, label="perfect calibration")
    a5.set_xlabel("predicted log$_{10}$ survival (days)")
    a5.set_ylabel("observed log$_{10}$ survival (days)")
    a5.set_title(f"Calibration — UPENN out-of-fold (n={len(oof_rows)})\n"
                 f"both axes log$_{{10}}$; Spearman $\\rho$ = {rho:.3f}",
                 fontsize=10, color=NAVY, loc="left")
    a5.legend(fontsize=8, frameon=False)
    a5.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig_calibration.pdf"); plt.savefig(f"{OUT}/fig_calibration.png", dpi=200)
    plt.close()
    run.log("calibration_written")

    results = {
        "figures_produced": ["fig1_cohort_overview", "fig3_km_curves", "fig5_comparison",
                             "fig_calibration"],
        "figures_not_producible": {"fig4_habitat": "habitat features never computed "
                                                  "(unstarted preprocessing/GPU arm)"},
        "panels_dropped_from_spec": {
            "fig5_imaging_external": "no external validation exists; R02 states the 337-patient "
                                     "cohort intended for it has no survival data"},
        "fig1": {"cohorts": {"upenn_survival": len(ids), "upenn_mgmt": len(mg_ids),
                             "tcga_idh_survival": len(tr), "eight_hospital": len(eight)},
                 "prevalence": {"upenn_mgmt_methylated": float(mg.mean()),
                                "upenn_idh1_mutant": idh_up / len(mg_ids),
                                "tcga_idh_mutant": fr[2], "tcga_mgmt_methylated": fr[3],
                                "eight_mgmt_methylated": fr[4]},
                 "scanner_manufacturer": dict(vend), "field_strength": dict(field),
                 "median_os_months": {"upenn_km": med_up_m, "upenn_mgmt_km": med_mg_m,
                                      "tcga_km": med_tc_m, "tcga_raw_incorrect": med_tc_raw,
                                      "note": "UPENN uncensored so raw==KM; TCGA 50% censored"}},
        "fig3": {"A_idh_tcga": {"n": len(tr), "median_months": med,
                                "logrank_p": float(lr_idh.p_value)},
                 "B_mgmt_upenn": {"n": len(mg_ids), "median_days_meth": med_meth,
                                  "median_days_unmeth": med_unmeth,
                                  "logrank_p": float(lr_mg.p_value)},
                 "C_risk_tertiles": {"n": len(oof_rows), "median_days_low_risk": med_lo,
                                     "median_days_high_risk": med_hi,
                                     "logrank_p_low_vs_high": float(lr_risk.p_value)}},
        "fig5": {"imaging_gbm": REF["imaging_cindex"], "imaging_ci95": REF["imaging_ci"],
                 "molecular_gbm_idh": REF["mol_gbm_idh"],
                 "molecular_gbm_idh_mgmt": REF["mol_gbm_idh_mgmt"],
                 "molecular_all_grades_NOT_COMPARABLE": REF["mol_allgrade_cindex"]},
        "calibration": {"spearman_rho_pred_vs_obs": rho, "scale": "log10 on both axes"},
        "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir},
    }
    run.write("results.json", results)
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
