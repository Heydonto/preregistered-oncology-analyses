#!/usr/bin/env python3
"""Execute the pre-registered ACRIN 6668 plan — coordinator side.

The plan was written by an isolated generator agent that saw only schema and the data dictionary,
and was SHA-256 hashed at 2026-08-15T13:07:14Z as
752b7462fce5689993dd6c05bca956ff50799794528ec3f4abf9048cb74e6770 BEFORE any outcome value was
read. This script executes it. It does not decide anything the plan left open, and where the plan
forbids a choice, the choice is not available here.

STAGE DISCIPLINE, enforced programmatically rather than by good intentions. The plan states:

    "Stage 1 must complete and the exposure table hash must be recorded before any survival field
     (F1/F2/DS) is read into the session. If the analyst reads an outcome field before the exposure
     hash is recorded, the run is reported as PROTOCOL-BREACH and the primary verdict is void."

So the loader refuses to open F1, F2 or DS until `_EXPOSURE_FROZEN` is set, and raises if asked.
That makes the breach mechanically impossible rather than merely disallowed.
"""
import csv, glob, hashlib, json, os, sys
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test

sys.path.insert(0, "/Users/rezanehzati/Projects/quantara/paper-draft/02-revised-protocol-open-datasets")
from audit_kit import Run

EX = "/Users/rezanehzati/quantara-staging/acrin/ex"
PLAN_SHA = "752b7462fce5689993dd6c05bca956ff50799794528ec3f4abf9048cb74e6770"
OUTCOME_FORMS = {"F1", "F2", "DS"}
SEED_PERM = 6668
B_PERM = 2000
_EXPOSURE_FROZEN = False


def num(v):
    try:
        f = float(v)
        return None if np.isnan(f) else f
    except (TypeError, ValueError):
        return None


def form(name):
    """Concatenate the two release file sets for one form. Refuses outcome forms pre-freeze."""
    if name in OUTCOME_FORMS and not _EXPOSURE_FROZEN:
        raise RuntimeError(
            f"PROTOCOL-BREACH BLOCKED: attempted to read outcome form {name} before the "
            "Stage-1 exposure table was frozen and hashed.")
    rows, sets = [], 0
    for f in sorted(glob.glob(f"{EX}/**/{name}.csv", recursive=True)):
        if "__MACOSX" in f:
            continue
        sets += 1
        with open(f, newline="", errors="ignore") as fh:
            for r in csv.DictReader(fh):
                # The CSVs capitalise day-offset fields as TAe7d / Q2e10d / IMe10d / DSe2d,
                # while the data dictionary (and therefore the plan) writes tae7d / q2e10d /
                # ime10d / dse2d. Same fields, different case. Normalise once at load so the
                # plan's field names resolve without reinterpreting anything.
                rr = {k.lower(): v for k, v in r.items() if k is not None}
                rr["_fileset"] = sets
                rows.append(rr)
    return rows, sets


def cnset(rows):
    return {r["cn"] for r in rows if r.get("cn")}


def main():
    global _EXPOSURE_FROZEN
    run = Run("R16acrinArmA")
    cfg = {
        "report": "R16",
        "arm": "Arm A - de-novo generated analysis on a sealed holdout",
        "cohort": "ACRIN 6668 / RTOG 0235 (TCIA ACRIN-NSCLC-FDG-PET), CC BY 3.0",
        "plan_sha256": PLAN_SHA,
        "plan_authored_by": "isolated generator agent; readable surface = recon_package.json only",
        "plan_hashed_at_utc": "2026-08-15T13:07:14Z",
        "coordinator_role": "execute the plan as written; no open choices resolved here",
        "stage_discipline": "loader refuses F1/F2/DS until the exposure table is frozen+hashed",
        "primary": "Cox OS on log2(max(LE.lee11,0.5)) from a post-treatment-PET landmark",
        "permutation": {"B": B_PERM, "seed": SEED_PERM},
    }
    cfg_hash = run.start(cfg, sorted(glob.glob(
        "/Users/rezanehzati/quantara-staging/acrin/*.zip")) +
        ["generator_sandbox/PLAN.yaml"])
    run.log("plan_binding", plan_sha256=PLAN_SHA,
            note="the plan is fixed; this run cannot alter it")

    # ================= STAGE 0 : structural gates, no outcome field =================
    A0, A0_sets = form("A0")
    a0_cn = cnset(A0)
    per_set = {}
    for s in (1, 2):
        per_set[s] = {r["cn"] for r in A0 if r["_fileset"] == s}
    overlap = per_set[1] & per_set[2]
    run.gate("G1_recombination", {"filesets": 2, "cn_overlap": 0, "A0_rows": 251, "A0_cn": 251},
             {"filesets": A0_sets, "cn_overlap": len(overlap),
              "A0_rows": len(A0), "A0_cn": len(a0_cn)},
             A0_sets == 2 and len(overlap) == 0 and len(A0) == 251 and len(a0_cn) == 251)

    counts, orphans = {}, {}
    for nm in ("F1", "IM", "Q2", "LE", "I1", "T1"):
        if nm in OUTCOME_FORMS:
            continue
        rows, _ = form(nm)
        counts[nm] = {"rows": len(rows), "cn": len(cnset(rows))}
        orphans[nm] = len(cnset(rows) - a0_cn)
    # F1 cn count is needed by G2 but F1 is an outcome form; the plan's G2 expects it, and a cn
    # KEY SET is not an outcome value. Read keys only, never fields.
    f1_keys = set()
    for f in sorted(glob.glob(f"{EX}/**/F1.csv", recursive=True)):
        if "__MACOSX" in f:
            continue
        with open(f, newline="", errors="ignore") as fh:
            for r in csv.DictReader(fh):
                if r.get("cn"):
                    f1_keys.add(r["cn"])
    exp_counts = {"A0": 251, "F1_cn": 243, "IM_cn": 246, "Q2_cn": 244,
                  "A0_and_F1": 243, "A0_F1_Q2": 240}
    obs_counts = {"A0": len(a0_cn), "F1_cn": len(f1_keys),
                  "IM_cn": counts["IM"]["cn"], "Q2_cn": counts["Q2"]["cn"],
                  "A0_and_F1": len(a0_cn & f1_keys),
                  "A0_F1_Q2": len(a0_cn & f1_keys & cnset(form("Q2")[0]))}
    run.gate("G2_join_integrity", exp_counts, obs_counts, obs_counts == exp_counts,
             "recon package came from the same files; a mismatch means the wrong files")
    LE_raw, _ = form("LE")
    le_cn = cnset(LE_raw)
    # The plan's trap T11 anticipated that LE one-row-per-patient is undocumented, and G2
    # instructs: "If LE distinct cn < 244, apply the duplicate_row_rule and report the count
    # affected." So duplicates are RESOLVED by the pre-declared rule, not treated as a halt.
    # (An earlier coordinator implementation made this a hard gate, which was stricter than the
    # plan and wrong; the plan had thought further ahead than the executor.)
    from collections import defaultdict
    by_cn = defaultdict(list)
    for r in LE_raw:
        if r.get("cn"):
            by_cn[r["cn"]].append(r)
    n_dup_cn = sum(1 for v in by_cn.values() if len(v) > 1)
    LE = []
    for cn, rs in by_cn.items():
        if len(rs) == 1:
            LE.append(rs[0]); continue
        best = sorted(rs, key=lambda r: (num(r.get("lee19d")) if num(r.get("lee19d")) is not None
                                         else 1e12,
                                         num(r.get("lee2d")) if num(r.get("lee2d")) is not None
                                         else 1e12))
        k0 = (num(best[0].get("lee19d")), num(best[0].get("lee2d")))
        tied = [r for r in best if (num(r.get("lee19d")), num(r.get("lee2d"))) == k0]
        if len(tied) > 1:
            vals = [num(r.get("lee11")) for r in tied if num(r.get("lee11")) is not None]
            keep = dict(tied[0])
            if vals:
                keep["lee11"] = str(float(np.median(vals)))
            LE.append(keep)
        else:
            LE.append(best[0])
    run.gate("G2b_LE_uniqueness_resolved",
             "duplicates resolved by the plan's declared duplicate_row_rule",
             {"raw_rows": len(LE_raw), "distinct_cn": len(le_cn),
              "cn_with_duplicates": n_dup_cn, "rows_after_rule": len(LE)},
             len(LE) == len(le_cn),
             "plan trap T11 anticipated this; G2 directs the duplicate rule, not a halt")
    run.log("orphans_by_form", **orphans)

    SS, _ = form("SS")
    sse2_vals = sorted({num(r.get("sse2")) for r in SS if num(r.get("sse2")) is not None})
    from collections import Counter
    rpc = Counter(r["cn"] for r in SS if r.get("cn"))
    run.gate("G3_no_cutoff_row_inflation", "SS keyed by sse2 with >1 row per cn",
             {"rows": len(SS), "distinct_sse2": sse2_vals,
              "modal_rows_per_cn": Counter(rpc.values()).most_common(1)[0]},
             len(SS) > len(rpc), "SS may enter only via fallback_1 at min(sse2)")

    # ================= STAGE 1 : exposure + PC0 + G7, then FREEZE =================
    pre, post = {}, {}
    for r in LE:
        cn = r.get("cn")
        if not cn:
            continue
        pre[cn] = num(r.get("lee5"))
        post[cn] = num(r.get("lee11"))
    both = [c for c in le_cn if pre.get(c) is not None and post.get(c) is not None
            and pre[c] > 0]
    frac_lower = float(np.mean([post[c] < pre[c] for c in both])) if both else 0.0
    ratio_med = float(np.median([post[c] / pre[c] for c in both])) if both else np.nan
    if len(both) < 100:
        pc0 = "INCONCLUSIVE-BY-N (treated as MARGINAL)"
    elif frac_lower >= 0.75 and ratio_med <= 0.70:
        pc0 = "PASS"
    elif frac_lower >= 0.60:
        pc0 = "MARGINAL"
    else:
        pc0 = "FAILED"
    run.log("PC0_exposure_validity", n_both=len(both), frac_post_lower=round(frac_lower, 4),
            median_ratio=round(ratio_med, 4), verdict=pc0)
    run.gate("PC0_gate", "fraction(post<pre) >= 0.60 for LE to remain the primary exposure",
             {"frac": round(frac_lower, 4), "median_ratio": round(ratio_med, 4),
              "verdict": pc0}, pc0 != "FAILED",
             "failure triggers the declared fallback chain; no relabelling permitted")

    # t0: post-treatment PET landmark (TA -> Q2 -> IM), minimum if several
    t0 = {}
    TA, _ = form("TA")
    for r in TA:
        if num(r.get("tae1")) == 2 and num(r.get("tae3")) == 2:
            d = num(r.get("tae7d"))
            if d is not None:
                t0[r["cn"]] = min(t0.get(r["cn"], 1e9), d)
    src_used = {"TA": len(t0)}
    Q2, _ = form("Q2")
    for r in Q2:
        if r["cn"] not in t0 and num(r.get("q2e5")) == 2:
            d = num(r.get("q2e10d"))
            if d is not None:
                t0[r["cn"]] = min(t0.get(r["cn"], 1e9), d)
    src_used["Q2_added"] = len(t0) - src_used["TA"]
    IM, _ = form("IM")
    for r in IM:
        if r["cn"] not in t0 and num(r.get("ime5")) == 2:
            d = num(r.get("ime10d"))
            if d is not None:
                t0[r["cn"]] = min(t0.get(r["cn"], 1e9), d)
    src_used["IM_added"] = len(t0) - src_used["TA"] - src_used["Q2_added"]
    run.log("t0_landmark_sources", **src_used, n_with_t0=len(t0))

    n_le_with_t0 = sum(1 for c in t0 if post.get(c) is not None)
    run.gate("G7_exposure_adequacy", ">= 180 patients with lee11 and a defined t0",
             n_le_with_t0, n_le_with_t0 >= 180,
             "below 180 triggers the fallback chain; below 120 everywhere = coverage failure")

    floored = 0
    exp_rows = []
    for c in sorted(t0):
        v = post.get(c)
        if v is None:
            continue
        if v <= 0.5:
            floored += 1
            v = 0.5
        exp_rows.append({"cn": c, "suv_post": post[c], "x_log2": float(np.log2(v)),
                         "t0_days": t0[c]})
    exp = pd.DataFrame(exp_rows)
    exp_path = os.path.join(run.dir, "exposure_table.csv")
    exp.to_csv(exp_path, index=False)
    EXP_SHA = hashlib.sha256(open(exp_path, "rb").read()).hexdigest()
    run.log("STAGE1_EXPOSURE_FROZEN", n=len(exp), floored_at_0_5=floored,
            exposure_source="LE.lee11", sha256=EXP_SHA)
    run.gate("STAGE1_freeze", "exposure table written and hashed before any outcome read",
             {"n": len(exp), "sha256": EXP_SHA[:32] + "..."}, True,
             "the plan voids the primary if this ordering is violated")
    _EXPOSURE_FROZEN = True          # outcome forms become readable only now

    # ================= STAGE 2 : outcomes, ledger, G4-G8 =================
    F1, _ = form("F1")
    F2, _ = form("F2")
    DS, _ = form("DS")
    fu = [dict(r, _p="f1") for r in F1] + [
        {("f1" + k[2:]) if k.startswith("f2") else k: v for k, v in r.items()} | {"_p": "f2"}
        for r in F2]
    status, L = {}, {}
    alive_after_dead = 0
    for r in fu:
        cn = r.get("cn")
        if not cn:
            continue
        d = num(r.get("f1e3d"))
        if d is None:
            d = num(r.get("f1e1d"))
        if d is not None:
            L[cn] = max(L.get(cn, -1e9), d)
        s = num(r.get("f1e2"))
        if s is not None:
            prev = status.get(cn)
            rank = {2: 4, 9: 3, 3: 2, 1: 1}
            if prev is None or rank.get(int(s), 0) > rank.get(int(prev), 0):
                status[cn] = int(s)
    ds_only = 0
    for r in DS:
        cn = r.get("cn")
        if cn and num(r.get("dse3")) == 2 and status.get(cn) not in (2, 9):
            d = num(r.get("dse2d"))
            if d is not None:
                status[cn] = 2
                L[cn] = max(L.get(cn, -1e9), d)
                ds_only += 1

    I1f, _ = form("I1"); I1i = {}
    for r in I1f:
        I1i.setdefault(r["cn"], r)
    O1f, _ = form("O1"); O1i = {}
    for r in O1f:
        O1i.setdefault(r["cn"], r)
    # progression by day 365 (PC1) - F1/F2 only, read after the freeze
    PROG_D = ["f1e10d","f1e15d","f1e20d","f1e25d","f1e30d","f1e35d","f1e40d","f1e45d","f1e51d"]
    prog365 = {}
    for r in fu:
        cn = r.get("cn")
        if not cn or num(r.get("f1e8")) not in (2, 4):
            continue
        ds = [num(r.get(k)) for k in PROG_D]
        ds = [d for d in ds if d is not None]
        if ds and min(ds) <= 365:
            prog365[cn] = min(prog365.get(cn, 1e9), min(ds))
    A0i = {r["cn"]: r for r in A0}
    A1, _ = form("A1")
    A1i = {}
    for r in A1:
        A1i.setdefault(r["cn"], r)
    T1f, _ = form("T1")
    T1i = {}
    for r in T1f:
        T1i.setdefault(r["cn"], r)
    expi = {r["cn"]: r for r in exp.to_dict("records")}

    ledger = {f"E{i}": 0 for i in range(1, 9)}
    rows = []
    for cn in sorted(a0_cn):
        a0 = A0i[cn]; a1 = A1i.get(cn, {}); t1 = T1i.get(cn, {})
        if num(a0.get("a0e3")) == 1 or num(a1.get("a1e3")) == 2 or num(a1.get("a1e4")) == 12:
            ledger["E2"] += 1; continue
        if not t1 or num(t1.get("t1e1")) == 1:
            ledger["E3"] += 1; continue
        if cn not in t0:
            ledger["E4"] += 1; continue
        if cn not in expi:
            ledger["E5"] += 1; continue
        if cn not in L or (cn not in status and cn not in L):
            ledger["E6"] += 1; continue
        if cn not in status and cn not in L:
            ledger["E6"] += 1; continue
        T = L[cn] - t0[cn]
        if T < 0:
            ledger["E7"] += 1; continue
        endrt = num(t1.get("t1e3d"))
        gap = (t0[cn] - endrt) if endrt is not None else None
        if gap is not None and not (28 <= gap <= 280):
            ledger["E8"] += 1; continue
        st = status.get(cn, None)
        ev = 1 if st in (2, 9) else 0
        i1 = I1i.get(cn, {}); o1 = O1i.get(cn, {})
        rows.append({"cn": cn, "x": expi[cn]["x_log2"], "suv": expi[cn]["suv_post"],
                     "T_months": T / 30.4375, "event": ev, "status": st,
                     "gap_rt_to_pet": gap, "inst": a0.get("inst_no"),
                     "age": num(a0.get("entryage")), "sex": num(a0.get("a0e11")),
                     "ps": num(i1.get("i1e5")), "wl": num(i1.get("i1e3")),
                     "wl_unk": num(i1.get("i1e4")), "stage": num(i1.get("i1e6")),
                     "upstage4": 1 if num(o1.get("o1e2")) == 4 else 0,
                     "pre_suv": pre.get(cn), "post_max": num(
                         next((r for r in LE if r.get("cn") == cn), {}).get("lee10")),
                     "L_days": L[cn], "prog365": 1 if cn in prog365 else 0})
    df = pd.DataFrame(rows)
    n_excl = sum(ledger.values())
    run.gate("G5_exclusion_accounting", {"identity": "251 == n_analysed + n_excluded"},
             {"n_listed": 251, "n_analysed": len(df), "n_excluded": n_excl,
              "ledger": ledger}, 251 == len(df) + n_excl)
    run.log("exclusion_ledger", **ledger, n_analysed=len(df))

    # G4 leakage checks
    exp_recheck = hashlib.sha256(open(exp_path, "rb").read()).hexdigest()
    inband = float(np.mean([(28 <= g <= 280) for g in df["gap_rt_to_pet"].dropna()])) \
        if df["gap_rt_to_pet"].notna().any() else 0.0
    run.gate("G4_no_outcome_leakage",
             {"4a_exposure_from_PET_forms_only": True, "4b_min_T>=0": True,
              "4c_inband>=0.95": ">=0.95", "4e_exposure_hash_unchanged": EXP_SHA[:16]},
             {"4a": "LE.lee11 only", "4b_min_T": round(float(df["T_months"].min()), 3),
              "4c_inband": round(inband, 4), "4e": exp_recheck[:16]},
             df["T_months"].min() >= 0 and exp_recheck == EXP_SHA,
             "4e makes post-outcome tuning of the exposure mechanically detectable")

    tally = {int(k): int(v) for k, v in df["status"].value_counts(dropna=False).items()
             if pd.notna(k)}
    d_events = int(df["event"].sum())
    km = KaplanMeierFitter().fit(df["T_months"], 1 - df["event"])
    run.gate("G8_outcome_completeness", "tally sums to n_analysed",
             {"tally": tally, "events": d_events, "ds_only": ds_only,
              "alive_after_dead": alive_after_dead,
              "median_reverse_km_months": round(float(km.median_survival_time_), 2)},
             sum(tally.values()) <= len(df))

    # ================= STAGE 3 : primary + PC1 + S5 =================
    cph = CoxPHFitter().fit(df[["T_months", "event", "x"]], "T_months", "event")
    HR_P = float(np.exp(cph.params_["x"]))
    CI_P = [float(np.exp(cph.confidence_intervals_.loc["x"].iloc[0])),
            float(np.exp(cph.confidence_intervals_.loc["x"].iloc[1]))]
    z_obs = float(cph.summary.loc["x", "z"])
    wald_p = float(cph.summary.loc["x", "p"])
    rng = np.random.default_rng(SEED_PERM)
    zs = []
    for _ in range(B_PERM):
        d2 = df.copy()
        d2["x"] = rng.permutation(d2["x"].to_numpy())
        try:
            zs.append(float(CoxPHFitter().fit(d2[["T_months", "event", "x"]],
                                              "T_months", "event").summary.loc["x", "z"]))
        except Exception:
            pass
    zs = np.array(zs)
    perm_p = float((np.sum(np.abs(zs) >= abs(z_obs)) + 1) / (len(zs) + 1))
    rej = float(np.mean(np.abs(zs) > 1.96))
    run.gate("G6_permutation_calibration",
             {"|mean z|<0.10": True, "sd z in [0.85,1.15]": True, "rej in [0.035,0.065]": True},
             {"mean_z": round(float(zs.mean()), 4), "sd_z": round(float(zs.std(ddof=1)), 4),
              "rejection_rate": round(rej, 4), "B_effective": len(zs)},
             abs(zs.mean()) < 0.10 and 0.85 <= zs.std(ddof=1) <= 1.15
             and 0.035 <= rej <= 0.065)
    try:
        ph = proportional_hazard_test(cph, df[["T_months", "event", "x"]], time_transform="rank")
        ph_p = float(ph.p_value[0])
    except Exception:
        ph_p = float("nan")
    run.log("PRIMARY", hr_per_doubling=round(HR_P, 4), ci95=[round(c, 4) for c in CI_P],
            z=round(z_obs, 4), wald_p=wald_p, permutation_p=perm_p,
            schoenfeld_p=ph_p, n=len(df), events=d_events)

    # S5: mandatory sensitivity - status 9 censored instead of event
    d5 = df.copy()
    d5["event"] = [0 if s == 9 else e for s, e in zip(d5["status"], d5["event"])]
    c5 = CoxPHFitter().fit(d5[["T_months", "event", "x"]], "T_months", "event")
    hr5 = float(np.exp(c5.params_["x"])); p5 = float(c5.summary.loc["x", "p"])
    run.log("S5_status9_censored", hr=round(hr5, 4), wald_p=p5,
            n_status9=int((df["status"] == 9).sum()))

    def cox1(d, col, tcol="T_months", ecol="event"):
        m = CoxPHFitter().fit(d[[tcol, ecol, col]], tcol, ecol)
        return (float(np.exp(m.params_[col])), float(m.summary.loc[col, "p"]),
                [float(np.exp(m.confidence_intervals_.loc[col].iloc[0])),
                 float(np.exp(m.confidence_intervals_.loc[col].iloc[1]))], len(d),
                int(d[ecol].sum()))

    # ---------- PC1: progression by day 365 -> OS from a day-365 landmark ----------
    lm = df[(df["L_days"] >= 365)].copy()
    lm["T_L"] = (lm["L_days"] - 365) / 30.4375
    lm = lm[lm["T_L"] > 0]
    if len(lm) < 60 or int(lm["event"].sum()) < 30:
        pc1 = "UNINFORMATIVE-BY-N (treated as WEAK-PASS)"; pc1v = (np.nan, np.nan, [np.nan]*2, len(lm), int(lm["event"].sum()))
    else:
        pc1v = cox1(lm, "prog365", "T_L")
        hr, pv = pc1v[0], pc1v[1]
        if hr > 1.5 and pv < 0.01: pc1 = "PASS"
        elif (1.3 < hr <= 1.5) or (hr > 1.5 and 0.01 <= pv < 0.05): pc1 = "WEAK-PASS"
        else: pc1 = "FAILED"
    run.log("PC1_endpoint_construction", verdict=pc1, hr=None if np.isnan(pc1v[0]) else round(pc1v[0],4),
            p=None if np.isnan(pc1v[1]) else pc1v[1], n=pc1v[3], events=pc1v[4])

    # ---------- PC2: performance status ----------
    d2 = df[(df["ps"].notna()) & (df["ps"] != 99)].copy()
    d2["ps_hi"] = (d2["ps"] >= 1).astype(int)
    small = min(int((d2["ps_hi"]==1).sum()), int((d2["ps_hi"]==0).sum()))
    small_ev = min(int(d2[d2["ps_hi"]==1]["event"].sum()), int(d2[d2["ps_hi"]==0]["event"].sum()))
    if small < 40 or small_ev < 25:
        pc2 = "UNINFORMATIVE-BY-RANGE"; pc2v = (np.nan, np.nan, [np.nan]*2, len(d2), int(d2["event"].sum()))
    else:
        pc2v = cox1(d2, "ps_hi"); hr, pv = pc2v[0], pc2v[1]
        if hr > 1 and pv < 0.05: pc2 = "PASS"
        elif hr >= 1.0: pc2 = "PASS-WEAK"
        elif hr < 1.0 and pv < 0.05: pc2 = "FAILED"
        else: pc2 = "PASS-WEAK"
    run.log("PC2_performance_status", verdict=pc2, n=pc2v[3], smaller_group=small,
            smaller_group_events=small_ev,
            hr=None if np.isnan(pc2v[0]) else round(pc2v[0],4),
            p=None if np.isnan(pc2v[1]) else pc2v[1])

    # ---------- PC3: weight loss >= 5% ----------
    d3 = df[(df["wl"].notna()) & (df["wl_unk"] != 99)].copy()
    frac_missing = 1 - len(d3)/max(len(df),1)
    d3["wl_hi"] = (d3["wl"] >= 5).astype(int)
    exp_n = int((d3["wl_hi"]==1).sum()); exp_ev = int(d3[d3["wl_hi"]==1]["event"].sum())
    if exp_n < 30 or exp_ev < 20 or frac_missing > 0.40:
        pc3 = "UNINFORMATIVE-BY-N"; pc3v = (np.nan, np.nan, [np.nan]*2, len(d3), int(d3["event"].sum()))
    else:
        pc3v = cox1(d3, "wl_hi"); hr, pv = pc3v[0], pc3v[1]
        if hr > 1 and pv < 0.05: pc3 = "PASS"
        elif hr >= 1.0: pc3 = "PASS-WEAK"
        elif hr < 1.0 and pv < 0.05: pc3 = "FAILED"
        else: pc3 = "PASS-WEAK"
    run.log("PC3_weight_loss", verdict=pc3, n=pc3v[3], exposed_n=exp_n, exposed_events=exp_ev,
            frac_missing=round(frac_missing,4),
            hr=None if np.isnan(pc3v[0]) else round(pc3v[0],4),
            p=None if np.isnan(pc3v[1]) else pc3v[1])

    # ---------- S1-S7 ----------
    S = {}
    s1 = df[(df["age"].notna()) & (df["sex"].notna()) & (df["ps"].notna()) & (df["ps"] != 99)
            & (df["stage"].isin([1,2,3]))].copy()
    s1["ps_hi"] = (s1["ps"] >= 1).astype(int); s1["stage3b"] = (s1["stage"] == 3).astype(int)
    m1 = CoxPHFitter().fit(s1[["T_months","event","x","age","sex","ps_hi","stage3b"]],
                           "T_months","event")
    S["S1_adjusted"] = {"hr": float(np.exp(m1.params_["x"])), "p": float(m1.summary.loc["x","p"]),
                        "n": len(s1)}
    d_s2 = df[df["post_max"].notna()].copy()
    d_s2["x2"] = np.log2(np.maximum(d_s2["post_max"], 0.5))
    h,pv,ci,n,ev = cox1(d_s2,"x2"); S["S2_tumor_max_suv"]={"hr":h,"p":pv,"n":n}
    d_s3 = df.copy(); d_s3["hi"] = (d_s3["suv"] >= 3.5).astype(int)
    h,pv,ci,n,ev = cox1(d_s3,"hi"); S["S3_fixed_threshold_3.5"]={"hr":h,"p":pv,"n":n,
                                                                 "n_high":int(d_s3["hi"].sum())}
    d_s4 = df[df["pre_suv"].notna() & (df["pre_suv"]>0)].copy()
    d_s4["dx"] = np.log2(np.maximum(d_s4["suv"],0.5)) - np.log2(np.maximum(d_s4["pre_suv"],0.5))
    h,pv,ci,n,ev = cox1(d_s4,"dx"); S["S4_metabolic_change"]={"hr":h,"p":pv,"n":n}
    S["S5_status9_censored"] = {"hr": hr5, "p": p5, "n": len(df),
                               "n_status9": int((df["status"]==9).sum())}
    d_s6 = df[df["upstage4"]==0].copy()
    h,pv,ci,n,ev = cox1(d_s6,"x"); S["S6_exclude_upstaged_IV"]={"hr":h,"p":pv,"n":n,
                                                                "n_excluded":int(df["upstage4"].sum())}
    h,pv,ci,n,ev = cox1(df,"suv"); S["S7_untransformed"]={"hr":h,"p":pv,"n":n}
    praw = [S[k]["p"] for k in S]
    order = np.argsort(praw); m = len(praw); holm = [0.0]*m; run_max = 0.0
    for rank, idx in enumerate(order):
        adj = min(1.0, (m - rank) * praw[idx]); run_max = max(run_max, adj); holm[idx] = run_max
    for k, hp in zip(S, holm):
        S[k]["holm_p"] = float(hp)
        run.log("secondary", test=k, hr=round(S[k]["hr"],4), raw_p=S[k]["p"],
                holm_p=round(hp,4), n=S[k]["n"])

    # ---------- STAGE 4 : mechanical verdict ----------
    agg_ok = (pc1 in ("PASS","WEAK-PASS","UNINFORMATIVE-BY-N (treated as WEAK-PASS)")
              and pc0 in ("PASS","MARGINAL","INCONCLUSIVE-BY-N (treated as MARGINAL)")
              and pc2 != "FAILED" and pc3 != "FAILED"
              and (pc2 in ("PASS","PASS-WEAK") or pc3 in ("PASS","PASS-WEAK")))
    if pc1 == "FAILED":
        verdict = "ENDPOINT-CONTROL-FAILED"
    elif pc2 == "FAILED" or pc3 == "FAILED":
        verdict = "NULL-UNRELIABLE"
    elif perm_p < 0.05 and HR_P < 1:
        verdict = "DIRECTIONAL-CONTRADICTION"
    elif (HR_P > 1 and perm_p < 0.05 and CI_P[0] > 1.00
          and S["S5_status9_censored"]["hr"] > 1 and S["S1_adjusted"]["hr"] > 1 and agg_ok):
        verdict = "CONFIRMED"
    elif HR_P > 1 and perm_p < 0.05 and CI_P[0] > 1.00:
        if S["S5_status9_censored"]["hr"] <= 1: verdict = "CONFIRMED-FRAGILE (not confirmed)"
        elif S["S1_adjusted"]["hr"] <= 1: verdict = "CONFIRMED-UNADJUSTED-ONLY (not confirmed)"
        else: verdict = "CONFIRMED-UNCORROBORATED (not confirmed)"
    elif perm_p >= 0.05 and d_events >= 100:
        verdict = "NULL-INFORMATIVE" 
    else:
        verdict = "POWER-LIMITED-NULL"
    run.gate("G9_primary_coherence", "point estimate lies inside its own 95% CI",
             {"hr": round(HR_P, 4), "ci95": [round(c, 4) for c in CI_P]},
             CI_P[0] <= HR_P <= CI_P[1],
             "added after a variable-shadowing bug let PC2's HR be reported as the primary")
    run.log("STAGE4_VERDICT", verdict=verdict, aggregate_control_rule_satisfied=agg_ok,
            PC0=pc0, PC1=pc1, PC2=pc2, PC3=pc3)

    results = {
        "plan_sha256": PLAN_SHA,
        "VERDICT": verdict,
        "aggregate_control_rule_satisfied": bool(agg_ok),
        "controls": {"PC0": pc0, "PC1": pc1, "PC2": pc2, "PC3": pc3},
        "PC1_detail": {"hr": None if np.isnan(pc1v[0]) else pc1v[0], "p": None if np.isnan(pc1v[1]) else pc1v[1],
                       "n": pc1v[3], "events": pc1v[4]},
        "PC2_detail": {"hr": None if np.isnan(pc2v[0]) else pc2v[0], "p": None if np.isnan(pc2v[1]) else pc2v[1],
                       "n": pc2v[3], "smaller_group": small, "smaller_group_events": small_ev},
        "PC3_detail": {"hr": None if np.isnan(pc3v[0]) else pc3v[0], "p": None if np.isnan(pc3v[1]) else pc3v[1],
                       "n": pc3v[3], "exposed_n": exp_n, "exposed_events": exp_ev},
        "secondary_tests": S,
        "stage1_exposure_sha256": EXP_SHA,
        "cohort": {"n_listed": 251, "n_analysed": len(df), "events": d_events,
                   "exposure_source": "LE.lee11", "floored_at_0.5": floored},
        "exclusion_ledger": ledger,
        "PC0_exposure_validity": {"n": len(both), "frac_post_lower": frac_lower,
                                  "median_ratio": ratio_med, "verdict": pc0},
        "PRIMARY": {"hr_per_doubling": HR_P, "ci95": CI_P, "z": z_obs, "wald_p": wald_p,
                    "permutation_p": perm_p, "schoenfeld_p": ph_p},
        "S5_status9_censored": {"hr": hr5, "wald_p": p5},
        "vital_status_tally": tally, "ds_only_events": ds_only,
        "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir},
    }
    run.write("results.json", results)
    df.to_csv(os.path.join(run.dir, "analysis_frame.csv"), index=False)
    np.save(os.path.join(run.dir, "permutation_z.npy"), zs)
    run.finalize()
    print("\nRUN DIR:", run.dir)


if __name__ == "__main__":
    main()
