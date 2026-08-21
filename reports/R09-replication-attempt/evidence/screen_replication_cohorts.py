#!/usr/bin/env python3
"""Systematic screen: is there ANY public cohort that can replicate R06?

R09 (TCGA) failed its own positive controls. R10's first candidate failed a setting gate
(median OS 93.8 months in a nominally metastatic cohort). Rather than stop after two tries,
this screens every lung study in cBioPortal against the four requirements a replication cohort
must meet, and reports the screen as a table so the search is auditable rather than anecdotal.

Requirements:
  1. overall-survival time and status present
  2. mutation profiling present (for KEAP1 / SMARCA4 / STK11 / TP53)
  3. enough events - R06's failed attempt had 45; require >= 100
  4. advanced-disease setting - median OS nearer R06's 28.8 months than TCGA's 50.3
  5. patient-disjoint from R06's cohort

Writes screen_results.json. No outcome (hazard ratio) is computed here; this only decides
eligibility, so the screen cannot be steered by results.
"""
import json, subprocess, sys
import numpy as np

R06_STUDY = "nsclc_ctdx_msk_2022"
MIN_EVENTS = 100
GENES = [9817, 6597, 8085, 6794, 7157]


def curl(url, body=None, timeout=180):
    cmd = ["curl", "-s", "--max-time", str(timeout), "--retry", "2", url]
    if body is not None:
        cmd += ["-X", "POST", "-H", "Content-Type: application/json", "-d", body]
    try:
        return json.loads(subprocess.run(cmd, capture_output=True, timeout=timeout + 30).stdout)
    except Exception:
        return None


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def km_median(T, E):
    """Kaplan-Meier median without pulling in lifelines for a screen."""
    order = np.argsort(T)
    T, E = np.asarray(T)[order], np.asarray(E)[order]
    n, s = len(T), 1.0
    for i, (t, e) in enumerate(zip(T, E)):
        if e:
            s *= 1 - 1 / (n - i)
            if s <= 0.5:
                return float(t)
    return None


def main():
    studies = curl("https://www.cbioportal.org/api/studies")
    lung = [s for s in studies
            if any(k in (s.get("name", "") + s.get("studyId", "")).lower()
                   for k in ("lung", "nsclc", "luad", "lusc", "sclc"))]
    print(f"screening {len(lung)} lung studies\n")

    r06 = curl(f"https://www.cbioportal.org/api/studies/{R06_STUDY}/patients?pageSize=100000")
    r06_pat = {x["patientId"] for x in r06} if r06 else set()

    rows = []
    for s in sorted(lung, key=lambda x: -(x.get("allSampleCount") or 0)):
        sid = s["studyId"]
        rec = {"study": sid, "name": s.get("name", "")[:60],
               "samples": s.get("allSampleCount")}
        cl = curl(f"https://www.cbioportal.org/api/studies/{sid}/clinical-data"
                  f"?clinicalDataType=PATIENT&projection=SUMMARY&pageSize=200000")
        if not cl:
            rec["verdict"] = "no clinical data returned"; rows.append(rec); continue
        om = {x["patientId"]: x["value"] for x in cl if x["clinicalAttributeId"] == "OS_MONTHS"}
        st = {x["patientId"]: x["value"] for x in cl if x["clinicalAttributeId"] == "OS_STATUS"}
        rec["has_os"] = bool(om and st)
        if not rec["has_os"]:
            rec["verdict"] = "no OS endpoint"; rows.append(rec)
            print(f"  {sid:38s} no OS endpoint"); continue

        usable = [p for p in om if st.get(p) and (num(om[p]) or 0) > 0]
        disj = [p for p in usable if p not in r06_pat]
        rec["n_os"] = len(usable)
        rec["n_disjoint"] = len(disj)
        rec["overlap_with_r06"] = len(usable) - len(disj)
        ev = sum(1 for p in disj if st[p].startswith("1"))
        rec["events_disjoint"] = ev
        med = km_median([num(om[p]) for p in disj], [st[p].startswith("1") for p in disj]) \
            if disj else None
        rec["median_os_months"] = round(med, 1) if med else None

        prof = curl(f"https://www.cbioportal.org/api/studies/{sid}/molecular-profiles")
        rec["has_mutations"] = bool(prof and any(
            p.get("molecularAlterationType") == "MUTATION_EXTENDED" for p in (prof or [])))

        fails = []
        if not rec["has_mutations"]:
            fails.append("no mutation profile")
        if ev < MIN_EVENTS:
            fails.append(f"only {ev} events")
        if med is None:
            fails.append("median OS not reached")
        elif abs(med - 28.8) >= abs(med - 50.3):
            fails.append(f"median OS {med:.1f}mo - not advanced-disease setting")
        if sid == R06_STUDY:
            fails.append("is R06's own cohort")
        if rec["n_disjoint"] == 0:
            fails.append("no disjoint patients")
        rec["fails"] = fails
        rec["eligible"] = not fails
        rec["verdict"] = "ELIGIBLE" if not fails else "; ".join(fails)
        rows.append(rec)
        print(f"  {sid:38s} n_os={rec['n_os']:>4} ev={ev:>4} "
              f"medOS={rec['median_os_months']} -> {rec['verdict'][:56]}")

    elig = [r for r in rows if r.get("eligible")]
    print(f"\n{'=' * 70}\nELIGIBLE COHORTS: {len(elig)}")
    for r in elig:
        print(f"  {r['study']}  n={r['n_disjoint']} events={r['events_disjoint']} "
              f"medOS={r['median_os_months']}")
    out = {"screened": len(rows), "min_events": MIN_EVENTS, "r06_study": R06_STUDY,
           "eligible": [r["study"] for r in elig], "rows": rows}
    json.dump(out, open("screen_results.json", "w"), indent=2, default=str)
    print("\nwrote screen_results.json")


if __name__ == "__main__":
    main()
