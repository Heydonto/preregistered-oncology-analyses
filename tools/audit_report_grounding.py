#!/usr/bin/env python3
"""Sweep every report for the defect R19 found: a number with no artefact behind it.

R19 discovered that two of Paper 1's numbers had never been written to a file. The obvious next
question is whether that is isolated or systemic across the sixteen reports, since the reports --
not the manuscripts -- are the IP record.

The test has to be SCOPED to have power. Searching all evidence at once is useless: with 1.5M
distinct numeric strings, every possible 4-decimal AUROC appears somewhere, and the global version
of this check scored 0% detection on fabricated values. So each report is checked against its OWN
`evidence/` directory first, which is small enough that a match means something.

Three tiers, because reports legitimately quote each other:

  OWN        the number appears in this report's own evidence/            -> fine
  CROSS      not in its own evidence, but present in the evidence of a
             report this one explicitly names (e.g. R17 quoting R15)      -> fine, but listed
  UNGROUNDED in neither                                                   -> the R19 defect

Detection power is measured per tier with fabricated values and printed. A tier whose power is low
is reported as uninformative rather than as a clean result -- that was the mistake the first
version of papers/verify_claims.py made.

Run:  python3 results/audit_report_grounding.py
"""
import json
import os
import random
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")

SKIP_CONTEXT = re.compile(
    r"(documentclass|usepackage|geometry|RGB|definecolor|includegraphics|width=|multirow|"
    r"fontsize|titleformat|vspace|hspace|linewidth|textwidth|rule\{|August 2026|July 2026|"
    r"cmidrule|label\{|ref\{|hline|\\\\\[|DOI |doi:|CC BY|CC0|arXiv)", re.I)
TRIVIAL = {str(i) for i in range(0, 13)} | {"2025", "2026", "95", "100", "1000", "05", "0.5"}


def numbers_in_tex(path):
    out = []
    for ln, line in enumerate(open(path, errors="ignore").read().split("\n"), 1):
        if SKIP_CONTEXT.search(line):
            continue
        s = line
        # scientific notation: "$p=1.3\\times10^{-13}$" leaves a bare "13" once \\times is
        # stripped, and the exponent is not a measurement. Remove exponents before harvesting.
        s = re.sub(r"\\times\s*10\s*\^\s*\{?-?\d+\}?", " ", s)
        s = re.sub(r"10\s*\^\s*\{?-?\d+\}?", " ", s)
        # \texttt{...} holds identifiers (slide IDs, file names, hashes), not measurements
        s = re.sub(r"\\texttt\{[^}]*\}", " ", s)
        s = re.sub(r"\\[a-zA-Z]+", " ", s)
        for m in re.finditer(
                r"(?<![\w.])(\d{1,3}(?:(?:\{,\}|,)\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?|\.\d+)(?![\w])",
                s):
            raw = m.group(1)
            c = raw.replace(",", "").replace("{", "").replace("}", "")
            if c in TRIVIAL:
                continue
            try:
                v = float(c)
            except ValueError:
                continue
            out.append((raw, v, ln, line.strip()[:120]))
    return out


def haystack_of(evdir):
    seen = set()
    if not os.path.isdir(evdir):
        return seen
    for dp, _, fs in os.walk(evdir):
        for f in fs:
            p = os.path.join(dp, f)
            try:
                if f.endswith(".json"):
                    def walk(o):
                        if isinstance(o, dict):
                            [walk(v) for v in o.values()]
                        elif isinstance(o, list):
                            [walk(v) for v in o]
                        elif isinstance(o, bool):
                            pass
                        elif isinstance(o, (int, float)):
                            for k in range(0, 7):
                                seen.add(f"{round(float(o), k):.{k}f}")
                            if float(o).is_integer():
                                seen.add(str(int(o)))
                        elif isinstance(o, str):
                            for mm in re.finditer(r"\d*\.?\d+", o):
                                seen.add(mm.group(0))
                    walk(json.load(open(p)))
                elif f.endswith((".csv", ".tsv", ".txt", ".jsonl", ".yaml", ".md", ".py")):
                    if os.path.getsize(p) > 8_000_000:
                        continue
                    for mm in re.finditer(r"\d*\.?\d+", open(p, errors="ignore").read()):
                        seen.add(mm.group(0))
                elif f.endswith(".npy"):
                    a = np.load(p, mmap_mode="r")
                    if a.size and a.dtype.kind in "fiu":
                        seen.add(str(int(a.size)))
                        if a.size < 3_000_000:
                            v = np.asarray(a).ravel()
                            for x in (v.min(), v.max(), np.median(v), v.mean()):
                                for k in range(0, 5):
                                    seen.add(f"{round(float(x), k):.{k}f}")
            except Exception:
                continue
    seen.discard("")
    return seen


def hit(v, hay, decimals=None):
    """Match at the precision the number is WRITTEN with, not at every precision.

    The first version generated candidates at 7 rounding precisions plus x100 and /100 variants
    -- about 28 strings per value -- so low-precision candidates like "0" or "1" matched every
    haystack and detection power sat at 33%. Matching at the stated precision (plus one either
    side for rounding, plus a percent variant at the corresponding precision) takes power to a
    level where a clean result means something.
    """
    if decimals is None:
        decimals = 0 if float(v).is_integer() else 6
    c = set()
    for k in {max(decimals - 1, 0), decimals, decimals + 1}:
        c.add(f"{round(v, k):.{k}f}")
    if float(v).is_integer():
        c.add(str(int(v)))
    # a fraction in evidence quoted as a percentage in prose, or the reverse, at like precision
    for alt, shift in ((v / 100.0, decimals + 2), (v * 100.0, max(decimals - 2, 0))):
        c.add(f"{round(alt, shift):.{shift}f}")
        if float(alt).is_integer():
            c.add(str(int(alt)))
    c |= {("-" + x) for x in list(c)}
    return bool(c & hay)


def decimals_of(raw):
    r = raw.replace(",", "").replace("{", "").replace("}", "")
    return len(r.split(".")[1]) if "." in r else 0


def power(hay, n=1500, seed=7):
    """Detection power against fabricated values of the kinds these reports quote."""
    rng = random.Random(seed)
    vals = ([(round(rng.uniform(0, 1), 4), 4) for _ in range(n // 3)]
            + [(round(rng.uniform(1, 3), 3), 3) for _ in range(n // 3)]
            + [(rng.randint(1000, 999999), 0) for _ in range(n // 3)])
    fp = sum(1 for v, d in vals if hit(float(v), hay, d)) / len(vals)
    return 100.0 * (1 - fp)


def main():
    reports = sorted(d for d in os.listdir(RES)
                     if os.path.isdir(os.path.join(RES, d))
                     and os.path.exists(os.path.join(RES, d, "REPORT.tex")))
    hays = {r: haystack_of(os.path.join(RES, r, "evidence")) for r in reports}
    rows, flagged = [], []
    for r in reports:
        tex = os.path.join(RES, r, "REPORT.tex")
        nums = numbers_in_tex(tex)
        own = hays[r]
        named = [o for o in reports if o != r
                 and re.search(r"\b" + re.escape(o.split("-")[0]) + r"\b",
                               open(tex, errors="ignore").read())]
        cross = set().union(*[hays[o] for o in named]) if named else set()
        o_ct = c_ct = u_ct = 0
        ung = {}
        for raw, v, ln, line in nums:
            dec = decimals_of(raw)
            if hit(v, own, dec):
                o_ct += 1
            elif hit(v, cross, dec):
                c_ct += 1
            else:
                u_ct += 1
                ung.setdefault(raw, (ln, line))
        rows.append((r, len(nums), o_ct, c_ct, u_ct, power(own), len(own), named))
        if ung:
            flagged.append((r, ung))

    print(f"{'report':38s} {'nums':>5s} {'own':>5s} {'cross':>6s} {'UNGR':>5s} "
          f"{'power%':>7s} {'|hay|':>7s}")
    for r, n, o, c, u, pw, hs, named in rows:
        mark = "  <-- CHECK" if u else ""
        print(f"{r:38s} {n:5d} {o:5d} {c:6d} {u:5d} {pw:7.1f} {hs:7d}{mark}")

    if flagged:
        print("\nnumbers not found in their own or any cross-referenced report's evidence:")
        for r, ung in flagged:
            print(f"\n  {r}")
            for raw, (ln, line) in sorted(ung.items()):
                print(f"    {raw:>12s}  L{ln}: {line}")
    else:
        print("\nno ungrounded numbers found in any report.")

    weak = [(r, pw) for r, n, o, c, u, pw, hs, nm in rows if pw < 60]
    if weak:
        print("\nreports where this test is UNINFORMATIVE (power < 60%), listed so the clean "
              "result above is not over-read:")
        for r, pw in weak:
            print(f"    {r:38s} power {pw:.1f}%")
    tot_u = sum(r[4] for r in rows)
    print(f"\ntotal ungrounded numbers across {len(rows)} reports: {tot_u}")
    return 1 if tot_u else 0


if __name__ == "__main__":
    sys.exit(main())
