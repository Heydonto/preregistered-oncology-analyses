#!/usr/bin/env python3
"""R20 - the external site probe across both encoders, on the full slide set.

TWO GAPS THIS CLOSES.

1. R18 established that subtype site-leakage is a property of the archive (it reproduces, larger,
   under an encoder trained on natural images) while methylation leakage was a property of
   Phikon-v2 features. That dissociation makes the *external* probe -- EAGLE vs HTMCP, which R19
   measured at balanced accuracy 1.000 -- an open question rather than a settled one. R18's audit
   lists it as un-run.

2. R19's probe used 78 slides: 49 EAGLE and 29 HTMCP. The archive holds **80** HTMCP slides. The
   78-slide file was built from a partial local staging, so **51 HTMCP slides were silently
   omitted**, and neither R19 nor Paper 1 says so. That is not a wrong number, but a
   non-representative subset described as "two external cohorts" is a defect of the same family as
   the H&E mislabelling R19 already found in the same file.

So this runs both encoders over three nested slide sets, with identical mean-pooled features,
identical probe and identical permutation scheme:

  SET_R19    49 EAGLE + 29 HTMCP = 78   the set R19 and Paper 1 actually used
  SET_HE     49 EAGLE + 27 HTMCP = 76   every H&E slide in the archive
  SET_ALL    49 EAGLE + 80 HTMCP = 129  everything, H&E and IHC together

SET_ALL is expected to separate trivially and is included as an upper bound rather than a result:
with IHC and H&E mixed, a probe can win on stain alone. SET_HE is the honest measurement.

Slide-id sets are identical between encoders (verified), so any difference is the encoder.

Run:  python3 m17_external_probe_two_encoders.py
"""
import collections
import json
import os

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from audit_kit import Run

PH = "/tmp/ph_ext_mean.npz"
DV = "/tmp/dv_ext_mean.npz"
R19_SET = "/Users/rezanehzati/quantara-staging/ext_feat/external_he_mean.npz"

SEED = 20260817
B_PERM = 200
C_REG = 1.0
N_SPLITS = 5

R19_PUBLISHED = 1.0000
CONFIRM = 0.95
WEAK = 0.75


def probe(X, y, seed):
    y = np.asarray(y)
    k = min(N_SPLITS, min(collections.Counter(y).values()))
    if k < 2:
        return float("nan")
    oof = np.empty(len(y), dtype=y.dtype)
    for tr, te in StratifiedKFold(n_splits=k, shuffle=True, random_state=seed).split(X, y):
        m = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=5000, C=C_REG)).fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    return float(balanced_accuracy_score(y, oof))


def with_null(X, y, seed):
    obs = probe(X, y, seed)
    rng = np.random.default_rng(seed)
    null = np.array([probe(X, rng.permutation(y), seed) for _ in range(B_PERM)])
    return {"balanced_accuracy": obs, "n": int(len(y)),
            "classes": {str(k): int(v) for k, v in sorted(collections.Counter(y).items())},
            "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)),
            "null_p95": float(np.quantile(null, 0.95)),
            "permutation_p": float((1 + int((null >= obs).sum())) / (B_PERM + 1))}


def main():
    run = Run("R20extprobe2enc")
    cfg_hash = run.start(
        {"question": "does the external EAGLE-vs-HTMCP separation reproduce under a second "
                     "encoder, and what happens on the 51 HTMCP slides R19 omitted?",
         "why": "R18 found subtype leakage archival but methylation leakage encoder-specific, so "
                "the external probe is an open question. Separately, R19's 78-slide input silently "
                "omitted 51 of the archive's 80 HTMCP slides.",
         "sets": {"SET_R19": "49 EAGLE + 29 HTMCP = 78, the set R19 and Paper 1 used",
                  "SET_HE": "49 EAGLE + 27 HTMCP = 76, every H&E slide in the archive",
                  "SET_ALL": "49 EAGLE + 80 HTMCP = 129, upper bound only -- mixing IHC and H&E "
                             "lets a probe win on stain alone, so this is not a result"},
         "encoders": ["owkin/phikon-v2", "facebook/dinov2-large"],
         "probe": {"model": "StandardScaler + LogisticRegression", "C": C_REG,
                   "cv": f"StratifiedKFold({N_SPLITS}, shuffled)",
                   "metric": "out-of-fold balanced accuracy"},
         "permutations": B_PERM, "p_formula": "(1 + #{null >= obs}) / (B + 1)", "seed": SEED,
         "primary_arm": "SET_HE under facebook/dinov2-large",
         "verdict_rules": [
             f"dinov2 SET_HE >= {CONFIRM} -> EXTERNAL_SIGNATURE_ENCODER_INDEPENDENT",
             f"dinov2 SET_HE >= {WEAK} -> EXTERNAL_SIGNATURE_WEAKER_UNDER_SECOND_ENCODER",
             f"dinov2 SET_HE < {WEAK} -> EXTERNAL_SIGNATURE_PHIKON_SPECIFIC"],
         "status": "HELD FOR IP - not for publication pending patent"},
        [PH, DV, R19_SET])

    E = {}
    for nm, p in (("phikon-v2", PH), ("dinov2-large", DV)):
        d = np.load(p, allow_pickle=True)
        E[nm] = {"sid": [str(x) for x in d["slide_ids"]], "X": d["mean_emb"].astype(np.float64),
                 "coh": np.array([str(x) for x in d["cohort"]], object),
                 "he": np.asarray(d["is_he"]), "enc": str(d["encoder"])}

    a, b = E["phikon-v2"], E["dinov2-large"]
    run.gate("G0_matched_slides",
             "both encoders cover the same 129 slides in the same order",
             {"phikon_n": len(a["sid"]), "dinov2_n": len(b["sid"]),
              "identical": a["sid"] == b["sid"]},
             a["sid"] == b["sid"] and len(a["sid"]) == 129,
             "any difference below is then attributable to the encoder alone")
    run.gate("G1_encoder_identity", "each file declares the intended encoder",
             {"phikon": a["enc"], "dinov2": b["enc"]},
             a["enc"] == "owkin/phikon-v2" and b["enc"] == "facebook/dinov2-large")

    r19 = {str(x) for x in np.load(R19_SET, allow_pickle=True)["slide_ids"]}
    in_r19 = np.array([s in r19 for s in a["sid"]])
    omitted = int(sum(1 for s, c in zip(a["sid"], a["coh"]) if c == "htmcp" and s not in r19))
    run.gate("G2_r19_subset_confirmed",
             "R19's 78-slide set is a strict subset, and the number of HTMCP slides it omitted",
             {"r19_n": int(in_r19.sum()), "archive_htmcp": int((a["coh"] == "htmcp").sum()),
              "htmcp_omitted_by_r19": omitted},
             int(in_r19.sum()) == 78 and omitted == 51,
             "asserts the defect is present and quantified; if this changes, the premise of this "
             "run needs re-checking")

    SETS = {"SET_R19": in_r19, "SET_HE": a["he"], "SET_ALL": np.ones(len(a["sid"]), bool)}
    res = {}
    for enc in ("phikon-v2", "dinov2-large"):
        d = E[enc]
        res[enc] = {}
        for sn, mask in SETS.items():
            res[enc][sn] = with_null(d["X"][mask], d["coh"][mask], SEED)
            r = res[enc][sn]
            run.log("probe", encoder=enc, set=sn, n=r["n"],
                    bal_acc=round(r["balanced_accuracy"], 4),
                    null_mean=round(r["null_mean"], 4), p=round(r["permutation_p"], 5))

    run.gate("G3_reproduces_r19",
             f"phikon on SET_R19 recovers R19's {R19_PUBLISHED}",
             round(res["phikon-v2"]["SET_R19"]["balanced_accuracy"], 4),
             abs(res["phikon-v2"]["SET_R19"]["balanced_accuracy"] - R19_PUBLISHED) <= 0.02)
    run.gate("G4_nulls_centred", "every null sits at chance 0.5 (within 0.06)",
             {f"{e}/{s}": round(res[e][s]["null_mean"], 4) for e in res for s in SETS},
             all(abs(res[e][s]["null_mean"] - 0.5) < 0.06 for e in res for s in SETS))

    prim = res["dinov2-large"]["SET_HE"]["balanced_accuracy"]
    if prim >= CONFIRM:
        verdict = "EXTERNAL_SIGNATURE_ENCODER_INDEPENDENT"
    elif prim >= WEAK:
        verdict = "EXTERNAL_SIGNATURE_WEAKER_UNDER_SECOND_ENCODER"
    else:
        verdict = "EXTERNAL_SIGNATURE_PHIKON_SPECIFIC"
    run.log("VERDICT", primary=verdict, dinov2_set_he=round(prim, 4),
            phikon_set_he=round(res["phikon-v2"]["SET_HE"]["balanced_accuracy"], 4))

    out = {"status": "HELD FOR IP - not for publication pending patent",
           "VERDICT": {"primary": verdict, "primary_arm": "dinov2-large / SET_HE",
                       "value": prim, "rules_pre_declared_in": "config.yaml, hash " + cfg_hash},
           "arms": res,
           "sets": {sn: {"n": int(m.sum()),
                         "cohorts": {c: int(((a["coh"] == c) & m).sum())
                                     for c in ("eagle", "htmcp")}}
                    for sn, m in SETS.items()},
           "r19_subset_defect": {
               "r19_n": int(in_r19.sum()), "archive_htmcp": int((a["coh"] == "htmcp").sum()),
               "htmcp_omitted_by_r19": omitted,
               "note": "R19's input was built from a partial local staging; neither R19 nor "
                       "Paper 1 records that 51 of 80 HTMCP slides were absent"},
           "set_all_caveat": "SET_ALL mixes IHC and H&E; a probe can separate on stain alone, so "
                             "it is an upper bound and not a measurement of site signal",
           "_provenance": {"config_sha256": cfg_hash, "run_dir": run.dir}}
    run.write("results.json", out)
    run.finalize()
    print("\nVERDICT:", verdict)
    for enc in res:
        for sn in SETS:
            r = res[enc][sn]
            print(f"  {enc:16s} {sn:8s} n={r['n']:3d}  bal.acc {r['balanced_accuracy']:.4f}  "
                  f"null {r['null_mean']:.4f}  p={r['permutation_p']:.5f}")
    print("RUN DIR:", run.dir)


if __name__ == "__main__":
    main()
