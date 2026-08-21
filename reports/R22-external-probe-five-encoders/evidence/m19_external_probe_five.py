#!/usr/bin/env python3
"""R22 - the external site probe across all five encoders.

R20 measured EAGLE vs HTMCP on two encoders and found the signature encoder-independent in
existence but not in magnitude: Phikon-v2 1.0000, dinov2-large 0.8768. R21 then showed that the
TCGA methylation dissociation splits by training corpus rather than by vendor, and that subtype
leakage, while universal, varies more than fourfold in relative magnitude across encoders.

That makes R20's two-encoder external result the last two-point claim in the series, and R21 is a
standing demonstration of what two points are worth here. This runs the same probe on all five.

Design is R20's, unchanged: mean-pooled slide features, standardise, logistic regression C=1.0,
5-fold stratified out-of-fold balanced accuracy, 200 label permutations,
p = (1 + #{null >= obs}) / (B + 1). Three nested slide sets:

  SET_R19    49 EAGLE + 29 HTMCP = 78   the set R19 and Paper 1 originally used
  SET_HE     49 EAGLE + 27 HTMCP = 76   every H&E slide in the archive  <- primary
  SET_ALL    49 EAGLE + 80 HTMCP = 129  upper bound only; mixing IHC and H&E lets a probe win
                                        on stain, and R20 showed it does

PRE-DECLARED, reusing R20's thresholds verbatim (hashed before R20's arms were read, so applying
them to new encoders is same-rule-new-data):

  every encoder's SET_HE >= 0.95            -> EXTERNAL_SIGNATURE_TOTAL_IN_ALL_ENCODERS
  every encoder's SET_HE clears its null    -> EXTERNAL_SIGNATURE_PRESENT_MAGNITUDE_VARIES
  any encoder fails to clear its null       -> EXTERNAL_SIGNATURE_NOT_UNIVERSAL

Feature width varies (1024/1536/1280); the probe standardises per feature so this is handled, but
it is recorded because a wider representation gives a linear probe more room and that is a real
confound on absolute magnitude.

Run:  python3 m19_external_probe_five.py
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

ENC = [("owkin/phikon-v2", "/tmp/ext_phikon-v2.npz", "histology"),
       ("facebook/dinov2-large", "/tmp/ext_dinov2-large.npz", "natural images"),
       ("MahmoodLab/UNI", "/tmp/ext_uni.npz", "histology"),
       ("bioptimus/H-optimus-0", "/tmp/ext_h-optimus-0.npz", "histology"),
       ("paige-ai/Virchow2", "/tmp/ext_virchow2.npz", "histology")]
R19_SET = "/Users/rezanehzati/quantara-staging/ext_feat/external_he_mean.npz"

# R20's seed, NOT a new one. The first attempt at this run used a fresh seed and the reproduction
# gate halted: dinov2-large came out 0.9240 against R20's 0.8768. The seed sets the StratifiedKFold
# split, and on 76 samples with ~15 per fold a couple of reassignments move balanced accuracy by
# ~0.05. So R20's single-seed 0.8768 was never a stable quantity, and neither would a new one be.
# This run reproduces R20 exactly on its own seed AND reports the seed distribution (N_SEEDS below),
# which is the number that should have been quoted in the first place.
SEED = 20260817
N_SEEDS = 25
B_PERM = 200
C_REG = 1.0
N_SPLITS = 5
TOTAL_BAR = 0.95           # R20's threshold, reused
R20_REF = {"owkin/phikon-v2": 1.0000, "facebook/dinov2-large": 0.8768}
TOL = 0.0001


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


def seed_distribution(X, y, seeds):
    """Balanced accuracy across many CV splits. One split on 76 samples is a noisy estimate; this
    is what R20 should have reported instead of a single number."""
    v = np.array([probe(X, y, s) for s in seeds])
    return {"mean": float(v.mean()), "sd": float(v.std(ddof=1)),
            "min": float(v.min()), "max": float(v.max()),
            "n_seeds": len(seeds), "values": [round(float(x), 4) for x in v]}


def with_null(X, y, seed):
    obs = probe(X, y, seed)
    rng = np.random.default_rng(seed)
    null = np.array([probe(X, rng.permutation(y), seed) for _ in range(B_PERM)])
    return {"balanced_accuracy": obs, "n": int(len(y)),
            "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)),
            "null_p95": float(np.quantile(null, 0.95)),
            "clears_null": bool(obs > np.quantile(null, 0.95)),
            "permutation_p": float((1 + int((null >= obs).sum())) / (B_PERM + 1))}


def main():
    run = Run("R22extfive")
    cfg = run.start(
        {"question": "does the external EAGLE-vs-HTMCP separation hold across all five encoders?",
         "why": "R20 measured it on two. R21 showed two-point claims in this series have been "
                "wrong twice; this is the last two-point claim left.",
         "encoders": {r: c for r, _, c in ENC},
         "sets": {"SET_R19": "78, the set R19 and Paper 1 used",
                  "SET_HE": "76, every H&E slide in the archive (PRIMARY)",
                  "SET_ALL": "129, upper bound only -- mixes IHC with H&E"},
         "probe": {"model": "StandardScaler + LogisticRegression", "C": C_REG,
                   "cv": f"StratifiedKFold({N_SPLITS}, shuffled)",
                   "metric": "out-of-fold balanced accuracy"},
         "permutations": B_PERM, "seed": SEED,
         "reused_pre_declared_threshold": {"total_bar": TOTAL_BAR,
                                          "source": "R20 config ab01b9847dedcfcf..., fixed before "
                                                    "R20's arms were read"},
         "known_confound": "feature width varies 1024/1536/1280; a wider representation gives a "
                           "linear probe more room, which bears on absolute magnitude",
         "verdict_rules": [
             f"all SET_HE >= {TOTAL_BAR} -> EXTERNAL_SIGNATURE_TOTAL_IN_ALL_ENCODERS",
             "all SET_HE clear their null p95 -> EXTERNAL_SIGNATURE_PRESENT_MAGNITUDE_VARIES",
             "any fails to clear its null -> EXTERNAL_SIGNATURE_NOT_UNIVERSAL"],
         "status": "peer-review audit copy"},
        [p for _, p, _ in ENC] + [R19_SET])

    D = {}
    for repo, path, corpus in ENC:
        d = np.load(path, allow_pickle=True)
        D[repo] = {"sid": [str(x) for x in d["slide_ids"]], "X": d["mean_emb"].astype(np.float64),
                   "coh": np.array([str(x) for x in d["cohort"]], object),
                   "he": np.asarray(d["is_he"]), "corpus": corpus,
                   "dim": int(d["mean_emb"].shape[1]), "declared": str(d["encoder"])}

    ref = D["owkin/phikon-v2"]["sid"]
    run.gate("G0_matched_slides", "all five cover the same 129 slides in the same order",
             {r: (v["sid"] == ref) for r, v in D.items()},
             all(v["sid"] == ref for v in D.values()) and len(ref) == 129)
    run.gate("G1_encoder_identity", "each matrix declares the encoder it should",
             {r: v["declared"] for r, v in D.items()},
             all(r == v["declared"] for r, v in D.items()))

    r19 = {str(x) for x in np.load(R19_SET, allow_pickle=True)["slide_ids"]}
    in_r19 = np.array([s in r19 for s in ref])
    SETS = {"SET_R19": in_r19, "SET_HE": D["owkin/phikon-v2"]["he"],
            "SET_ALL": np.ones(len(ref), bool)}
    run.gate("G2_set_sizes", "78 / 76 / 129",
             {k: int(v.sum()) for k, v in SETS.items()},
             [int(SETS[k].sum()) for k in ("SET_R19", "SET_HE", "SET_ALL")] == [78, 76, 129])

    res = {}
    for repo, v in D.items():
        res[repo] = {"corpus": v["corpus"], "dim": v["dim"], "arms": {}}
        for sn, mask in SETS.items():
            res[repo]["arms"][sn] = with_null(v["X"][mask], v["coh"][mask], SEED)
            a = res[repo]["arms"][sn]
            run.log("probe", encoder=repo, set=sn, n=a["n"],
                    bal_acc=round(a["balanced_accuracy"], 4),
                    null_mean=round(a["null_mean"], 4), clears=a["clears_null"],
                    p=round(a["permutation_p"], 5))

    # seed distribution on the primary set -- the honest magnitude
    seeds = list(range(SEED, SEED + N_SEEDS))
    for repo, v in D.items():
        m = SETS["SET_HE"]
        res[repo]["set_he_seed_distribution"] = seed_distribution(v["X"][m], v["coh"][m], seeds)
        sd = res[repo]["set_he_seed_distribution"]
        run.log("seed_dist", encoder=repo, mean=round(sd["mean"], 4), sd=round(sd["sd"], 4),
                lo=round(sd["min"], 4), hi=round(sd["max"], 4), n_seeds=sd["n_seeds"])

    # R20 reproduction check on the two encoders it measured
    rep = {r: round(res[r]["arms"]["SET_HE"]["balanced_accuracy"], 4) for r in R20_REF}
    run.gate("G3_reproduces_r20", f"R20's SET_HE values recovered: {R20_REF}", rep,
             all(abs(res[r]["arms"]["SET_HE"]["balanced_accuracy"] - v) <= TOL
                 for r, v in R20_REF.items()),
             "the probe is deterministic given the seed, so this should be exact")
    run.gate("G4_nulls_centred", "every null sits at chance 0.5 (within 0.06)",
             {f"{r}/{s}": round(res[r]["arms"][s]["null_mean"], 4) for r in res for s in SETS},
             all(abs(res[r]["arms"][s]["null_mean"] - 0.5) < 0.06 for r in res for s in SETS))

    he = {r: res[r]["arms"]["SET_HE"] for r in res}
    if all(v["balanced_accuracy"] >= TOTAL_BAR for v in he.values()):
        verdict = "EXTERNAL_SIGNATURE_TOTAL_IN_ALL_ENCODERS"
    elif all(v["clears_null"] for v in he.values()):
        verdict = "EXTERNAL_SIGNATURE_PRESENT_MAGNITUDE_VARIES"
    else:
        verdict = "EXTERNAL_SIGNATURE_NOT_UNIVERSAL"
    vals = sorted(v["balanced_accuracy"] for v in he.values())
    run.log("VERDICT", primary=verdict, he_range=[round(vals[0], 4), round(vals[-1], 4)])

    out = {"status": "peer-review audit copy",
           "VERDICT": {"primary": verdict, "primary_arm": "SET_HE (76 H&E slides)",
                       "he_balanced_accuracy_range": [vals[0], vals[-1]],
                       "n_encoders": len(res),
                       "rules_reused_from": "R20 config ab01b9847dedcfcf..., pre-declared"},
           "encoders": res,
           "by_corpus": {c: sorted(round(res[r]["arms"]["SET_HE"]["balanced_accuracy"], 4)
                                   for r in res if res[r]["corpus"] == c)
                         for c in ("histology", "natural images")},
           "r20_reproduction": rep,
           "seed_sensitivity_note": (
               "R20 quoted single-seed balanced accuracies. A first attempt at this run with a "
               "different seed produced 0.9240 for dinov2-large against R20's 0.8768 and the "
               "reproduction gate halted. On 76 samples with ~15 per fold, one split is a noisy "
               "estimate. set_he_seed_distribution over %d splits is the quantity to quote; the "
               "single-seed values are retained only to demonstrate reproduction." % N_SEEDS),
           "set_all_caveat": "SET_ALL mixes IHC with H&E; higher numbers there indicate a stain "
                             "shortcut, not more site signal",
           "_provenance": {"config_sha256": cfg, "run_dir": run.dir}}
    run.write("results.json", out)
    run.finalize()
    print("\nVERDICT:", verdict)
    print("  SET_HE across %d CV splits (mean +- sd, range):" % N_SEEDS)
    for r in res:
        d = res[r]["set_he_seed_distribution"]
        print(f"    {r:24s} {d['mean']:.4f} +- {d['sd']:.4f}   [{d['min']:.4f}, {d['max']:.4f}]")
    for r in res:
        a = res[r]["arms"]
        print(f"  {r:24s} d={res[r]['dim']:5d}  R19 {a['SET_R19']['balanced_accuracy']:.4f}  "
              f"HE {a['SET_HE']['balanced_accuracy']:.4f}  ALL {a['SET_ALL']['balanced_accuracy']:.4f}"
              f"  (null {a['SET_HE']['null_mean']:.4f})")
    print("RUN DIR:", run.dir)


if __name__ == "__main__":
    main()
