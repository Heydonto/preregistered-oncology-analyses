#!/usr/bin/env python3
"""R08 acquisition: TCGA-LUAD/LUSC 450k methylation -> one compact beta matrix.

Downloads the 919 open-access SeSAMe level-3 beta files named in MANIFEST_450k.tsv, verifies
each against the MD5 GDC published for it, and assembles them into a single float32 matrix.
Raw files are deleted immediately after parsing, so scratch never holds more than the
concurrent working set.

Only the matrix comes back off this VM (~1.8 GB), not the 11.23 GiB of raw input. GDC is a
permanent archive keyed by file_id, so the manifest plus MD5s is sufficient provenance to
re-acquire the raw data; storing it would be redundant.

Integrity rules, all fatal:
  * every download must match its published MD5
  * every file must carry the identical probe vector, verified by hash against the first file
  * fewer than 919 successfully parsed files is a failure, not a smaller matrix
"""
import concurrent.futures as cf, csv, hashlib, os, subprocess, sys, time
import numpy as np

BUCKET = "gs://heydonto-quantara-lungcdx"
RUN = f"{BUCKET}/data-request-2026-08/_run/r08"
OUT = f"{BUCKET}/data-request-2026-08/lung/tcga_methylation_450k"
SCRATCH = "/mnt/scratch/r08"
WORKERS = int(os.environ.get("WORKERS", "12"))
EXPECT_FILES = 919
EXPECT_PROBES = 486427


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def fetch(rec):
    """Download one file, verify MD5, return (index, probes_hash, beta_vector)."""
    fid, dst = rec["file_id"], os.path.join(SCRATCH, rec["file_id"])
    for attempt in range(4):
        rc = subprocess.run(["curl", "-sS", "--max-time", "300", "--retry", "2",
                             "-o", dst, f"https://api.gdc.cancer.gov/data/{fid}"]).returncode
        if rc == 0 and os.path.exists(dst):
            h = hashlib.md5()
            with open(dst, "rb") as fh:
                for b in iter(lambda: fh.read(1 << 20), b""):
                    h.update(b)
            if h.hexdigest() == rec["md5sum"]:
                break
            log(f"  MD5 MISMATCH {fid} attempt {attempt+1}")
        time.sleep(2 * (attempt + 1))
    else:
        os.path.exists(dst) and os.remove(dst)
        return rec["_i"], None, None

    probes, beta = [], np.empty(EXPECT_PROBES, np.float32)
    n = 0
    with open(dst) as fh:
        for line in fh:
            p, _, v = line.rstrip("\n").partition("\t")
            probes.append(p)
            beta[n] = np.nan if v in ("NA", "", "NaN") else float(v)
            n += 1
    os.remove(dst)
    if n != EXPECT_PROBES:
        return rec["_i"], f"BADROWS:{n}", None
    ph = hashlib.sha256("\n".join(probes).encode()).hexdigest()
    return rec["_i"], ph, (beta, probes)


def main():
    os.makedirs(SCRATCH, exist_ok=True)
    subprocess.run(["gcloud", "storage", "cp", f"{RUN}/MANIFEST_450k.tsv",
                    "/tmp/MANIFEST_450k.tsv"], check=True)
    recs = list(csv.DictReader(open("/tmp/MANIFEST_450k.tsv"), delimiter="\t"))
    for i, r in enumerate(recs):
        r["_i"] = i
    log(f"manifest: {len(recs)} files")
    assert len(recs) == EXPECT_FILES, f"expected {EXPECT_FILES}, manifest has {len(recs)}"

    M = np.full((EXPECT_PROBES, len(recs)), np.nan, np.float32)
    canon_hash, canon_probes = None, None
    ok, bad, done = 0, [], 0
    with cf.ThreadPoolExecutor(WORKERS) as ex:
        for i, ph, payload in ex.map(fetch, recs):
            done += 1
            if payload is None:
                bad.append((recs[i]["file_id"], ph)); log(f"  FAILED {recs[i]['file_id']} {ph}")
            else:
                beta, probes = payload
                if canon_hash is None:
                    canon_hash, canon_probes = ph, probes
                if ph != canon_hash:
                    bad.append((recs[i]["file_id"], "PROBE_ORDER_MISMATCH"))
                    log(f"  PROBE MISMATCH {recs[i]['file_id']}")
                else:
                    M[:, i] = beta; ok += 1
            if done % 50 == 0:
                log(f"  {done}/{len(recs)} ok={ok} bad={len(bad)}")

    log(f"parsed ok={ok} bad={len(bad)}")
    if bad:
        log("FATAL: " + "; ".join(f"{a}:{b}" for a, b in bad[:20]))
        sys.exit(3)
    if ok != EXPECT_FILES:
        log(f"FATAL: {ok} of {EXPECT_FILES} parsed"); sys.exit(3)

    np.save("/mnt/scratch/beta_450k.npy", M)
    open("/mnt/scratch/probes.txt", "w").write("\n".join(canon_probes) + "\n")
    with open("/mnt/scratch/samples.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["col", "file_id", "sample", "case", "project", "sample_type"])
        for r in recs:
            w.writerow([r["_i"], r["file_id"], r["sample"], r["case"],
                        r["project"], r["sample_type"]])
    frac_nan = float(np.isnan(M).mean())
    open("/mnt/scratch/ACQUISITION.txt", "w").write(
        f"files={ok}\nprobes={EXPECT_PROBES}\nprobe_vector_sha256={canon_hash}\n"
        f"matrix_shape={M.shape}\nfraction_nan={frac_nan:.6f}\n"
        f"dtype=float32\ncolumn_order=see samples.tsv\n")
    log(f"matrix {M.shape}  NaN fraction {frac_nan:.4f}")

    for f in ("beta_450k.npy", "probes.txt", "samples.tsv", "ACQUISITION.txt"):
        subprocess.run(["gcloud", "storage", "cp", f"/mnt/scratch/{f}", f"{OUT}/{f}"], check=True)
    subprocess.run(["gcloud", "storage", "cp", "/tmp/MANIFEST_450k.tsv",
                    f"{OUT}/MANIFEST_450k.tsv"], check=True)
    log("uploaded; done")


if __name__ == "__main__":
    main()
