"""Process GSE68379 lung cell-line IDATs into a beta matrix, on Modal.

SOURCE CHANGED 2026-08-20. This originally read GSE68379_RAW.tar from our GCS mirror. That copy
is byte-complete against both GCS and GEO's declared size (9,020,579,840) and ends with a valid
tar terminator, but enumerates only 1,302 of the 2,056 members that GEO's own filelist.txt
declares -- 649 of 1,028 cell lines, 127 of 198 lung. tarfile raises nothing; it simply returns
early. Two extraction modes gave the identical truncation, and a diagnostic run confirmed the
archive itself is the problem.

So the tar is abandoned and the 396 lung IDATs are fetched per sample from GEO instead, using the
filenames in filelist.txt. Slower per file, but each download is independently verifiable and a
corrupt one fails loudly.

GSE68379 ships 1,028 cell lines as raw 450k IDATs (2,056 files, 8.83 GB) with no processed
matrix. 198 of those lines are lung, and 187 of them also have GDSC2 dose-response, which is the
substrate for the only genuinely on-topic analysis available: whether DNA methylation predicts
drug sensitivity in lung cancer cell lines.

This extracts only the lung samples' IDATs from the tar, runs methylprep, and writes a beta
matrix. Selecting by GSM keeps ~1.7 GB of the 8.83 GB, though the tar must still be streamed in
full because tar is sequential.

Run:  modal run drugresist/process_idats.py::main
"""
import modal

BUCKET = "heydonto-quantara-lungcdx"
FILELIST = "data-request-2026-08/labels/lung/GSE68379/suppl/filelist.txt"
MATRIX = "data-request-2026-08/labels/lung/GSE68379/matrix/GSE68379_series_matrix.txt.gz"
OUT_PREFIX = "data-request-2026-08/labels/lung/GSE68379/processed/"

app = modal.App("quantara-gse68379-idat")
gcs = modal.Secret.from_name("gcs-key")
image = (modal.Image.debian_slim(python_version="3.11")
         # methylprep 1.7.1 calls DataFrame.append, removed in pandas 2.x, so pandas is held
         # below 2 and numpy pinned to a version that pandas 1.5.3 builds against.
         .pip_install("methylprep==1.7.1", "numpy==1.24.4", "pandas==1.5.3",
                      "scipy==1.10.1", "statsmodels==0.14.0",
                      "google-cloud-storage==2.18.2"))


def _client():
    import json, os
    key = os.environ["GCP_SA_KEY"]
    open("/tmp/sa.json", "w").write(key)
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/tmp/sa.json"
    from google.cloud import storage
    return storage.Client(project=json.loads(key).get("project_id", "heydonto-425716"))


def _lung_gsms(bkt):
    """GSM -> (cell_line, histology, cosmic_id) for every lung sample in the series matrix."""
    import gzip, io
    raw = bkt.blob(MATRIX).download_as_bytes()
    lines = gzip.open(io.BytesIO(raw), "rt", errors="ignore").read().split("\n")
    ch = [l for l in lines if l.startswith("!Sample_characteristics_ch1")]
    def vals(i):
        return [x.strip('"').split(": ", 1)[-1] for x in ch[i].split("\t")[1:]]
    gsm = [x.strip('"') for x in
           [l for l in lines if l.startswith("!Sample_geo_accession")][0].split("\t")[1:]]
    cell, site, hist, cid = vals(0), vals(1), vals(2), vals(3)
    return {g: (c, h, i) for g, c, s, h, i in zip(gsm, cell, site, hist, cid) if s == "lung"}


@app.function(image=image, secrets=[gcs], timeout=14400, cpu=8, memory=32768)
def run() -> dict:
    import io, os, shutil
    import numpy as np
    import pandas as pd

    cli = _client()
    bkt = cli.bucket(BUCKET)
    lung = _lung_gsms(bkt)
    print(f"lung samples in series matrix: {len(lung)}", flush=True)

    work = "/tmp/idats"
    os.makedirs(work, exist_ok=True)

    # exact filenames from GEO's own listing
    fl = bkt.blob(FILELIST).download_as_bytes().decode("utf-8", "ignore").split("\n")
    want = []
    for line in fl[2:]:
        p = line.split("\t")
        if len(p) > 1 and ".idat" in p[1]:
            g = p[1].split("_", 1)[0]
            if g in lung:
                want.append((g, p[1].strip()))
    print(f"lung idats to fetch from GEO: {len(want)}", flush=True)

    import urllib.request
    from concurrent.futures import ThreadPoolExecutor

    gsm_of, bad = {}, []

    def fetch(item):
        g, fn = item
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{g[:-3]}nnn/{g}/suppl/{fn}"
        stripped = fn.split("_", 1)[1]
        if stripped.endswith(".gz"):
            stripped = stripped[:-3]
        dest = os.path.join(work, stripped)
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    data = r.read()
                if data[:2] == b"\x1f\x8b":
                    import gzip as gz
                    data = gz.decompress(data)
                if len(data) < 1000:
                    raise ValueError(f"too small: {len(data)}")
                open(dest, "wb").write(data)
                return ("_".join(stripped.split("_")[:2]), g)
            except Exception as e:
                if attempt == 2:
                    bad.append((fn, f"{type(e).__name__}:{e}"))
                    return None
        return None

    with ThreadPoolExecutor(max_workers=12) as pool:
        for i, res in enumerate(pool.map(fetch, want), 1):
            if res:
                gsm_of[res[0]] = res[1]
            if i % 50 == 0:
                print(f"  fetched {i}/{len(want)}", flush=True)

    kept = len([f for f in os.listdir(work) if f.endswith(".idat")])
    print(f"fetched {kept} idats for {len(gsm_of)} arrays; {len(bad)} failed", flush=True)
    for b in bad[:10]:
        print("   failed:", b, flush=True)
    if len(gsm_of) < len(lung) - 5:
        raise RuntimeError(
            f"TRUNCATED INPUT: {len(gsm_of)} arrays from {len(lung)} lung samples. "
            "Refusing to process a silently-truncated input.")

    import methylprep
    print("running methylprep ...", flush=True)
    methylprep.run_pipeline(work, betas=True, make_sample_sheet=True,
                            export=False, batch_size=100)
    # methylprep writes beta_values*.pkl into the data dir
    import glob
    pkls = sorted(glob.glob(os.path.join(work, "beta_values*.pkl")))
    print("beta pickles:", pkls, flush=True)
    betas = pd.concat([pd.read_pickle(p) for p in pkls], axis=1) if len(pkls) > 1 \
        else pd.read_pickle(pkls[0])
    print("beta matrix:", betas.shape, flush=True)

    # relabel columns sentrix -> GSM, then attach cell line / histology / cosmic
    cols = {}
    for c in betas.columns:
        s = str(c)
        cols[c] = gsm_of.get(s, gsm_of.get(s.replace(".", "_"), s))
    betas = betas.rename(columns=cols)
    meta = pd.DataFrame(
        [{"gsm": g, "cell_line": lung[g][0], "histology": lung[g][1], "cosmic_id": lung[g][2]}
         for g in betas.columns if g in lung])
    print(f"labelled {len(meta)} of {betas.shape[1]} columns", flush=True)

    betas = betas[[c for c in betas.columns if c in lung]]
    buf = io.BytesIO()
    np.savez_compressed(buf,
                        betas=betas.to_numpy(dtype=np.float32),
                        probes=np.array(betas.index.astype(str), dtype=object),
                        gsm=np.array(betas.columns.astype(str), dtype=object))
    buf.seek(0)
    bkt.blob(OUT_PREFIX + "lung_betas.npz").upload_from_file(buf)
    bkt.blob(OUT_PREFIX + "lung_meta.csv").upload_from_string(
        meta.to_csv(index=False), content_type="text/csv")
    shutil.rmtree(work, ignore_errors=True)
    out = {"arrays": int(betas.shape[1]), "probes": int(betas.shape[0]),
           "labelled": int(len(meta)), "idats_skipped": len(bad)}
    print(out, flush=True)
    return out


@app.local_entrypoint()
def main():
    print(run.remote())
