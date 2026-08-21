"""Copy dinov2-large features from GCS into the Modal volume the MIL jobs read from.

The MIL scripts mount `nsclc-rwpr-features-tcga` and read `/vol/<dir>/<slide_id>.npz`. The
Phikon-v2 features already live under `/vol/tcga/`. This stages the second encoder alongside
them, so the comparison scripts differ only in which directory they point at:

    /vol/tcga/              Phikon-v2, TCGA          (1,053 slides, already present)
    /vol/dinov2-tcga/       dinov2-large, TCGA
    /vol/dinov2-external/   dinov2-large, EAGLE + HTMCP
    /vol/phikon-external/   Phikon-v2, EAGLE + HTMCP (needed for the site-probe control)

Runs in-cloud, so ~30 GB never crosses the local network. Idempotent: existing files are skipped.

Run:
  modal run spine/stage_dinov2_volume.py::stage --spec dinov2-tcga
  modal run spine/stage_dinov2_volume.py::stage --spec dinov2-external
  modal run spine/stage_dinov2_volume.py::stage --spec phikon-external
  modal run spine/stage_dinov2_volume.py::report
"""

import json
import os

import modal

BUCKET = "heydonto-quantara-lungcdx"

SPECS = {
    "dinov2-tcga": ("features/dinov2-large/tcga/", "dinov2-tcga"),
    "dinov2-external": ("features/dinov2-large/external/", "dinov2-external"),
    "phikon-external": ("features/phikon-v2/external/", "phikon-external"),
}

app = modal.App("nsclc-rwpr-stage-features")
vol = modal.Volume.from_name("nsclc-rwpr-features-tcga", create_if_missing=False)
gcs_secret = modal.Secret.from_name("gcs-key")
image = modal.Image.debian_slim(python_version="3.11").pip_install(
    "google-cloud-storage==2.18.2", "numpy==1.26.4",
)


def _gcs_client():
    from google.cloud import storage

    key = os.environ["GCP_SA_KEY"]
    p = "/tmp/sa.json"
    open(p, "w").write(key)
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = p
    return storage.Client(project=json.loads(key).get("project_id", "heydonto-425716"))


@app.function(image=image, volumes={"/vol": vol}, secrets=[gcs_secret], timeout=14400,
              cpu=8, memory=8192)
def stage_one(spec: str) -> dict:
    from concurrent.futures import ThreadPoolExecutor

    prefix, dest = SPECS[spec]
    out_dir = f"/vol/{dest}"
    os.makedirs(out_dir, exist_ok=True)

    cli = _gcs_client()
    bkt = cli.bucket(BUCKET)
    blobs = [b for b in bkt.list_blobs(prefix=prefix) if b.name.endswith(".npz")]
    have = {f for f in os.listdir(out_dir) if f.endswith(".npz")}
    todo = [b for b in blobs if os.path.basename(b.name) not in have]
    print(f"{spec}: {len(blobs)} in GCS, {len(have)} already staged, {len(todo)} to copy",
          flush=True)

    done = {"n": 0, "bytes": 0}

    def pull(b):
        fp = os.path.join(out_dir, os.path.basename(b.name))
        tmp = fp + ".part"
        b.download_to_filename(tmp)
        os.replace(tmp, fp)          # never leave a half-written npz at its final name
        done["n"] += 1
        done["bytes"] += os.path.getsize(fp)
        if done["n"] % 100 == 0:
            print(f"  {done['n']}/{len(todo)}", flush=True)

    with ThreadPoolExecutor(max_workers=16) as pool:
        list(pool.map(pull, todo))

    vol.commit()
    final = len([f for f in os.listdir(out_dir) if f.endswith(".npz")])
    print(f"{spec}: staged {done['n']}, {done['bytes']/2**30:.2f} GiB, "
          f"{final} npz total in {out_dir}", flush=True)
    return {"spec": spec, "copied": done["n"], "total": final,
            "gib": round(done["bytes"] / 2**30, 2)}


@app.function(image=image, volumes={"/vol": vol}, timeout=1800)
def report_fn() -> dict:
    out = {}
    for d in sorted(os.listdir("/vol")):
        p = f"/vol/{d}"
        if not os.path.isdir(p):
            continue
        fs = [f for f in os.listdir(p) if f.endswith(".npz")]
        out[d] = {"npz": len(fs),
                  "gib": round(sum(os.path.getsize(os.path.join(p, f)) for f in fs) / 2**30, 2)}
    return out


@app.local_entrypoint()
def stage(spec: str):
    print(stage_one.remote(spec))


@app.local_entrypoint()
def report():
    for k, v in report_fn.remote().items():
        print(f"  {k:24s} {v['npz']:5d} npz  {v['gib']:8.2f} GiB")
