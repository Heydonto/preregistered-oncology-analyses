"""Count what the GSE68379 tar actually contains. Two extraction modes both stopped at 127 of
198 lung arrays, so the question is whether the archive itself ends early."""
import modal

BUCKET = "heydonto-quantara-lungcdx"
TAR = "data-request-2026-08/labels/lung/GSE68379/suppl/GSE68379_RAW.tar"
app = modal.App("quantara-gse68379-diag")
gcs = modal.Secret.from_name("gcs-key")
image = modal.Image.debian_slim(python_version="3.11").pip_install("google-cloud-storage==2.18.2")


def _client():
    import json, os
    open("/tmp/sa.json", "w").write(os.environ["GCP_SA_KEY"])
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/tmp/sa.json"
    from google.cloud import storage
    return storage.Client(project=json.loads(os.environ["GCP_SA_KEY"]).get("project_id"))


@app.function(image=image, secrets=[gcs], timeout=7200, cpu=4, memory=16384)
def diag() -> dict:
    import os, tarfile
    p = "/tmp/t.tar"
    bkt = _client().bucket(BUCKET)
    blob = bkt.blob(TAR)
    blob.reload()
    bkt.blob(TAR).download_to_filename(p)
    size = os.path.getsize(p)
    out = {"gcs_size": blob.size, "local_size": size, "sizes_match": blob.size == size}

    names, err = [], None
    try:
        with tarfile.open(p, "r") as tf:
            for m in tf:
                if m.isfile():
                    names.append(os.path.basename(m.name))
    except Exception as e:
        err = f"{type(e).__name__}: {e}"
    out["members_enumerated"] = len(names)
    out["enumeration_error"] = err
    out["distinct_gsms"] = len({n.split('_', 1)[0] for n in names if n.startswith('GSM')})
    out["last_member"] = names[-1] if names else None
    # is the archive terminated properly? tar ends with >=1024 zero bytes
    with open(p, "rb") as f:
        f.seek(-1024, 2)
        out["ends_with_zero_block"] = f.read(1024) == b"\x00" * 1024
    print(out, flush=True)
    return out


@app.local_entrypoint()
def main():
    print(diag.remote())
