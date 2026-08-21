"""Second-encoder arm: re-encode every TCGA + external slide with facebook/dinov2-large.

WHY. R12-R15 measured a large site signature in whole-slide features: subtype AUROC inflates
by +0.171 and the six methylation targets by +0.085 in mean rho when folds share tissue-source
sites instead of being site-disjoint. Every one of those numbers came from Phikon-v2 features.
An obvious objection is that the signature is a property of THAT encoder -- something about
pan-cancer histology pretraining that happens to encode institutional stain and scanner
character -- rather than a property of the archive. If so, the finding is about Owkin's model,
not about federated pathology, and Paper 1's claim would be far narrower than stated.

THE CONTRAST. facebook/dinov2-large is the tightest available control:

                        Phikon-v2                 facebook/dinov2-large
  architecture          Dinov2Model               Dinov2Model
  width / depth         1024 / 24                 1024 / 24
  SSL algorithm         DINOv2                    DINOv2
  readout               CLS token                 CLS token
  preprocessing         ImageNet mean/std         ImageNet mean/std
  training corpus       pan-cancer histology      LVD-142M natural images
  patch size            16                        14

Same class, same capacity, same self-supervised recipe, same readout, same normalisation. The
training corpus is what changes. So if the site signature survives here it cannot be a quirk of
histology pretraining, and if it vanishes it was never about the archive. Neither model is gated,
so this needs no licence acceptance.

WHAT IS HELD FIXED. Comparability is the whole point, so the tiling is not re-designed: the
tissue grid, the Otsu threshold, the 224 px tile at ~0.5 um/px, the 40% tissue floor, the MPP
fallback and the batch size are copied verbatim from encode_tcga.py / encode_external.py. That
code is deterministic given a slide, so the tile coordinates should come out bit-identical to the
Phikon-v2 run. `verify` checks that against the stored coords rather than assuming it -- if the
grids differ, the comparison is confounded by tiling and the numbers must not be compared.

Cost basis, measured from the Phikon-v2 manifests rather than guessed: 1,053 TCGA slides took
38.09 GPU-h and 129 external slides 3.01 GPU-h. dinov2-large sees 256 tokens per tile against
Phikon-v2's 196, so expect roughly 1.35x, i.e. ~55 GPU-h on L4.

Run:
  modal run spine/encode_dinov2.py::verify              # tile-grid equality check, 3 slides
  modal run --detach spine/encode_dinov2.py::run_blocking
Monitor:
  modal app logs nsclc-rwpr-encode-dinov2
"""

import modal

BUCKET = "heydonto-quantara-lungcdx"
WSI_PREFIXES = ["wsi/TCGA-LUAD/", "wsi/TCGA-LUSC/", "wsi_external/eagle/", "wsi_external/htmcp/"]
OUT_ROOT = "features/dinov2-large/"
PHIKON_ROOT = "features/phikon-v2/"

ENCODER_NAME = "facebook/dinov2-large"
ENCODER_REVISION = "47b73eefe95e8d44ec3623f8890bd894b6ea2d6c"
EMB_DIM = 1024

# ---- copied verbatim from the Phikon-v2 spine; do not tune independently ----
TILE_PX = 224
TARGET_MPP = 0.5
ASSUMED_MPP = 0.25
TISSUE_FRAC_MIN = 0.40
BATCH_SIZE = 128  # 256 tokens/tile vs 196; smaller batch keeps L4 memory headroom
READ_THREADS = 8
GPU_KIND = "L4"
MAX_CONTAINERS = 12

app = modal.App("nsclc-rwpr-encode-dinov2")
gcs_secret = modal.Secret.from_name("gcs-key")


def _cohort_of(rel_path: str) -> str:
    """wsi/TCGA-* -> tcga ; wsi_external/* -> external."""
    return "tcga" if rel_path.startswith("wsi/") else "external"


def _gcs_client():
    import json
    import os

    key = os.environ["GCP_SA_KEY"]
    path = "/tmp/sa.json"
    with open(path, "w") as f:
        f.write(key)
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = path
    from google.cloud import storage

    proj = json.loads(key).get("project_id", "heydonto-425716")
    return storage.Client(project=proj)


def _bake_model():
    """Image-build step: pull the pinned dinov2-large revision into /model."""
    import os

    from huggingface_hub import snapshot_download

    p = snapshot_download(
        ENCODER_NAME,
        revision=ENCODER_REVISION,
        local_dir="/model",
        allow_patterns=["config.json", "model.safetensors", "preprocessor_config.json"],
    )
    print("model staged:", sorted(os.listdir(p)))


image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("openslide-tools")
    .pip_install(
        "openslide-python==1.4.1",
        "pillow==10.4.0",
        "numpy==1.26.4",
        "torch==2.4.1",
        "transformers==4.44.2",
        "huggingface-hub==0.24.6",
        "google-cloud-storage==2.18.2",
    )
    .run_function(_bake_model)
)


# ------------------------------------------------------- tiling (verbatim copy)
def _otsu_threshold(channel):
    import numpy as np

    hist = np.bincount(channel.ravel(), minlength=256).astype(np.float64)
    total = hist.sum()
    if total == 0:
        return 0
    w = np.cumsum(hist)
    m = np.cumsum(hist * np.arange(256))
    w0 = w / total
    w1 = 1.0 - w0
    with np.errstate(divide="ignore", invalid="ignore"):
        mu0 = m / w
        mu1 = (m[-1] - m) / (total - w)
    var_between = w0 * w1 * (mu0 - mu1) ** 2
    var_between[~np.isfinite(var_between)] = 0
    return int(np.argmax(var_between))


def _tissue_grid(slide, tile_l0):
    import numpy as np

    w0, h0 = slide.dimensions
    thumb = slide.get_thumbnail((2048, 2048)).convert("RGB")
    arr = np.asarray(thumb)
    mx = arr.max(axis=2).astype(np.float32)
    mn = arr.min(axis=2).astype(np.float32)
    sat = np.where(mx > 0, (mx - mn) / np.maximum(mx, 1.0), 0.0)
    sat_u8 = (sat * 255).astype(np.uint8)
    thr = _otsu_threshold(sat_u8)
    mask = sat_u8 > thr
    mask &= mx > 40

    sx = mask.shape[1] / w0
    sy = mask.shape[0] / h0
    coords = []
    for y in range(0, h0 - tile_l0 + 1, tile_l0):
        my0, my1 = int(y * sy), max(int((y + tile_l0) * sy), int(y * sy) + 1)
        row = mask[my0:my1]
        for x in range(0, w0 - tile_l0 + 1, tile_l0):
            mx0, mx1 = int(x * sx), max(int((x + tile_l0) * sx), int(x * sx) + 1)
            if row[:, mx0:mx1].mean() >= TISSUE_FRAC_MIN:
                coords.append((x, y))
    return coords


def _read_tile(slide, x, y, tile_l0, level, level_ds):
    from PIL import Image

    size_at_level = max(1, int(round(tile_l0 / level_ds)))
    img = slide.read_region((x, y), level, (size_at_level, size_at_level)).convert("RGB")
    if img.size != (TILE_PX, TILE_PX):
        img = img.resize((TILE_PX, TILE_PX), Image.BILINEAR)
    import numpy as np

    return np.asarray(img)


# ------------------------------------------------------------------- worker
@app.function(
    image=image,
    gpu=GPU_KIND,
    cpu=4,
    memory=12288,
    timeout=10800,
    retries=modal.Retries(max_retries=2, initial_delay=10.0),
    max_containers=MAX_CONTAINERS,
    secrets=[gcs_secret],
)
def encode_slide(rel_path: str) -> dict:
    import io
    import os
    import time
    from concurrent.futures import ThreadPoolExecutor

    import numpy as np
    import openslide
    import torch

    t0 = time.time()
    slide_id = os.path.basename(rel_path).split(".")[0]
    cohort = _cohort_of(rel_path)
    row = {
        "slide_id": slide_id,
        "cohort": cohort,
        "slide_path": f"gs://{BUCKET}/{rel_path}",
        "n_tiles": 0,
        "status": "",
        "wall_s": 0.0,
        "mpp": "",
        "mpp_assumed": False,
    }

    client = _gcs_client()
    bucket = client.bucket(BUCKET)
    out_blob = bucket.blob(f"{OUT_ROOT}{cohort}/{slide_id}.npz")
    if out_blob.exists():
        row.update(status="skipped_exists", wall_s=round(time.time() - t0, 1))
        return row

    local = f"/tmp/{os.path.basename(rel_path)}"
    try:
        bucket.blob(rel_path).download_to_filename(local)
        slide = openslide.OpenSlide(local)

        mpp_raw = slide.properties.get(openslide.PROPERTY_NAME_MPP_X)
        if mpp_raw:
            mpp = float(mpp_raw)
        else:
            mpp = ASSUMED_MPP
            row["mpp_assumed"] = True
        row["mpp"] = mpp

        tile_l0 = int(round(TILE_PX * TARGET_MPP / mpp))
        level = slide.get_best_level_for_downsample(tile_l0 / TILE_PX + 0.01)
        level_ds = slide.level_downsamples[level]

        coords = _tissue_grid(slide, tile_l0)
        n = len(coords)
        row["n_tiles"] = n
        if n == 0:
            row["status"] = "empty_no_tissue"
            return row

        device = "cuda"
        from transformers import AutoModel

        model = AutoModel.from_pretrained("/model", torch_dtype=torch.float16).eval().to(device)
        mean = torch.tensor([0.485, 0.456, 0.406], device=device).view(1, 3, 1, 1)
        std = torch.tensor([0.229, 0.224, 0.225], device=device).view(1, 3, 1, 1)

        embs = np.empty((n, EMB_DIM), dtype=np.float16)
        pool = ThreadPoolExecutor(max_workers=READ_THREADS)
        with torch.inference_mode():
            for b0 in range(0, n, BATCH_SIZE):
                batch = coords[b0 : b0 + BATCH_SIZE]
                tiles = list(
                    pool.map(lambda c: _read_tile(slide, c[0], c[1], tile_l0, level, level_ds), batch)
                )
                px = torch.from_numpy(np.stack(tiles)).to(device)
                px = px.permute(0, 3, 1, 2).float().div_(255.0).sub_(mean).div_(std).half()
                out = model(pixel_values=px)
                embs[b0 : b0 + len(batch)] = out.last_hidden_state[:, 0].cpu().numpy().astype(np.float16)
        pool.shutdown()

        if np.isnan(embs).any():
            row["status"] = "error_nan_embeddings"
            return row

        buf = io.BytesIO()
        np.savez_compressed(
            buf,
            embeddings=embs,
            coords=np.asarray(coords, dtype=np.int32),
            mpp=np.float64(mpp),
            mpp_assumed=np.bool_(row["mpp_assumed"]),
            tile_count=np.int64(n),
            tile_size_level0=np.int64(tile_l0),
            target_mpp=np.float64(TARGET_MPP),
            slide_path=np.str_(row["slide_path"]),
            encoder=np.str_(ENCODER_NAME),
            encoder_revision=np.str_(ENCODER_REVISION),
        )
        buf.seek(0)
        out_blob.upload_from_file(buf, content_type="application/octet-stream")
        row["status"] = "ok"
    except Exception as e:
        row["status"] = f"error:{type(e).__name__}:{str(e)[:200]}"
    finally:
        row["wall_s"] = round(time.time() - t0, 1)
        if os.path.exists(local):
            os.remove(local)
    return row


# --------------------------------------------------- tile-grid equality check
@app.function(image=image, gpu=GPU_KIND, cpu=4, memory=12288, timeout=3600,
              secrets=[gcs_secret])
def check_coords(rel_path: str) -> dict:
    """Recompute the tissue grid and compare to the coords stored by the Phikon-v2 run.

    The comparison in R18 only isolates the encoder if both runs saw the same tiles. This
    proves it on real slides instead of trusting that the copied code is deterministic."""
    import io
    import os

    import numpy as np
    import openslide

    slide_id = os.path.basename(rel_path).split(".")[0]
    cohort = _cohort_of(rel_path)
    out = {"slide_id": slide_id, "cohort": cohort, "match": False, "note": ""}

    client = _gcs_client()
    bucket = client.bucket(BUCKET)
    ref = bucket.blob(f"{PHIKON_ROOT}{cohort}/{slide_id}.npz")
    if not ref.exists():
        out["note"] = "no phikon-v2 reference npz"
        return out

    local = f"/tmp/{os.path.basename(rel_path)}"
    try:
        bucket.blob(rel_path).download_to_filename(local)
        slide = openslide.OpenSlide(local)
        mpp_raw = slide.properties.get(openslide.PROPERTY_NAME_MPP_X)
        mpp = float(mpp_raw) if mpp_raw else ASSUMED_MPP
        tile_l0 = int(round(TILE_PX * TARGET_MPP / mpp))
        mine = np.asarray(_tissue_grid(slide, tile_l0), dtype=np.int32)

        d = np.load(io.BytesIO(ref.download_as_bytes()), allow_pickle=True)
        theirs = d["coords"]
        out["n_mine"] = int(mine.shape[0])
        out["n_theirs"] = int(theirs.shape[0])
        out["tile_size_level0_match"] = int(d["tile_size_level0"]) == tile_l0
        out["match"] = bool(mine.shape == theirs.shape and (mine == theirs).all())
        if not out["match"]:
            inter = len({tuple(t) for t in mine} & {tuple(t) for t in theirs})
            out["note"] = f"jaccard_intersection={inter}"
    except Exception as e:
        out["note"] = f"error:{type(e).__name__}:{str(e)[:160]}"
    finally:
        if os.path.exists(local):
            os.remove(local)
    return out


# ------------------------------------------------------------------- driver
@app.function(image=image, cpu=2, memory=4096, timeout=86400, secrets=[gcs_secret])
def driver(limit: int = 0, only: list = None, cohort: str = "") -> list:
    import csv
    import io

    client = _gcs_client()
    bucket = client.bucket(BUCKET)

    slides = []
    for prefix in WSI_PREFIXES:
        slides += [b.name for b in bucket.list_blobs(prefix=prefix) if b.name.endswith(".svs")]
    slides.sort()
    if cohort:
        slides = [s for s in slides if _cohort_of(s) == cohort]
    if only:
        slides = [s for s in slides if any(o in s for o in only)]

    done = set()
    for c in ("tcga", "external"):
        done |= {
            b.name.split("/")[-1][: -len(".npz")]
            for b in bucket.list_blobs(prefix=f"{OUT_ROOT}{c}/")
            if b.name.endswith(".npz")
        }
    todo = [s for s in slides if s.split("/")[-1].split(".")[0] not in done]
    if limit:
        todo = todo[:limit]
    print(f"slides={len(slides)} done={len(done)} todo={len(todo)}", flush=True)

    fields = ["slide_id", "cohort", "slide_path", "n_tiles", "status", "wall_s", "mpp",
              "mpp_assumed"]
    rows, done_ct = [], 0
    for row in encode_slide.map(todo, order_outputs=False):
        rows.append(row)
        done_ct += 1
        print(f"[{done_ct}/{len(todo)}] {row['slide_id']} {row['cohort']} {row['status']} "
              f"n_tiles={row['n_tiles']} wall_s={row['wall_s']}", flush=True)
        if done_ct % 25 == 0 or done_ct == len(todo):
            _write_manifests(bucket, fields, rows)
    if rows:
        _write_manifests(bucket, fields, rows)
    return rows


def _write_manifests(bucket, fields, new_rows):
    """One _manifest.csv per cohort, merged with whatever is already there."""
    import csv
    import io

    for cohort in ("tcga", "external"):
        blob = bucket.blob(f"{OUT_ROOT}{cohort}/_manifest.csv")
        existing = {}
        if blob.exists():
            for r in csv.DictReader(io.StringIO(blob.download_as_text())):
                existing[r["slide_id"]] = r
        touched = False
        for r in new_rows:
            if r["cohort"] != cohort:
                continue
            if r["status"] != "skipped_exists" or r["slide_id"] not in existing:
                existing[r["slide_id"]] = {k: str(r.get(k, "")) for k in fields}
                touched = True
        if not touched and existing:
            continue
        out = io.StringIO()
        w = csv.DictWriter(out, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for sid in sorted(existing):
            w.writerow({k: existing[sid].get(k, "") for k in fields})
        blob.upload_from_string(out.getvalue(), content_type="text/csv")


# --------------------------------------------------------------- entrypoints
VERIFY_SLIDES = [
    "TCGA-05-4244-01Z-00-DX1",  # LUAD, MPP from metadata
    "TCGA-18-3406-01Z-00-DX1",  # LUSC
]


@app.function(image=image, cpu=2, memory=4096, timeout=7200, secrets=[gcs_secret])
def verify_driver() -> list:
    """Pick 2 TCGA + 2 external slides and check tile-grid equality on each."""
    client = _gcs_client()
    bucket = client.bucket(BUCKET)
    ext = sorted(
        b.name for p in ("wsi_external/eagle/", "wsi_external/htmcp/")
        for b in bucket.list_blobs(prefix=p) if b.name.endswith(".svs")
    )
    tcga = sorted(
        b.name for p in ("wsi/TCGA-LUAD/", "wsi/TCGA-LUSC/")
        for b in bucket.list_blobs(prefix=p)
        if b.name.endswith(".svs") and any(v in b.name for v in VERIFY_SLIDES)
    )
    picks = tcga[:2] + ext[:1] + ext[-1:]
    print("verifying:", picks, flush=True)
    return list(check_coords.map(picks))


@app.local_entrypoint()
def verify():
    rows = verify_driver.remote()
    ok = 0
    for r in rows:
        print(r)
        ok += bool(r.get("match"))
    print(f"\ntile-grid identical on {ok}/{len(rows)} slides")
    if ok != len(rows):
        print("STOP: tiling differs from the Phikon-v2 run. The encoder comparison would be "
              "confounded by tile selection -- fix before encoding.")


@app.local_entrypoint()
def smoke():
    rows = driver.remote(only=VERIFY_SLIDES)
    for r in rows:
        print(r)


@app.local_entrypoint()
def run_blocking(cohort: str = "", limit: int = 0):
    """Blocking: driver.spawn() does not survive `modal run`'s ephemeral app teardown."""
    rows = driver.remote(limit=limit, cohort=cohort)
    ok = sum(1 for r in rows if r.get("status") == "ok")
    print(f"encoded ok={ok} of {len(rows)}")
    for r in rows:
        if r.get("status") not in ("ok", "skipped_exists"):
            print("  NOT OK:", r)
