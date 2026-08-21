"""Third, fourth and fifth encoders: UNI, H-optimus-0 and Virchow2 over all 1,186 slides.

WHY. R18 showed that subtype site-leakage reproduces under facebook/dinov2-large -- an encoder
trained on natural images -- so the signature is a property of the archive rather than of Phikon-v2.
R20 showed the external EAGLE-vs-HTMCP separation is encoder-independent in existence but not in
magnitude (1.000 vs 0.877). Both reports carry the same limitation: two encoders is not a survey,
and one of the two is not a pathology model at all. Nothing so far says whether the pattern holds
across *pathology* foundation models.

These three answer that. UNI is the sharpest addition:

                     arch          CLS dim  patch  reg  normalisation      corpus
  Phikon-v2          ViT-L/16      1024     16     0    ImageNet           pan-cancer histology
  dinov2-large       ViT-L/14      1024     14     0    ImageNet           LVD-142M natural
  UNI                ViT-L/16      1024     16     0    ImageNet           Mass-100K histology
  H-optimus-0        ViT-g/14      1536     14     4    its own            histology
  Virchow2           ViT-H/14      1280     14     4    ImageNet           histology

UNI matches Phikon-v2 on architecture, patch size, width and normalisation, differing only in
training corpus -- a tighter contrast than dinov2-large, which differed in patch size too.

TWO THINGS THAT ARE NO LONGER HELD CONSTANT, stated because R18's contrast was tighter.

  1. Normalisation. H-optimus-0 prescribes mean (0.707, 0.579, 0.704) / std (0.212, 0.230, 0.178),
     not ImageNet. Forcing ImageNet stats on it would cripple it and make the comparison
     meaningless, so each model gets its own prescribed preprocessing, read from its own
     pretrained_cfg. Normalisation is therefore a variable across these arms, not a constant.
  2. Capacity. 1024 / 1536 / 1280-dimensional outputs from ViT-L, ViT-g and ViT-H. The MIL projects
     whatever it is given down to 512 before pooling, so the head is identical, but the encoders are
     not the same size and a difference in leakage could in principle be a difference in capacity.

WHAT IS STILL HELD CONSTANT. Tiles. The tissue grid, Otsu threshold, 224 px at ~0.5 um/px, 40%
tissue floor, MPP fallback and resize filter are copied verbatim from the Phikon-v2 spine, so the
tile coordinates come out identical -- `verify` checks that against the stored coords rather than
assuming it. Also constant: the readout is the CLS token for every encoder.

ON THE CLS-ONLY READOUT. Virchow2's authors recommend concatenating the class token with the mean
of the patch tokens (2560-d). We take CLS only, because holding the readout constant is what makes
five encoders comparable, and because the quantity of interest is the *gap* between fold regimes
rather than absolute capability. This may understate Virchow2's absolute performance and the report
says so.

LICENCES. H-optimus-0 is Apache-2.0. UNI and Virchow2 are CC-BY-NC-ND-4.0, non-commercial academic
use, with a clause directing commercial entities to the corresponding author. Every arm records its
encoder's licence in the output npz so that any downstream use can be checked against it, and the
Apache-2.0 arm is the one to lead with in anything that could be read as commercial.

Run:
  modal run spine/encode_pathfm.py::verify --model uni
  modal run --detach spine/encode_pathfm.py::run_blocking --model uni
  ... then h-optimus-0, then virchow2
Monitor:
  modal app logs nsclc-rwpr-encode-pathfm
"""

import modal

BUCKET = "heydonto-quantara-lungcdx"
WSI_PREFIXES = ["wsi/TCGA-LUAD/", "wsi/TCGA-LUSC/", "wsi_external/eagle/", "wsi_external/htmcp/"]
PHIKON_ROOT = "features/phikon-v2/"

# key -> (hf repo, output subdir, CLS dim, licence, extra timm kwargs)
MODELS = {
    # UNI's config.json carries init_values at TOP level, not under model_args, so timm does not
    # apply it and builds the ViT without LayerScale -- then the checkpoint's blocks.N.ls1.gamma
    # keys have nowhere to go and load_state_dict raises. Passing it explicitly creates the module.
    # The value itself is irrelevant (the checkpoint overwrites it); only not-None matters.
    "uni": ("MahmoodLab/UNI", "uni", 1024, "CC-BY-NC-ND-4.0",
            {"init_values": 1.0, "dynamic_img_size": True}),
    "h-optimus-0": ("bioptimus/H-optimus-0", "h-optimus-0", 1536, "Apache-2.0",
                    {"init_values": 1e-5, "dynamic_img_size": False}),
    "virchow2": ("paige-ai/Virchow2", "virchow2", 1280, "CC-BY-NC-ND-4.0", {}),
}
OUT_ROOT = "features/"

# ---- copied verbatim from the Phikon-v2 spine; do not tune independently ----
TILE_PX = 224
TARGET_MPP = 0.5
ASSUMED_MPP = 0.25
TISSUE_FRAC_MIN = 0.40
READ_THREADS = 8
GPU_KIND = "L4"
MAX_CONTAINERS = 12
# batch shrinks with model size to keep L4 memory headroom; affects speed, not results
BATCH = {"uni": 128, "virchow2": 64, "h-optimus-0": 48}

app = modal.App("nsclc-rwpr-encode-pathfm")
gcs_secret = modal.Secret.from_name("gcs-key")
hf_secret = modal.Secret.from_name("hf-token")


def _cohort_of(rel_path: str) -> str:
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


def _hf_login():
    """The secret's key name is not guaranteed; accept any of the usual spellings."""
    import os

    for k in ("HF_TOKEN", "HUGGING_FACE_HUB_TOKEN", "HUGGINGFACE_TOKEN", "HF_API_TOKEN"):
        v = os.environ.get(k)
        if v:
            os.environ["HF_TOKEN"] = v
            os.environ["HUGGING_FACE_HUB_TOKEN"] = v
            return True
    return False


image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("openslide-tools")
    .pip_install(
        "openslide-python==1.4.1",
        "pillow==10.4.0",
        "numpy==1.26.4",
        "torch==2.4.1",
        "timm==1.0.11",
        "huggingface-hub==0.26.2",
        "google-cloud-storage==2.18.2",
    )
)


def _build(model_key):
    """Create the encoder and return (model, mean, std, cls_dim). Normalisation comes from the
    model's own pretrained_cfg, never assumed."""
    import timm
    import torch

    repo, _, dim, _, kwargs = MODELS[model_key]
    _hf_login()
    if model_key == "virchow2":
        from timm.layers import SwiGLUPacked

        m = timm.create_model(f"hf-hub:{repo}", pretrained=True,
                              mlp_layer=SwiGLUPacked, act_layer=torch.nn.SiLU, **kwargs)
    else:
        m = timm.create_model(f"hf-hub:{repo}", pretrained=True, **kwargs)
    cfg = timm.data.resolve_model_data_config(m)
    mean, std = cfg["mean"], cfg["std"]
    print(f"{model_key}: {type(m).__name__} mean={mean} std={std}", flush=True)
    return m, mean, std, dim


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
    vb = w0 * w1 * (mu0 - mu1) ** 2
    vb[~np.isfinite(vb)] = 0
    return int(np.argmax(vb))


def _tissue_grid(slide, tile_l0):
    import numpy as np

    w0, h0 = slide.dimensions
    thumb = slide.get_thumbnail((2048, 2048)).convert("RGB")
    arr = np.asarray(thumb)
    mx = arr.max(axis=2).astype(np.float32)
    mn = arr.min(axis=2).astype(np.float32)
    sat = np.where(mx > 0, (mx - mn) / np.maximum(mx, 1.0), 0.0)
    sat_u8 = (sat * 255).astype(np.uint8)
    mask = (sat_u8 > _otsu_threshold(sat_u8)) & (mx > 40)

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

    size = max(1, int(round(tile_l0 / level_ds)))
    img = slide.read_region((x, y), level, (size, size)).convert("RGB")
    if img.size != (TILE_PX, TILE_PX):
        img = img.resize((TILE_PX, TILE_PX), Image.BILINEAR)
    import numpy as np

    return np.asarray(img)


# ------------------------------------------------------------------- worker
@app.function(image=image, gpu=GPU_KIND, cpu=4, memory=16384, timeout=14400,
              retries=modal.Retries(max_retries=2, initial_delay=10.0),
              max_containers=MAX_CONTAINERS, secrets=[gcs_secret, hf_secret])
def encode_slide(rel_path: str, model_key: str) -> dict:
    import io
    import os
    import time
    from concurrent.futures import ThreadPoolExecutor

    import numpy as np
    import openslide
    import torch

    t0 = time.time()
    repo, subdir, dim, lic, _ = MODELS[model_key]
    slide_id = os.path.basename(rel_path).split(".")[0]
    cohort = _cohort_of(rel_path)
    row = {"slide_id": slide_id, "cohort": cohort, "model": model_key,
           "slide_path": f"gs://{BUCKET}/{rel_path}", "n_tiles": 0, "status": "",
           "wall_s": 0.0, "mpp": "", "mpp_assumed": False}

    cli = _gcs_client()
    bkt = cli.bucket(BUCKET)
    out_blob = bkt.blob(f"{OUT_ROOT}{subdir}/{cohort}/{slide_id}.npz")
    if out_blob.exists():
        row.update(status="skipped_exists", wall_s=round(time.time() - t0, 1))
        return row

    local = f"/tmp/{os.path.basename(rel_path)}"
    try:
        bkt.blob(rel_path).download_to_filename(local)
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
        model, mean_v, std_v, dim = _build(model_key)
        model = model.eval().to(device).half()
        mean = torch.tensor(mean_v, device=device).view(1, 3, 1, 1)
        std = torch.tensor(std_v, device=device).view(1, 3, 1, 1)

        bs = BATCH[model_key]
        embs = np.empty((n, dim), dtype=np.float16)
        pool = ThreadPoolExecutor(max_workers=READ_THREADS)
        with torch.inference_mode():
            for b0 in range(0, n, bs):
                batch = coords[b0:b0 + bs]
                tiles = list(pool.map(
                    lambda c: _read_tile(slide, c[0], c[1], tile_l0, level, level_ds), batch))
                px = torch.from_numpy(np.stack(tiles)).to(device)
                px = px.permute(0, 3, 1, 2).float().div_(255.0).sub_(mean).div_(std).half()
                out = model(px)
                # global_pool="token" models return (B, D); Virchow2 (global_pool="") returns
                # (B, tokens, D) with layout [CLS, reg x4, patches...]. CLS is index 0 either way.
                if out.dim() == 3:
                    out = out[:, 0]
                if out.shape[1] != dim:
                    raise RuntimeError(f"{model_key}: got dim {out.shape[1]}, expected {dim}")
                embs[b0:b0 + len(batch)] = out.float().cpu().numpy().astype(np.float16)
        pool.shutdown()

        if np.isnan(embs).any():
            row["status"] = "error_nan_embeddings"
            return row

        buf = io.BytesIO()
        np.savez_compressed(
            buf, embeddings=embs, coords=np.asarray(coords, dtype=np.int32),
            mpp=np.float64(mpp), mpp_assumed=np.bool_(row["mpp_assumed"]),
            tile_count=np.int64(n), tile_size_level0=np.int64(tile_l0),
            target_mpp=np.float64(TARGET_MPP), slide_path=np.str_(row["slide_path"]),
            encoder=np.str_(repo), encoder_licence=np.str_(lic),
            readout=np.str_("cls_token"),
            norm_mean=np.asarray(mean_v, np.float64), norm_std=np.asarray(std_v, np.float64))
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
@app.function(image=image, gpu=GPU_KIND, cpu=4, memory=16384, timeout=3600,
              secrets=[gcs_secret, hf_secret])
def check_one(rel_path: str, model_key: str) -> dict:
    """Recompute the grid, compare to the Phikon-v2 coords, and confirm the encoder's real
    output width and normalisation on one real batch."""
    import io
    import os

    import numpy as np
    import openslide
    import torch

    repo, _, dim, lic, _ = MODELS[model_key]
    slide_id = os.path.basename(rel_path).split(".")[0]
    cohort = _cohort_of(rel_path)
    out = {"slide_id": slide_id, "model": model_key, "repo": repo, "licence": lic,
           "coords_match": False, "note": ""}
    cli = _gcs_client()
    bkt = cli.bucket(BUCKET)
    ref = bkt.blob(f"{PHIKON_ROOT}{cohort}/{slide_id}.npz")
    if not ref.exists():
        out["note"] = "no phikon reference"
        return out

    local = f"/tmp/{os.path.basename(rel_path)}"
    try:
        bkt.blob(rel_path).download_to_filename(local)
        slide = openslide.OpenSlide(local)
        mpp_raw = slide.properties.get(openslide.PROPERTY_NAME_MPP_X)
        mpp = float(mpp_raw) if mpp_raw else ASSUMED_MPP
        tile_l0 = int(round(TILE_PX * TARGET_MPP / mpp))
        mine = np.asarray(_tissue_grid(slide, tile_l0), dtype=np.int32)
        d = np.load(io.BytesIO(ref.download_as_bytes()), allow_pickle=True)
        theirs = d["coords"]
        out["n_mine"] = int(mine.shape[0])
        out["n_theirs"] = int(theirs.shape[0])
        out["coords_match"] = bool(mine.shape == theirs.shape and (mine == theirs).all())

        model, mean_v, std_v, dim = _build(model_key)
        model = model.eval().to("cuda").half()
        mean = torch.tensor(mean_v, device="cuda").view(1, 3, 1, 1)
        std = torch.tensor(std_v, device="cuda").view(1, 3, 1, 1)
        tiles = np.stack([_read_tile(slide, *mine[i], tile_l0, level=0,
                                    level_ds=slide.level_downsamples[0]) for i in range(2)])
        with torch.inference_mode():
            px = torch.from_numpy(tiles).to("cuda").permute(0, 3, 1, 2).float()
            px = px.div_(255.0).sub_(mean).div_(std).half()
            o = model(px)
        out["raw_out_shape"] = list(o.shape)
        if o.dim() == 3:
            o = o[:, 0]
        out["cls_dim"] = int(o.shape[1])
        out["cls_dim_matches_config"] = int(o.shape[1]) == dim
        out["norm_mean"] = [round(float(x), 6) for x in mean_v]
        out["finite"] = bool(torch.isfinite(o).all())
    except Exception as e:
        out["note"] = f"error:{type(e).__name__}:{str(e)[:200]}"
    finally:
        if os.path.exists(local):
            os.remove(local)
    return out


# ------------------------------------------------------------------- driver
@app.function(image=image, cpu=2, memory=4096, timeout=86400, secrets=[gcs_secret, hf_secret])
def driver(model_key: str, limit: int = 0, only: list = None) -> list:
    import csv
    import io

    repo, subdir, dim, lic, _ = MODELS[model_key]
    cli = _gcs_client()
    bkt = cli.bucket(BUCKET)
    slides = []
    for p in WSI_PREFIXES:
        slides += [b.name for b in bkt.list_blobs(prefix=p) if b.name.endswith(".svs")]
    slides.sort()
    if only:
        slides = [s for s in slides if any(o in s for o in only)]
    done = set()
    for c in ("tcga", "external"):
        done |= {b.name.split("/")[-1][:-4] for b in
                 bkt.list_blobs(prefix=f"{OUT_ROOT}{subdir}/{c}/") if b.name.endswith(".npz")}
    todo = [s for s in slides if s.split("/")[-1].split(".")[0] not in done]
    if limit:
        todo = todo[:limit]
    print(f"[{model_key}] slides={len(slides)} done={len(done)} todo={len(todo)}", flush=True)

    fields = ["slide_id", "cohort", "model", "slide_path", "n_tiles", "status", "wall_s",
              "mpp", "mpp_assumed"]
    rows, ct = [], 0
    for row in encode_slide.starmap([(t, model_key) for t in todo], order_outputs=False):
        rows.append(row)
        ct += 1
        print(f"[{model_key} {ct}/{len(todo)}] {row['slide_id']} {row['cohort']} "
              f"{row['status']} n_tiles={row['n_tiles']} wall_s={row['wall_s']}", flush=True)
        if ct % 25 == 0 or ct == len(todo):
            _write_manifests(bkt, subdir, fields, rows)
    if rows:
        _write_manifests(bkt, subdir, fields, rows)
    return rows


def _write_manifests(bkt, subdir, fields, new_rows):
    import csv
    import io

    for cohort in ("tcga", "external"):
        blob = bkt.blob(f"{OUT_ROOT}{subdir}/{cohort}/_manifest.csv")
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
VERIFY = ["TCGA-05-4244-01Z-00-DX1", "TCGA-18-3406-01Z-00-DX1"]


@app.function(image=image, cpu=2, memory=4096, timeout=7200, secrets=[gcs_secret, hf_secret])
def verify_driver(model_key: str) -> list:
    cli = _gcs_client()
    bkt = cli.bucket(BUCKET)
    tcga = sorted(b.name for p in ("wsi/TCGA-LUAD/", "wsi/TCGA-LUSC/")
                  for b in bkt.list_blobs(prefix=p)
                  if b.name.endswith(".svs") and any(v in b.name for v in VERIFY))
    ext = sorted(b.name for b in bkt.list_blobs(prefix="wsi_external/eagle/")
                 if b.name.endswith(".svs"))[:1]
    picks = tcga[:2] + ext
    print("verifying:", picks, flush=True)
    return list(check_one.starmap([(p, model_key) for p in picks]))


@app.local_entrypoint()
def verify(model: str):
    rows = verify_driver.remote(model)
    ok = 0
    for r in rows:
        print(r)
        ok += bool(r.get("coords_match")) and bool(r.get("cls_dim_matches_config"))
    print(f"\n{model}: tiles identical AND width as configured on {ok}/{len(rows)} slides")
    if ok != len(rows):
        print("STOP: fix before encoding -- the cross-encoder comparison would be confounded.")


@app.local_entrypoint()
def run_blocking(model: str, limit: int = 0):
    rows = driver.remote(model, limit=limit)
    ok = sum(1 for r in rows if r.get("status") == "ok")
    print(f"[{model}] encoded ok={ok} of {len(rows)}")
    for r in rows:
        if r.get("status") not in ("ok", "skipped_exists"):
            print("  NOT OK:", r)
