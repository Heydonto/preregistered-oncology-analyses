#!/usr/bin/env python3
"""R18 Modal job — is the site signature a property of the archive or of Phikon-v2?

R12-R15 measured site leakage in whole-slide features: subtype AUROC inflates by +0.171 and the
six methylation targets by +0.0849 in mean rho when folds share tissue-source sites instead of
being site-disjoint. Every one of those numbers came from Phikon-v2. If the effect is really a
property of that encoder -- pan-cancer histology pretraining happening to encode institutional
stain and scanner character -- then Paper 1's claim is about Owkin's model, not about federated
pathology.

This script runs the identical measurement over a second encoder, facebook/dinov2-large: same
Dinov2Model class, same 1024-d width, same 24 layers, same DINOv2 self-supervision, same CLS
readout, same ImageNet normalisation, trained on natural images instead of histology. Tiles are
byte-identical between the two feature sets (verified by spine/encode_dinov2.py::verify), so the
encoder is the only thing that changes.

Nothing else moves either. The architecture, TILE_CAP, LR, WD, EPOCHS, BATCH and SEED are copied
from mil_embed.py / mil_meth.py, and the fold assignments are the same JSON files the published
runs consumed -- grouped folds hold out 11-14 whole sites per fold, random folds mix 39-49.

So the Phikon-v2 arms here are a REPRODUCTION CONTROL: they must recover the published numbers
before the dinov2 arms mean anything. If they do not, the harness changed and the comparison is
void.

  mode=subtype   BCE, one head per target, reports AUROC for subtype and KEAP1
  mode=meth      MSE over 6 targets, train-fold z-scoring, reports Spearman per target

Run (four arms per encoder pair):
  modal run armD-meth/mil_encoder_compare.py::main --feat-dir tcga --folds grouped --mode subtype
  modal run armD-meth/mil_encoder_compare.py::main --feat-dir tcga --folds random  --mode subtype
  modal run armD-meth/mil_encoder_compare.py::main --feat-dir dinov2-tcga --folds grouped --mode subtype
  ... and the same four with --mode meth
"""
import json
import os

import modal

BUCKET = "heydonto-quantara-lungcdx"
IN_GROUPED = "nsclc-rwpr-study/armD-meth/mil_input.json"
IN_RANDOM = "nsclc-rwpr-study/armD-meth/mil_input_random.json"
TGT_BLOB = "nsclc-rwpr-study/armD-meth/meth_targets.json"
OUT_PREFIX = "nsclc-rwpr-study/armD-meth/r18/"

# Feature width is a property of the encoder, not a constant. Phikon-v2, dinov2-large and UNI are
# 1024-d; H-optimus-0 is 1536 (ViT-g) and Virchow2 is 1280 (ViT-H). The width is read from the
# features and checked against this table, so a mis-staged directory fails loudly instead of
# training a projection layer against the wrong input size.
EXPECTED_DIM = {"tcga": 1024, "dinov2-tcga": 1024, "uni-tcga": 1024,
                "hopt-tcga": 1536, "virchow2-tcga": 1280}
TILE_CAP = 4096
LR = 2e-4
WD = 1e-5
EPOCHS = 20
BATCH = 32
SEED = 0
TARGET_NAMES = ["keap1_sig", "global", "island", "opensea", "tss200", "body"]

# published Phikon-v2 values the control arms must recover (R15 / m13)
PUBLISHED = {
    "subtype_grouped": 0.799, "subtype_random": 0.9703,
    "keap1_grouped": 0.664, "meth_mean_inflation": 0.0849,
}

app = modal.App("nsclc-rwpr-r18-encoder-compare")
vol = modal.Volume.from_name("nsclc-rwpr-features-tcga", create_if_missing=False)
gcs_secret = modal.Secret.from_name("gcs-key")
image = modal.Image.debian_slim(python_version="3.11").pip_install(
    "numpy==1.26.4", "torch==2.4.1", "scipy==1.14.1", "scikit-learn==1.5.2",
    "google-cloud-storage==2.18.2",
)


def _gcs_client():
    from google.cloud import storage

    key = os.environ["GCP_SA_KEY"]
    p = "/tmp/sa.json"
    open(p, "w").write(key)
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = p
    return storage.Client(project=json.loads(key).get("project_id", "heydonto-425716"))


def build_model(n_out, emb_dim):
    import torch
    import torch.nn as nn

    class GatedMIL(nn.Module):
        """Verbatim from armC/fl_train.py via mil_embed.py; n_out is the only parameter."""

        def __init__(self):
            super().__init__()
            self.proj = nn.Sequential(nn.Linear(emb_dim, 512), nn.ReLU(), nn.Dropout(0.25))
            self.att_v = nn.Linear(512, 256)
            self.att_u = nn.Linear(512, 256)
            self.att_w = nn.Linear(256, 1)
            self.head = nn.Linear(512, n_out)

        def pooled(self, x, mask=None):
            single = x.dim() == 2
            if single:
                x = x.unsqueeze(0)
            h = self.proj(x)
            a = self.att_w(torch.tanh(self.att_v(h)) * torch.sigmoid(self.att_u(h)))
            if mask is not None:
                a = a.masked_fill(~mask.unsqueeze(-1), float("-inf"))
            w = torch.softmax(a, dim=1)
            z = (w * h).sum(dim=1)
            return z.squeeze(0) if single else z

        def forward(self, x, mask=None):
            z = self.pooled(x, mask)
            single = z.dim() == 1
            out = self.head(z.unsqueeze(0) if single else z)
            if n_out == 1:
                out = out.squeeze(-1)
            return out.squeeze(0) if single else out

    return GatedMIL()


def _pad(bags, device, emb_dim):
    import numpy as np
    import torch

    L = max(b.shape[0] for b in bags)
    x = torch.zeros(len(bags), L, emb_dim, dtype=torch.float32, device=device)
    m = torch.zeros(len(bags), L, dtype=torch.bool, device=device)
    for i, b in enumerate(bags):
        n = b.shape[0]
        x[i, :n] = torch.from_numpy(b.astype("float32")).to(device)
        m[i, :n] = True
    return x, m


@app.function(image=image, gpu="L4", volumes={"/vol": vol}, secrets=[gcs_secret],
              timeout=14400, memory=32768)
def run(feat_dir: str, folds_mode: str, mode: str) -> dict:
    import numpy as np
    import torch
    from scipy.stats import spearmanr
    from sklearn.metrics import roc_auc_score

    assert folds_mode in ("grouped", "random")
    assert mode in ("subtype", "meth")

    cli = _gcs_client()
    bkt = cli.bucket(BUCKET)
    spec = json.loads(bkt.blob(IN_GROUPED if folds_mode == "grouped" else IN_RANDOM)
                      .download_as_bytes())
    pats, folds = spec["patients"], spec["folds"]

    TG = None
    if mode == "meth":
        tgt = json.loads(bkt.blob(TGT_BLOB).download_as_bytes())
        TG = {p: [tgt["targets"][k][i] for k in TARGET_NAMES]
              for i, p in enumerate(tgt["pids"])}
        pats = [p for p in pats if p["pid"] in TG]

    vdir = f"/vol/{feat_dir}"
    need = sorted({s for p in pats for s in p["slides"]})
    cache, enc_seen, dims = {}, set(), set()
    for sid in need:
        fp = f"{vdir}/{sid}.npz"
        if not os.path.exists(fp):
            continue
        with np.load(fp, allow_pickle=True) as z:
            e = np.ascontiguousarray(z["embeddings"])
            if "encoder" in z.files:
                enc_seen.add(str(z["encoder"]))
        if dims and e.shape[1] not in dims:
            raise RuntimeError(f"{sid}: embedding dim {e.shape[1]} disagrees with {dims}")
        dims.add(e.shape[1])
        if e.shape[0] > TILE_CAP:
            idx = np.random.default_rng(0).choice(e.shape[0], TILE_CAP, replace=False)
            e = e[np.sort(idx)]
        cache[sid] = e
    if len(dims) != 1:
        raise RuntimeError(f"mixed feature widths in {vdir}: {sorted(dims)}")
    emb_dim = dims.pop()
    exp = EXPECTED_DIM.get(feat_dir)
    if exp is not None and emb_dim != exp:
        raise RuntimeError(f"{feat_dir}: width {emb_dim} != expected {exp}")
    print(f"[{feat_dir}/{folds_mode}/{mode}] cached {len(cache)}/{len(need)} slides, "
          f"{sum(v.nbytes for v in cache.values())/2**30:.1f} GiB, dim={emb_dim}, "
          f"encoders={sorted(enc_seen)}", flush=True)
    if len(enc_seen) > 1:
        raise RuntimeError(f"mixed encoders in {vdir}: {sorted(enc_seen)}")

    pats = [p for p in pats if any(s in cache for s in p["slides"])]
    pid_of = {p["pid"]: p for p in pats}
    pids_all = [p["pid"] for p in pats]
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"  {len(pats)} patients with >=1 cached slide", flush=True)

    n_out = 6 if mode == "meth" else 1
    targets = TARGET_NAMES if mode == "meth" else ["subtype", "keap1"]

    def slide_items(pids, target=None):
        out = []
        for pid in pids:
            p = pid_of.get(pid)
            if p is None:
                continue
            if mode == "meth":
                if pid not in TG:
                    continue
                y = TG[pid]
            else:
                if p[target] is None:
                    continue
                y = float(p[target])
            for s in p["slides"]:
                if s in cache:
                    out.append((s, y, pid))
        return out

    def train(items, seed):
        torch.manual_seed(seed)
        model = build_model(n_out, emb_dim).to(device).train()
        opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
        rng = np.random.default_rng(seed)
        if mode == "meth":
            lossf = torch.nn.MSELoss()
            Yt = np.array([it[1] for it in items], dtype=np.float64)
            mu, sd = Yt.mean(0), Yt.std(0) + 1e-8   # TRAIN-fold statistics only
        else:
            lossf = torch.nn.BCEWithLogitsLoss()
            mu = sd = None
        for ep in range(EPOCHS):
            order = rng.permutation(len(items))
            tot = 0.0
            for i in range(0, len(order), BATCH):
                idx = order[i:i + BATCH]
                bags = [cache[items[j][0]] for j in idx]
                if mode == "meth":
                    y = torch.tensor((np.array([items[j][1] for j in idx]) - mu) / sd,
                                     dtype=torch.float32, device=device)
                else:
                    y = torch.tensor([items[j][1] for j in idx], dtype=torch.float32,
                                     device=device)
                x, m = _pad(bags, device, emb_dim)
                opt.zero_grad()
                loss = lossf(model(x, m), y)
                loss.backward()
                opt.step()
                tot += float(loss) * len(idx)
            if ep % 5 == 0 or ep == EPOCHS - 1:
                print(f"    epoch {ep} loss {tot/max(len(order),1):.4f}", flush=True)
        return model.eval(), mu, sd

    @torch.no_grad()
    def predict(model, mu, sd, pids):
        out = {}
        for pid in pids:
            p = pid_of.get(pid)
            if p is None:
                continue
            vals = []
            for s in p["slides"]:
                if s not in cache:
                    continue
                z = model.pooled(torch.from_numpy(cache[s].astype("float32")).to(device))
                h = model.head(z.unsqueeze(0)).squeeze(0)
                if mode == "meth":
                    vals.append(h.cpu().numpy() * sd + mu)
                else:
                    vals.append(float(torch.sigmoid(h.squeeze())))
            if vals:
                out[pid] = np.mean(vals, axis=0)
        return out

    oof = {t: {} for t in targets}
    for fi, test_pids in enumerate(folds):
        test_pids = [p for p in test_pids if p in pid_of]
        train_pids = [p for p in pids_all if p not in set(test_pids)]
        print(f"  fold {fi}: train {len(train_pids)} test {len(test_pids)}", flush=True)
        if mode == "meth":
            items = slide_items(train_pids)
            model, mu, sd = train(items, SEED + fi)
            pr = predict(model, mu, sd, test_pids)
            for pid, v in pr.items():
                for k, t in enumerate(TARGET_NAMES):
                    oof[t][pid] = float(v[k])
            del model
        else:
            for target in targets:
                items = slide_items(train_pids, target)
                model, mu, sd = train(items, SEED + fi)
                pr = predict(model, mu, sd, test_pids)
                oof[target].update({k: float(v) for k, v in pr.items()})
                del model
        torch.cuda.empty_cache()

    res = {"feat_dir": feat_dir, "folds": folds_mode, "mode": mode,
           "encoder": sorted(enc_seen), "emb_dim": int(emb_dim), "n_patients": len(pats),
           "n_slides_cached": len(cache), "metrics": {}}
    if mode == "subtype":
        for t in targets:
            keep = [p for p in pids_all if p in oof[t] and pid_of[p][t] is not None]
            y = np.array([pid_of[p][t] for p in keep], float)
            s = np.array([oof[t][p] for p in keep], float)
            res["metrics"][t] = {"auroc": float(roc_auc_score(y, s)), "n": len(keep),
                                 "n_pos": int(y.sum())}
    else:
        for k, t in enumerate(TARGET_NAMES):
            keep = [p for p in pids_all if p in oof[t] and p in TG]
            obs = np.array([TG[p][k] for p in keep], float)
            pr = np.array([oof[t][p] for p in keep], float)
            r = spearmanr(obs, pr)
            res["metrics"][t] = {"rho": float(r.statistic), "p": float(r.pvalue),
                                 "n": len(keep)}
        res["metrics"]["_mean_rho"] = float(np.mean(
            [res["metrics"][t]["rho"] for t in TARGET_NAMES]))

    blob = f"{OUT_PREFIX}{feat_dir}__{folds_mode}__{mode}.json"
    bkt.blob(blob).upload_from_string(json.dumps(res, indent=2),
                                     content_type="application/json")
    print(f"wrote gs://{BUCKET}/{blob}", flush=True)
    print(json.dumps(res["metrics"], indent=2), flush=True)
    return res


@app.local_entrypoint()
def main(feat_dir: str, folds: str = "grouped", mode: str = "subtype"):
    r = run.remote(feat_dir, folds, mode)
    print(json.dumps(r["metrics"], indent=2))
