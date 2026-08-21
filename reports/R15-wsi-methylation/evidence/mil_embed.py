#!/usr/bin/env python3
"""R15 Modal job — attention-MIL over tile features, grouped by tissue-source site.

Why this exists: the first R15 attempt used the armC mean-pooled patient embeddings and its
positive control failed (subtype AUROC 0.737 under grouped CV, against R12/R13's 0.97-0.98).
Averaging a whole slide dilutes tumour with stroma. This runs the same gated attention-MIL
architecture R12 validated, over the same tile features, so the representation is one whose
capability is already established.

Contract:
  in   gs://.../nsclc-rwpr-study/armD-meth/mil_input.json
       {patients: [{pid, slides: [...], site, subtype, keap1}], folds: [[test_pids], ...]}
  out  gs://.../nsclc-rwpr-study/armD-meth/mil_output.npz
       pids, oof_subtype, oof_keap1, Z (out-of-fold attention-pooled 512-d), fold_of

Every patient's prediction AND embedding come from a model that never saw their site, because
the folds are GroupKFold on tissue-source site. Z is therefore usable downstream without
leakage.

GatedMIL is copied verbatim from armC/fl_train.py so the architecture is identical.
"""
import json
import os

import modal

BUCKET = "heydonto-quantara-lungcdx"
IN_BLOB = "nsclc-rwpr-study/armD-meth/mil_input.json"
OUT_BLOB = "nsclc-rwpr-study/armD-meth/mil_output.npz"
VOL_DIR = "/vol/tcga"
EMB_DIM = 1024
TILE_CAP = 4096
LR = 2e-4
WD = 1e-5
EPOCHS = 20
BATCH = 32
SEED = 0

app = modal.App("nsclc-rwpr-armD-meth")
vol = modal.Volume.from_name("nsclc-rwpr-features-tcga", create_if_missing=False)
gcs_secret = modal.Secret.from_name("gcs-key")
image = modal.Image.debian_slim(python_version="3.11").pip_install(
    "numpy==1.26.4", "torch==2.4.1", "google-cloud-storage==2.18.2",
)


def _gcs_client():
    from google.cloud import storage
    key = os.environ["GCP_SA_KEY"]
    p = "/tmp/sa.json"
    open(p, "w").write(key)
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = p
    return storage.Client(project=json.loads(key).get("project_id", "heydonto-425716"))


def build_model():
    import torch.nn as nn

    class GatedMIL(nn.Module):
        def __init__(self):
            super().__init__()
            self.proj = nn.Sequential(nn.Linear(EMB_DIM, 512), nn.ReLU(), nn.Dropout(0.25))
            self.att_v = nn.Linear(512, 256)
            self.att_u = nn.Linear(512, 256)
            self.att_w = nn.Linear(256, 1)
            self.head = nn.Linear(512, 1)

        def pooled(self, x, mask=None):
            import torch
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
            out = self.head(z.unsqueeze(0) if single else z).squeeze(-1)
            return out.squeeze(0) if single else out

    return GatedMIL()


def _pad(bags, device):
    import torch
    L = max(b.shape[0] for b in bags)
    B = len(bags)
    x = torch.zeros(B, L, EMB_DIM, dtype=torch.float32, device=device)
    m = torch.zeros(B, L, dtype=torch.bool, device=device)
    for i, b in enumerate(bags):
        n = b.shape[0]
        x[i, :n] = torch.from_numpy(b.astype("float32")).to(device)
        m[i, :n] = True
    return x, m


@app.function(image=image, gpu="L4", volumes={"/vol": vol}, secrets=[gcs_secret],
              timeout=10800)
def run():
    import numpy as np
    import torch

    cli = _gcs_client()
    bkt = cli.bucket(BUCKET)
    spec = json.loads(bkt.blob(IN_BLOB).download_as_bytes())
    pats = spec["patients"]
    folds = spec["folds"]
    print(f"{len(pats)} patients, {len(folds)} folds", flush=True)

    need = sorted({s for p in pats for s in p["slides"]})
    cache = {}
    for sid in need:
        fp = f"{VOL_DIR}/{sid}.npz"
        if not os.path.exists(fp):
            continue
        with np.load(fp) as z:
            e = np.ascontiguousarray(z["embeddings"])
        if e.shape[0] > TILE_CAP:
            idx = np.random.default_rng(0).choice(e.shape[0], TILE_CAP, replace=False)
            e = e[np.sort(idx)]
        cache[sid] = e
    print(f"cached {len(cache)}/{len(need)} slides, "
          f"{sum(v.nbytes for v in cache.values())/2**30:.1f} GiB", flush=True)

    pats = [p for p in pats if any(s in cache for s in p["slides"])]
    pid_of = {p["pid"]: p for p in pats}
    print(f"{len(pats)} patients with at least one cached slide", flush=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    def slide_items(pids, target):
        out = []
        for pid in pids:
            p = pid_of.get(pid)
            if p is None or p[target] is None:
                continue
            for s in p["slides"]:
                if s in cache:
                    out.append((s, float(p[target]), pid))
        return out

    def train(items, seed):
        torch.manual_seed(seed)
        model = build_model().to(device).train()
        opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
        lossf = torch.nn.BCEWithLogitsLoss()
        rng = np.random.default_rng(seed)
        for ep in range(EPOCHS):
            order = rng.permutation(len(items))
            tot = 0.0
            for i in range(0, len(order), BATCH):
                idx = order[i:i + BATCH]
                bags = [cache[items[j][0]] for j in idx]
                y = torch.tensor([items[j][1] for j in idx], device=device)
                x, m = _pad(bags, device)
                opt.zero_grad()
                out = model(x, m)
                loss = lossf(out, y)
                loss.backward()
                opt.step()
                tot += float(loss) * len(idx)
            if ep % 5 == 0 or ep == EPOCHS - 1:
                print(f"    epoch {ep} loss {tot/max(len(order),1):.4f}", flush=True)
        return model.eval()

    @torch.no_grad()
    def predict(model, pids):
        prob, emb = {}, {}
        for pid in pids:
            p = pid_of.get(pid)
            if p is None:
                continue
            ps, zs = [], []
            for s in p["slides"]:
                if s not in cache:
                    continue
                x = torch.from_numpy(cache[s].astype("float32")).to(device)
                z = model.pooled(x)
                ps.append(float(torch.sigmoid(model.head(z.unsqueeze(0)).squeeze())))
                zs.append(z.cpu().numpy())
            if ps:
                prob[pid] = float(np.mean(ps))
                emb[pid] = np.mean(zs, axis=0)
        return prob, emb

    pids_all = [p["pid"] for p in pats]
    oof = {"subtype": {}, "keap1": {}}
    Z, fold_of = {}, {}
    for fi, test_pids in enumerate(folds):
        test_pids = [p for p in test_pids if p in pid_of]
        train_pids = [p for p in pids_all if p not in set(test_pids)]
        print(f"fold {fi}: train {len(train_pids)} test {len(test_pids)}", flush=True)
        for target in ("subtype", "keap1"):
            items = slide_items(train_pids, target)
            print(f"  {target}: {len(items)} training slides", flush=True)
            model = train(items, SEED + fi)
            prob, emb = predict(model, test_pids)
            oof[target].update(prob)
            if target == "subtype":          # tumour-focusing pretext for the embedding
                Z.update(emb)
                for pid in prob:
                    fold_of[pid] = fi
            del model
            torch.cuda.empty_cache()

    keep = [p for p in pids_all if p in Z]
    out = {
        "pids": np.array(keep, dtype=object),
        "oof_subtype": np.array([oof["subtype"].get(p, np.nan) for p in keep], np.float32),
        "oof_keap1": np.array([oof["keap1"].get(p, np.nan) for p in keep], np.float32),
        "Z": np.stack([Z[p] for p in keep]).astype(np.float32),
        "fold_of": np.array([fold_of.get(p, -1) for p in keep], np.int32),
    }
    p = "/tmp/mil_output.npz"
    np.savez_compressed(p, **out)
    bkt.blob(OUT_BLOB).upload_from_filename(p)
    print(f"uploaded {OUT_BLOB}: {len(keep)} patients, Z {out['Z'].shape}", flush=True)
    return {"n": len(keep), "z_dim": int(out["Z"].shape[1])}


@app.local_entrypoint()
def main():
    print(run.remote())
