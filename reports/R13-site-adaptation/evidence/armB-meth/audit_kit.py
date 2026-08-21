#!/usr/bin/env python3
"""Shared audit harness for the report series. Same contract as R01-R03:

  * config written and hashed BEFORE any outcome is read
  * every gate records expected / observed / PASS-FAIL, and is kept on failure
  * failed runs are never deleted
  * write-once run directory with a full file manifest
"""
import hashlib, json, os, platform, subprocess, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Run:
    def __init__(self, tag):
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                             capture_output=True).stdout.decode().strip() or "nogit"
        self.dir = os.path.join(ROOT, "runs",
                               time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()) + f"-{tag}-{sha}")
        os.makedirs(self.dir, exist_ok=True)
        self.gates, self.logf = [], open(os.path.join(self.dir, "log.jsonl"), "a")

    def log(self, e, **kw):
        r = {"t": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "event": e, **kw}
        self.logf.write(json.dumps(r, default=str) + "\n")
        self.logf.flush()
        print(f"[{r['t']}] {e}: " + " ".join(f"{k}={v}" for k, v in kw.items()), flush=True)

    def gate(self, name, expected, observed, ok, note=""):
        g = {"gate": name, "expected": expected, "observed": observed,
             "result": "PASS" if ok else "FAIL", "note": note,
             "t": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        self.gates.append(g)
        self.log("gate", gate=name, expected=expected, observed=observed, result=g["result"])
        json.dump(self.gates, open(os.path.join(self.dir, "gates.json"), "w"),
                  indent=2, default=str)
        if not ok:
            self.log("HALT", reason=f"gate {name} failed")
            self.finalize()
            sys.exit(2)

    def write(self, n, o):
        p = os.path.join(self.dir, n)
        os.makedirs(os.path.dirname(p), exist_ok=True) if os.path.dirname(n) else None
        with open(p, "w") as fh:
            json.dump(o, fh, indent=2, default=str) if n.endswith(".json") else fh.write(o)
        return p

    def start(self, cfg, inputs):
        """Write config + env + input hashes. Returns the config hash."""
        h = sha256(self.write("config.yaml", json.dumps(cfg, indent=2)))
        self.log("config_written", sha256=h)
        self.write("env.txt", f"python {platform.python_version()}\n{platform.platform()}\n"
                   + subprocess.run([sys.executable, "-m", "pip", "freeze"],
                                    capture_output=True).stdout.decode())
        self.write("inputs.json", {p: {"bytes": os.path.getsize(p), "sha256": sha256(p)}
                                   for p in inputs if os.path.exists(p)})
        return h

    def finalize(self):
        L = []
        for dp, _, fs in os.walk(self.dir):
            for f in sorted(fs):
                if f != "MANIFEST.sha256":
                    fp = os.path.join(dp, f)
                    L.append(f"{sha256(fp)}  {os.path.relpath(fp, self.dir)}")
        open(os.path.join(self.dir, "MANIFEST.sha256"), "w").write("\n".join(sorted(L)) + "\n")


def bh(pvals):
    """Benjamini-Hochberg q-values, order preserved."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        i = n - rank
        prev = min(prev, p[idx] * n / i)
        q[idx] = prev
    return q


def boot_ci(fn, n, seed, B=2000):
    rng = np.random.default_rng(seed)
    v = []
    for _ in range(B):
        try:
            v.append(fn(rng.integers(0, n, n)))
        except Exception:
            pass
    a = np.array(v)
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))
