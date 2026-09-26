"""Shared experiment runner: configs, provenance, sharded execution, resume.

Layout of one run (never overwritten; run_id is a UTC timestamp):

    research/results/<experiment_id>/<run_id>/
        manifest.json      config, config hash, git commit/dirty flag, package
                           versions, paper reference, shard status table
        progress.json      done/failed/pending counts, elapsed, ETA (heartbeat)
        shards/<id>.parquet  one atomically written file per completed shard
        failed/<id>.txt    traceback of a failed shard

A shard is a pure function of (config, shard spec): it derives all randomness
from integer seed parts, so re-running a shard reproduces it bit for bit.
`--resume <run_dir>` skips shards already present; `--rerun-failed` retries
only shards recorded as failed.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import importlib
import json
import os
import platform
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

import pandas as pd
import yaml

from research.awqpe import CODE_VERSION, PAPER_REFERENCE

REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS_ROOT = REPO_ROOT / "research" / "results"
CONFIG_ROOT = REPO_ROOT / "research" / "configs"


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    if "experiment_id" not in cfg:
        raise ValueError(f"{path}: config must define experiment_id.")
    return cfg


def config_hash(cfg: dict) -> str:
    canonical = json.dumps(cfg, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]


def git_state() -> dict:
    def run(*args):
        try:
            return subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, timeout=30).stdout.strip()
        except Exception:  # git missing is recorded, not fatal
            return ""

    return {"commit": run("rev-parse", "HEAD"), "branch": run("rev-parse", "--abbrev-ref", "HEAD"), "dirty": bool(run("status", "--porcelain", "--", "research"))}


def environment() -> dict:
    from importlib.metadata import PackageNotFoundError, version

    pkgs = {}
    for name in ("numpy", "scipy", "pandas", "pyarrow", "qiskit", "qiskit-aer", "PyYAML"):
        try:
            pkgs[name] = version(name)
        except PackageNotFoundError:
            pkgs[name] = None
    return {"python": sys.version.split()[0], "platform": platform.platform(), "executable": sys.executable, "packages": pkgs}


def _atomic_write_json(path: Path, obj) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")
    os.replace(tmp, path)


def _atomic_write_parquet(path: Path, df: pd.DataFrame) -> None:
    tmp = path.with_suffix(".parquet.tmp")
    df.to_parquet(tmp, index=False)
    os.replace(tmp, path)


def common_arguments(description: str, default_config: str) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--config", default=str(CONFIG_ROOT / default_config), help="YAML experiment configuration")
    ap.add_argument("--seed", type=int, default=None, help="override master_seed from the config")
    ap.add_argument("--replicate", type=int, default=None, help="override the number of replicates per phase")
    ap.add_argument("--resume", default=None, help="existing run directory to continue")
    ap.add_argument("--rerun-failed", action="store_true", help="with --resume: retry only failed shards")
    ap.add_argument("--output", default=str(RESULTS_ROOT), help="results root directory")
    ap.add_argument("--max-workers", type=int, default=max(1, (os.cpu_count() or 2) - 4))
    ap.add_argument("--pilot", action="store_true", help="use the config's pilot overrides (small scale)")
    return ap


def resolve_config(args) -> dict:
    cfg = load_config(args.config)
    if args.pilot:
        cfg = {**cfg, **cfg.get("pilot_overrides", {}), "is_pilot": True}
    cfg.pop("pilot_overrides", None)
    if args.seed is not None:
        cfg["master_seed"] = args.seed
    if args.replicate is not None:
        cfg["replicates"] = args.replicate
    return cfg


def _worker(module: str, func: str, cfg: dict, spec: dict, shard_path: str) -> tuple[str, float]:
    t0 = time.time()
    fn = getattr(importlib.import_module(module), func)
    df = fn(cfg, spec)
    _atomic_write_parquet(Path(shard_path), df)
    return spec["shard_id"], time.time() - t0


def run_sharded(args, cfg: dict, shard_specs: list[dict], shard_fn, extra_manifest: dict | None = None) -> Path:
    """Execute shard_fn(cfg, spec) -> DataFrame for every spec, with resume."""
    ids = [s["shard_id"] for s in shard_specs]
    if len(set(ids)) != len(ids):
        raise ValueError("shard ids must be unique.")
    if args.resume:
        run_dir = Path(args.resume)
        manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
        if manifest["config_hash"] != config_hash(cfg):
            raise SystemExit("refusing to resume: config differs from the run's manifest.")
    else:
        run_id = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + ("_pilot" if cfg.get("is_pilot") else "")
        run_dir = Path(args.output) / cfg["experiment_id"] / run_id
        run_dir.mkdir(parents=True, exist_ok=False)
        manifest = {
            "experiment_id": cfg["experiment_id"],
            "run_id": run_id,
            "created_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
            "code_version": CODE_VERSION,
            "paper_reference": PAPER_REFERENCE,
            "git": git_state(),
            "environment": environment(),
            "command": " ".join(sys.argv),
            "config": cfg,
            "config_hash": config_hash(cfg),
            "shards": {sid: "pending" for sid in ids},
            **(extra_manifest or {}),
        }
    (run_dir / "shards").mkdir(exist_ok=True)
    (run_dir / "failed").mkdir(exist_ok=True)
    for sid in ids:
        if (run_dir / "shards" / f"{sid}.parquet").exists():
            manifest["shards"][sid] = "done"
    todo = [s for s in shard_specs if manifest["shards"].get(s["shard_id"]) != "done"]
    if args.resume and args.rerun_failed:
        todo = [s for s in todo if manifest["shards"].get(s["shard_id"]) == "failed"]
    for s in todo:  # shards about to run are pending, whatever their previous status
        manifest["shards"][s["shard_id"]] = "pending"
    _atomic_write_json(run_dir / "manifest.json", manifest)

    module, func = shard_fn.__module__, shard_fn.__name__
    t_start, done_now, durations = time.time(), 0, []

    def heartbeat():
        states = list(manifest["shards"].values())
        remaining = states.count("pending")
        mean = sum(durations) / len(durations) if durations else None
        eta = mean * remaining / max(1, args.max_workers) if mean else None
        _atomic_write_json(run_dir / "progress.json", {
            "updated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(), "total": len(states),
            "done": states.count("done"), "failed": states.count("failed"), "pending": remaining,
            "elapsed_s": round(time.time() - t_start, 1), "eta_s": None if eta is None else round(eta, 1)})

    heartbeat()
    print(f"[{cfg['experiment_id']}] run_dir={run_dir} shards: {len(todo)} to run, {len(ids) - len(todo)} already done", flush=True)
    workers, total = args.max_workers, len(todo)
    for attempt in range(3):
        if not todo:
            break
        broken = False
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_worker, module, func, cfg, s, str(run_dir / "shards" / f"{s['shard_id']}.parquet")): s for s in todo}
            for fut in as_completed(futures):
                sid = futures[fut]["shard_id"]
                try:
                    _, dt = fut.result()
                    manifest["shards"][sid] = "done"
                    durations.append(dt)
                    done_now += 1
                except BrokenProcessPool:
                    broken = True  # a worker died (usually out of memory); resubmit below
                except Exception:
                    manifest["shards"][sid] = "failed"
                    (run_dir / "failed" / f"{sid}.txt").write_text(traceback.format_exc(), encoding="utf-8")
                    print(f"  shard {sid} FAILED (see failed/{sid}.txt)", flush=True)
                _atomic_write_json(run_dir / "manifest.json", manifest)
                heartbeat()
                if done_now and done_now % max(1, total // 20) == 0:
                    print(f"  {done_now}/{total} shards done", flush=True)
        todo = [s for s in todo if manifest["shards"].get(s["shard_id"]) == "pending"]
        if not broken:
            break
        workers = max(1, workers // 2)
        print(f"  worker pool broke (likely out of memory); retrying {len(todo)} shards with {workers} workers", flush=True)
    for s in todo:
        manifest["shards"][s["shard_id"]] = "failed"
        (run_dir / "failed" / f"{s['shard_id']}.txt").write_text("worker pool broke repeatedly (out of memory?)", encoding="utf-8")
    _atomic_write_json(run_dir / "manifest.json", manifest)
    heartbeat()
    states = list(manifest["shards"].values())
    print(f"[{cfg['experiment_id']}] finished: {states.count('done')} done, {states.count('failed')} failed -> {run_dir}", flush=True)
    return run_dir


def load_results(run_dir: str | Path) -> pd.DataFrame:
    files = sorted((Path(run_dir) / "shards").glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no shards in {run_dir}")
    return pd.concat((pd.read_parquet(f) for f in files), ignore_index=True)


def latest_run(experiment_id: str, root: str | Path = RESULTS_ROOT, include_pilot: bool = True) -> Path:
    runs = sorted(p for p in (Path(root) / experiment_id).glob("*") if p.is_dir() and (include_pilot or not p.name.endswith("_pilot")))
    if not runs:
        raise FileNotFoundError(f"no runs for {experiment_id}")
    return runs[-1]
