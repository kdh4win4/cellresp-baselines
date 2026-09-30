"""Fetch raw sources (when a URL is configured) and write a checksum manifest.

    python -m cellresp.download --config configs/sources.yaml

Files already present in data/raw are never re-downloaded, so manually
downloaded DepMap files are handled the same way. The manifest records
release label, URL, size, sha256 and timestamp for every file, which is the
minimum needed to say exactly which data a result was computed on.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests

from .config import load_config, raw_path


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for block in r.iter_content(1 << 20):
                fh.write(block)
    tmp.rename(dest)


def build_manifest(cfg: dict, download: bool = True) -> dict:
    entries, missing = {}, []
    for name, src in cfg["sources"].items():
        path = raw_path(cfg, name)
        if not path.exists() and download and src.get("url"):
            print(f"[download] {name} <- {src['url']}")
            try:
                fetch(src["url"], path)
            except Exception as exc:  # keep going; report at the end
                print(f"[download] FAILED {name}: {exc}")
        if not path.exists():
            if not src.get("optional"):
                missing.append(name)
            continue
        entries[name] = {
            "filename": path.name,
            "release": src.get("release"),
            "url": src.get("url"),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "recorded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
    manifest = {"sources": entries, "missing_required": missing}
    out = Path(cfg["raw_dir"]) / "manifest.json"
    out.write_text(json.dumps(manifest, indent=2))
    print(f"[manifest] {len(entries)} files -> {out}")
    if missing:
        print(f"[manifest] MISSING required sources: {missing} "
              f"(place them in {cfg['raw_dir']} and re-run)")
    return manifest


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/sources.yaml")
    ap.add_argument("--no-download", action="store_true")
    args = ap.parse_args()
    build_manifest(load_config(args.config), download=not args.no_download)


if __name__ == "__main__":
    main()
