from __future__ import annotations

from pathlib import Path

import yaml


def load_config(path: str | Path) -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    for key in ("raw_dir", "processed_dir", "reports_dir", "results_dir"):
        Path(cfg[key]).mkdir(parents=True, exist_ok=True)
    return cfg


def raw_path(cfg: dict, source: str) -> Path:
    return Path(cfg["raw_dir"]) / cfg["sources"][source]["filename"]
