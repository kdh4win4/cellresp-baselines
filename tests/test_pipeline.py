"""End-to-end on synthetic files written to a temp dir."""
import json

import yaml

from cellresp.config import load_config
from cellresp.download import build_manifest
from cellresp.pipeline import evaluate_cmd, harmonize
from synthetic import make


def test_end_to_end(tmp_path):
    gdsc, model, expr, hgnc = make()
    raw = tmp_path / "raw"; raw.mkdir()
    gdsc.to_csv(raw / "gdsc.csv", index=False)
    model.to_csv(raw / "Model.csv", index=False)
    expr.to_csv(raw / "expr.csv", index=False)
    hgnc.to_csv(raw / "hgnc.txt", sep="\t", index=False)
    cfg = {
        "raw_dir": str(raw), "processed_dir": str(tmp_path / "proc"),
        "reports_dir": str(tmp_path / "rep"), "results_dir": str(tmp_path / "res"),
        "sources": {
            "gdsc2_response": {"filename": "gdsc.csv", "url": None},
            "gdsc_compounds": {"filename": "none.csv", "url": None, "optional": True},
            "depmap_model": {"filename": "Model.csv", "url": None},
            "depmap_expression": {"filename": "expr.csv", "url": None},
            "hgnc": {"filename": "hgnc.txt", "url": None, "optional": True},
        },
        "harmonize": {"use_pubchem": False, "require_expression": True},
        "evaluate": {"target": "LN_IC50", "n_folds": 3, "seed": 0,
                     "split_kinds": ["random", "unseen_drug", "unseen_cell"],
                     "min_group_size": 5, "additive_shrinkage": 0.0},
    }
    p = tmp_path / "cfg.yaml"; p.write_text(yaml.safe_dump(cfg))
    cfg = load_config(p)
    man = build_manifest(cfg, download=False)
    assert man["missing_required"] == []
    rep = harmonize(cfg)
    assert rep["cells"]["by_method"]["unmatched"] == 1
    summary = evaluate_cmd(cfg, use_wandb=False)
    assert len(summary) == 3 * 4
    assert (tmp_path / "rep" / "harmonization.json").exists()
