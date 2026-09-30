"""End-to-end entry points.

    python -m cellresp.pipeline harmonize --config configs/sources.yaml
    python -m cellresp.pipeline evaluate  --config configs/sources.yaml [--wandb]

harmonize writes data/processed/{responses,cells,drugs,genes}.parquet,
expression.parquet and reports/harmonization.json.
evaluate writes results/baselines_by_fold.csv, results/baselines_summary.csv
and optionally logs both (plus the raw-data manifest as an artifact) to W&B.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import io
from .baselines import BASELINES
from .config import load_config, raw_path
from .harmonize.assemble import assemble
from .harmonize.cells import map_cells
from .harmonize.drugs import build_drug_table, resolve_pubchem
from .harmonize.genes import standardize_genes
from .metrics import evaluate as score
from .splits import assign_folds, check_no_leakage


def harmonize(cfg: dict) -> dict:
    proc, rep_dir = Path(cfg["processed_dir"]), Path(cfg["reports_dir"])
    hcfg = cfg["harmonize"]

    gdsc = io.load_gdsc(raw_path(cfg, "gdsc2_response"))
    comp_path = raw_path(cfg, "gdsc_compounds")
    compounds = io.load_compounds(comp_path) if comp_path.exists() else None
    model = io.load_model(raw_path(cfg, "depmap_model"))
    expr = io.load_expression(raw_path(cfg, "depmap_expression"))
    hgnc = io.load_hgnc(raw_path(cfg, "hgnc"))

    cells, cell_rep = map_cells(gdsc, model)

    pubchem = None
    if hcfg.get("use_pubchem"):
        names = sorted(gdsc["DRUG_NAME"].dropna().unique())
        pubchem = resolve_pubchem(names, hcfg["pubchem_cache"])
    drugs, drug_rep = build_drug_table(gdsc, compounds, pubchem)

    expr_std, genes, gene_rep = standardize_genes(expr, hgnc)

    expr_ids = set(expr_std.index) if hcfg.get("require_expression", True) else None
    responses, asm_rep = assemble(gdsc, cells, drugs, expr_ids)

    responses.to_parquet(proc / "responses.parquet", index=False)
    cells.to_parquet(proc / "cells.parquet", index=False)
    drugs.to_parquet(proc / "drugs.parquet", index=False)
    genes.to_parquet(proc / "genes.parquet", index=False)
    expr_std.loc[expr_std.index.isin(responses["depmap_id"])].to_parquet(proc / "expression.parquet")

    # Review lists: the rows a human should look at before trusting the table.
    cells[cells["match_method"].isin(["name", "ambiguous_name", "unmatched"])
          | cells["name_conflict"] | cells["many_to_one"]] \
        .to_csv(rep_dir / "cells_to_review.csv", index=False)
    multi = drugs.groupby("name_norm")["DRUG_ID"].transform("nunique") > 1
    drugs[multi | (drugs["pubchem_status"] != "ok")] \
        .to_csv(rep_dir / "drugs_to_review.csv", index=False)

    report = {"cells": cell_rep, "drugs": drug_rep, "genes": gene_rep, "assembled": asm_rep}
    (rep_dir / "harmonization.json").write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(report, indent=2, default=str))
    return report


def run_baselines(responses: pd.DataFrame, ecfg: dict) -> pd.DataFrame:
    df = responses.rename(columns={ecfg["target"]: "y"})[["depmap_id", "drug_key", "y"]]
    df = df.dropna(subset=["y"]).reset_index(drop=True)
    rows = []
    for kind in ecfg["split_kinds"]:
        folds = assign_folds(df, kind, ecfg["n_folds"], ecfg["seed"])
        check_no_leakage(df, folds, kind)
        for f in range(ecfg["n_folds"]):
            train, test = df[folds != f], df[folds == f].copy()
            for name, fn in BASELINES.items():
                kw = {"shrinkage": ecfg.get("additive_shrinkage", 0.0)} if name == "additive" else {}
                test["pred"] = fn(train, test, **kw)
                m = score(test, ecfg["min_group_size"])
                rows.append({"split": kind, "fold": f, "baseline": name, **m})
    return pd.DataFrame(rows)


def summarize(by_fold: pd.DataFrame) -> pd.DataFrame:
    cols = ["rmse", "pearson", "spearman", "per_drug_pearson_mean", "per_cell_pearson_mean"]
    g = by_fold.groupby(["split", "baseline"])[cols]
    out = g.mean().round(3).astype(str) + " ± " + g.std().round(3).astype(str)
    return out.reset_index()


def evaluate_cmd(cfg: dict, use_wandb: bool, include_extrapolated: bool = True,
                 target: str | None = None) -> pd.DataFrame:
    ecfg = dict(cfg["evaluate"])
    if target:
        ecfg["target"] = target
    responses = pd.read_parquet(Path(cfg["processed_dir"]) / "responses.parquet")
    subset = "all"
    if not include_extrapolated:
        responses = responses[~responses["ic50_extrapolated"]]
        subset = "measured_only"
    tag = f"{ecfg['target']}_{subset}"
    by_fold = run_baselines(responses, ecfg)
    summary = summarize(by_fold)
    res = Path(cfg["results_dir"])
    by_fold.to_csv(res / f"baselines_by_fold_{tag}.csv", index=False)
    summary.to_csv(res / f"baselines_summary_{tag}.csv", index=False)
    print(summary.to_string(index=False))

    if use_wandb:
        import wandb
        run = wandb.init(project="cellresp-baselines", job_type="baselines",
                         config={**ecfg, "subset": subset}, name=tag)
        manifest = Path(cfg["raw_dir"]) / "manifest.json"
        if manifest.exists():
            art = wandb.Artifact("raw-data-manifest", type="dataset")
            art.add_file(str(manifest))
            run.log_artifact(art)
        run.log({"baselines_by_fold": wandb.Table(dataframe=by_fold),
                 "baselines_summary": wandb.Table(dataframe=summary)})
        run.finish()
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["harmonize", "evaluate"])
    ap.add_argument("--config", default="configs/sources.yaml")
    ap.add_argument("--wandb", action="store_true")
    ap.add_argument("--target", choices=["LN_IC50", "AUC"], default=None,
                    help="override evaluate.target from the config")
    ap.add_argument("--measured-only", action="store_true",
                    help="exclude pairs whose IC50 is extrapolated beyond the tested range")
    a = ap.parse_args()
    cfg = load_config(a.config)
    if a.step == "harmonize":
        harmonize(cfg)
    else:
        evaluate_cmd(cfg, a.wandb, include_extrapolated=not a.measured_only, target=a.target)


if __name__ == "__main__":
    main()
