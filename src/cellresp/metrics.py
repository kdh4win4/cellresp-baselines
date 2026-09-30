"""Metrics.

Pooled Pearson over all pairs is dominated by between-drug variance (drugs
differ in potency by orders of magnitude), so a drug-mean baseline already
scores highly. The informative number is the per-drug correlation: for a given
drug, does the model rank cell lines correctly? That is the question a drug
discovery scientist actually asks. Per-cell correlation is the mirror question.

A constant prediction within a group has undefined correlation; such groups
are counted as nan and excluded, and the count is reported.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def _corr(a, b, fn):
    if len(a) < 3 or np.ptp(a) < 1e-12 or np.ptp(b) < 1e-12:
        return np.nan
    return fn(a, b)[0]


def grouped_corr(df: pd.DataFrame, by: str, min_size: int) -> dict:
    vals = []
    for _, g in df.groupby(by):
        if len(g) >= min_size:
            vals.append(_corr(g["y"].to_numpy(), g["pred"].to_numpy(), stats.pearsonr))
    vals = np.array(vals, dtype=float)
    ok = vals[~np.isnan(vals)]
    return {
        f"per_{by}_pearson_mean": float(ok.mean()) if len(ok) else np.nan,
        f"per_{by}_pearson_median": float(np.median(ok)) if len(ok) else np.nan,
        f"per_{by}_n_groups": int(len(vals)),
        f"per_{by}_n_undefined": int(np.isnan(vals).sum()),
    }


def evaluate(df: pd.DataFrame, min_size: int = 10) -> dict:
    y, p = df["y"].to_numpy(), df["pred"].to_numpy()
    out = {
        "n": int(len(df)),
        "rmse": float(np.sqrt(np.mean((y - p) ** 2))),
        "pearson": float(_corr(y, p, stats.pearsonr)),
        "spearman": float(_corr(y, p, stats.spearmanr)),
    }
    out.update(grouped_corr(df.rename(columns={"drug_key": "drug"}), "drug", min_size))
    out.update(grouped_corr(df.rename(columns={"depmap_id": "cell"}), "cell", min_size))
    return out
