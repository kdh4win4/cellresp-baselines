"""Statistical baselines. Each takes train/test frames with columns
depmap_id, drug_key, y and returns predictions for test.

When a group is absent from training (by design in the unseen splits) each
baseline falls back to the next coarser estimate, so every baseline always
predicts. The fallback is part of the baseline's definition and should be
reported with it.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def global_mean(train: pd.DataFrame, test: pd.DataFrame) -> np.ndarray:
    return np.full(len(test), train["y"].mean())


def drug_mean(train, test):
    mu = train["y"].mean()
    d = train.groupby("drug_key")["y"].mean()
    return test["drug_key"].map(d).fillna(mu).to_numpy()


def cell_mean(train, test):
    mu = train["y"].mean()
    c = train.groupby("depmap_id")["y"].mean()
    return test["depmap_id"].map(c).fillna(mu).to_numpy()


def additive(train, test, shrinkage: float = 0.0):
    """y_hat = drug mean + cell effect, where the cell effect is the cell's mean
    residual after removing drug means. Unlike cell_mean this is not biased by
    which drugs a cell happened to be screened against. `shrinkage` is a
    pseudo-count pulling small-sample cell effects toward 0."""
    mu = train["y"].mean()
    d = train.groupby("drug_key")["y"].mean()
    resid = train["y"] - train["drug_key"].map(d)
    g = resid.groupby(train["depmap_id"])
    c = g.sum() / (g.size() + shrinkage)
    return (test["drug_key"].map(d).fillna(mu)
            + test["depmap_id"].map(c).fillna(0.0)).to_numpy()


BASELINES = {
    "global_mean": global_mean,
    "drug_mean": drug_mean,
    "cell_mean": cell_mean,
    "additive": additive,
}
