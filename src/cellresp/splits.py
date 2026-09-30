"""Fold assignment for three generalization questions.

random       : held-out (cell, drug) pairs; both cell and drug seen in training.
               Measures matrix completion, the easiest setting.
unseen_drug  : every test drug is absent from training (new compound).
unseen_cell  : every test cell line is absent from training (new cell context).

Group splits assign whole groups to folds, so leakage is impossible by
construction; check_no_leakage() verifies it anyway.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

GROUP_COL = {"unseen_drug": "drug_key", "unseen_cell": "depmap_id"}


def assign_folds(df: pd.DataFrame, kind: str, n_folds: int = 5, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    if kind == "random":
        folds = np.arange(len(df)) % n_folds
        rng.shuffle(folds)
        return folds
    if kind not in GROUP_COL:
        raise ValueError(f"unknown split kind {kind!r}")
    groups = df[GROUP_COL[kind]].to_numpy()
    uniq = np.unique(groups)
    rng.shuffle(uniq)
    fold_of = {g: i % n_folds for i, g in enumerate(uniq)}
    return np.array([fold_of[g] for g in groups])


def check_no_leakage(df: pd.DataFrame, folds: np.ndarray, kind: str) -> None:
    if kind == "random":
        return
    col = GROUP_COL[kind]
    for f in np.unique(folds):
        train = set(df.loc[folds != f, col])
        test = set(df.loc[folds == f, col])
        if train & test:
            raise AssertionError(f"{kind} fold {f}: {len(train & test)} {col} in train and test")
