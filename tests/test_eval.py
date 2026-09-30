import numpy as np
import pandas as pd
import pytest

from cellresp.baselines import BASELINES, additive
from cellresp.metrics import evaluate
from cellresp.splits import assign_folds, check_no_leakage
from synthetic import make


def _table():
    gdsc, *_ = make(n_cells=80, n_drugs=20)
    return (gdsc.rename(columns={"SANGER_MODEL_ID": "depmap_id", "LN_IC50": "y"})
                .assign(drug_key=lambda d: d["DRUG_ID"].astype(str))
                [["depmap_id", "drug_key", "y"]].reset_index(drop=True))


@pytest.mark.parametrize("kind", ["random", "unseen_drug", "unseen_cell"])
def test_splits_have_no_leakage(kind):
    df = _table()
    folds = assign_folds(df, kind, 5, 0)
    check_no_leakage(df, folds, kind)
    assert set(folds) == set(range(5))


def test_additive_recovers_additive_truth():
    df = _table()
    folds = assign_folds(df, "random", 5, 0)
    tr, te = df[folds != 0], df[folds == 0].copy()
    te["pred"] = additive(tr, te)
    m = evaluate(te, min_size=5)
    assert m["per_drug_pearson_mean"] > 0.7   # ranks cells within drug


def test_drug_mean_has_undefined_per_drug_correlation():
    df = _table()
    folds = assign_folds(df, "random", 5, 0)
    tr, te = df[folds != 0], df[folds == 0].copy()
    te["pred"] = BASELINES["drug_mean"](tr, te)
    m = evaluate(te, min_size=5)
    assert m["pearson"] > 0.8                          # pooled looks great...
    assert m["per_drug_n_undefined"] == m["per_drug_n_groups"]  # ...but says nothing within drug


def test_unseen_drug_fallback_predicts_everywhere():
    df = _table()
    folds = assign_folds(df, "unseen_drug", 5, 0)
    tr, te = df[folds != 0], df[folds == 0].copy()
    for fn in BASELINES.values():
        assert np.isfinite(fn(tr, te)).all()
