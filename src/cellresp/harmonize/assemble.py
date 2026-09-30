"""Join responses, cell map and drug table into one modelling table.

Output grain: one row per (depmap_id, drug_key). Replicate measurements
(duplicate DRUG_IDs for the same compound, or repeated screens) are averaged
and the count is kept in n_measurements.

IC50 extrapolation: GDSC fits a curve and reports LN_IC50 even when the IC50
lies outside the tested range. Those values are model extrapolations, not
measurements. They are flagged (not dropped) so evaluation can be run with
and without them.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def assemble(gdsc: pd.DataFrame, cells: pd.DataFrame, drugs: pd.DataFrame,
             expression_ids: set[str] | None = None) -> tuple[pd.DataFrame, dict]:
    df = gdsc.merge(cells[["SANGER_MODEL_ID", "COSMIC_ID", "CELL_LINE_NAME",
                           "depmap_id", "match_method", "many_to_one"]],
                    on=["SANGER_MODEL_ID", "COSMIC_ID", "CELL_LINE_NAME"], how="left")
    df = df.merge(drugs[["DRUG_ID", "drug_key"]], on="DRUG_ID", how="left")
    n0 = len(df)
    steps = {"gdsc_rows": n0}

    df = df[df["depmap_id"].notna()]
    steps["after_cell_mapping"] = len(df)
    df = df[~df["many_to_one"]]
    steps["after_many_to_one_removed"] = len(df)
    if expression_ids is not None:
        df = df[df["depmap_id"].isin(expression_ids)]
        steps["after_require_expression"] = len(df)
    df = df[df["LN_IC50"].notna()]
    steps["after_nonmissing_target"] = len(df)

    df["ic50_extrapolated"] = (df["LN_IC50"] > np.log(df["MAX_CONC"])) | \
                              (df["LN_IC50"] < np.log(df["MIN_CONC"]))

    out = (df.groupby(["depmap_id", "drug_key"], as_index=False)
             .agg(LN_IC50=("LN_IC50", "mean"),
                  AUC=("AUC", "mean"),
                  LN_IC50_sd=("LN_IC50", "std"),
                  n_measurements=("LN_IC50", "size"),
                  drug_ids=("DRUG_ID", lambda s: ";".join(map(str, sorted(set(s))))),
                  ic50_extrapolated=("ic50_extrapolated", "any"),
                  max_conc=("MAX_CONC", "max")))
    steps["final_pairs"] = len(out)

    report = {
        "row_flow": steps,
        "cells": int(out["depmap_id"].nunique()),
        "drugs": int(out["drug_key"].nunique()),
        "pairs_with_replicates": int((out["n_measurements"] > 1).sum()),
        "fraction_extrapolated": round(float(out["ic50_extrapolated"].mean()), 4),
        "matrix_density": round(len(out) / (out["depmap_id"].nunique()
                                            * out["drug_key"].nunique()), 4),
    }
    return out, report
