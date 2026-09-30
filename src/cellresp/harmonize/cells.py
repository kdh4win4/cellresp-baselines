"""Map GDSC cell lines to DepMap ModelIDs.

Priority (first hit wins, each hit records its method):
  1. SANGER_MODEL_ID  == Model.SangerModelID
  2. COSMIC_ID        == Model.COSMICID
  3. normalized name  == normalized Model.CellLineName / StrippedCellLineName
     (only when the normalized name is unique on the DepMap side)

Name matching is the fallback because it is where silent errors happen:
'HEL' vs 'HEL-92-1-7', 'KM12' vs 'KM12-C/KM12-SM', case and punctuation
differences. Every name-based match is therefore flagged for review, and
ID-based matches whose names disagree are reported as conflicts.
"""
from __future__ import annotations

import re

import pandas as pd


def norm_name(s) -> str:
    if pd.isna(s):
        return ""
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def _id_str(s: pd.Series) -> pd.Series:
    """COSMIC IDs arrive as int, float ('906826.0') or str depending on reader."""
    num = pd.to_numeric(s, errors="coerce")
    return num.astype("Int64").astype(str).replace("<NA>", pd.NA)


def map_cells(gdsc: pd.DataFrame, model: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    cells = (gdsc[["SANGER_MODEL_ID", "COSMIC_ID", "CELL_LINE_NAME"]]
             .drop_duplicates().reset_index(drop=True))
    cells["cosmic"] = _id_str(cells["COSMIC_ID"])
    cells["name_norm"] = cells["CELL_LINE_NAME"].map(norm_name)

    m = model.copy()
    m["cosmic"] = _id_str(m["COSMICID"])
    by_sanger = (m.dropna(subset=["SangerModelID"])
                  .drop_duplicates("SangerModelID")
                  .set_index("SangerModelID")["ModelID"])
    by_cosmic = (m.dropna(subset=["cosmic"])
                  .drop_duplicates("cosmic")
                  .set_index("cosmic")["ModelID"])

    names = pd.concat([
        m[["ModelID", "CellLineName"]].rename(columns={"CellLineName": "n"}),
        m[["ModelID", "StrippedCellLineName"]].rename(columns={"StrippedCellLineName": "n"}),
    ])
    names["n"] = names["n"].map(norm_name)
    names = names[names["n"] != ""].drop_duplicates()
    counts = names.groupby("n")["ModelID"].nunique()
    unique_names = names[names["n"].isin(counts[counts == 1].index)]
    by_name = unique_names.drop_duplicates("n").set_index("n")["ModelID"]
    ambiguous_names = set(counts[counts > 1].index)

    depmap_id, method = [], []
    for row in cells.itertuples(index=False):
        if pd.notna(row.SANGER_MODEL_ID) and row.SANGER_MODEL_ID in by_sanger.index:
            depmap_id.append(by_sanger[row.SANGER_MODEL_ID]); method.append("sanger_id")
        elif pd.notna(row.cosmic) and row.cosmic in by_cosmic.index:
            depmap_id.append(by_cosmic[row.cosmic]); method.append("cosmic_id")
        elif row.name_norm in by_name.index:
            depmap_id.append(by_name[row.name_norm]); method.append("name")
        else:
            depmap_id.append(pd.NA)
            method.append("ambiguous_name" if row.name_norm in ambiguous_names else "unmatched")
    cells["depmap_id"] = depmap_id
    cells["match_method"] = method

    # Name disagreement on ID-based matches -> review list
    name_of = m.set_index("ModelID")["CellLineName"].map(norm_name)
    matched = cells["depmap_id"].notna()
    cells["depmap_name_norm"] = cells["depmap_id"].map(name_of)
    cells["name_conflict"] = matched & (cells["name_norm"] != cells["depmap_name_norm"]) \
        & cells["match_method"].isin(["sanger_id", "cosmic_id"])

    # Two GDSC lines landing on one DepMap model is a mapping error, not a merge.
    dup = cells[matched].groupby("depmap_id")["SANGER_MODEL_ID"].nunique()
    many_to_one = dup[dup > 1].index
    cells["many_to_one"] = cells["depmap_id"].isin(many_to_one)

    report = {
        "gdsc_cell_lines": int(len(cells)),
        "by_method": cells["match_method"].value_counts().to_dict(),
        "name_conflicts": int(cells["name_conflict"].sum()),
        "many_to_one_depmap_ids": int(len(many_to_one)),
    }
    return cells.drop(columns=["cosmic"]), report
