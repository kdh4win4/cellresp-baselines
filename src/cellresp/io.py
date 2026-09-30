"""Raw-file loaders. They validate columns and nothing else; all mapping
decisions live in cellresp.harmonize so they can be tested in isolation."""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

GENE_COL = re.compile(r"^(?P<symbol>\S+) \((?P<entrez>\d+)\)$")

GDSC_REQUIRED = [
    "CELL_LINE_NAME", "SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME",
    "MIN_CONC", "MAX_CONC", "LN_IC50", "AUC",
]
MODEL_REQUIRED = ["ModelID"]
MODEL_OPTIONAL = ["SangerModelID", "COSMICID", "RRID", "CellLineName",
                  "StrippedCellLineName", "OncotreeLineage"]


def _read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix in (".xlsx", ".xls"):
        return pd.read_excel(path)
    sep = "\t" if path.suffix in (".tsv", ".txt") else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def _require(df: pd.DataFrame, cols: list[str], what: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{what}: missing columns {missing}; got {list(df.columns)[:15]}...")


def load_gdsc(path) -> pd.DataFrame:
    df = _read_table(path)
    _require(df, GDSC_REQUIRED, "GDSC response")
    if "COSMIC_ID" not in df.columns:
        print("[gdsc] COSMIC_ID absent; matching by Sanger Model ID and name only")
        df["COSMIC_ID"] = pd.NA
    return df


def load_compounds(path) -> pd.DataFrame:
    df = _read_table(path)
    _require(df, ["DRUG_ID", "DRUG_NAME"], "GDSC compounds")
    return df


def load_model(path) -> pd.DataFrame:
    df = _read_table(path)
    _require(df, MODEL_REQUIRED, "DepMap Model")
    for c in MODEL_OPTIONAL:
        if c not in df.columns:
            df[c] = pd.NA
    return df


def load_expression(path) -> pd.DataFrame:
    """Return a ModelID x gene-column frame (float32).

    Handles both layouts seen across DepMap releases: ModelID as the first
    (unnamed) column, or long-format metadata columns (SequencingID, ModelID,
    IsDefaultEntryForModel, ...) followed by gene columns.
    """
    df = pd.read_csv(path, low_memory=False)
    if "ModelID" in df.columns:
        if "IsDefaultEntryForModel" in df.columns:
            flag = df["IsDefaultEntryForModel"].astype(str).str.lower()
            df = df[flag.isin(["yes", "true", "1"])]
        idx = df["ModelID"]
    else:
        idx = df.iloc[:, 0]
    gene_cols = [c for c in df.columns if GENE_COL.match(str(c))]
    if not gene_cols:
        raise ValueError("expression: no columns of the form 'SYMBOL (ENTREZ)'")
    out = df[gene_cols].astype("float32")
    out.index = pd.Index(idx.astype(str), name="ModelID")
    if out.index.duplicated().any():
        n = int(out.index.duplicated().sum())
        print(f"[expression] {n} duplicate ModelID rows; keeping first")
        out = out[~out.index.duplicated()]
    return out


def load_hgnc(path) -> pd.DataFrame | None:
    if path is None or not Path(path).exists():
        return None
    return pd.read_csv(path, sep="\t", low_memory=False,
                       usecols=lambda c: c in {"hgnc_id", "symbol", "entrez_id",
                                               "ensembl_gene_id", "status"})
