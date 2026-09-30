"""Gene identifiers for DepMap expression.

DepMap columns are 'SYMBOL (ENTREZ)'. Entrez is the stable key, so columns are
re-keyed to Entrez IDs. If an HGNC table is supplied, each gene gets its
current approved symbol and Ensembl ID, and symbols that changed since the
DepMap release are counted (these are the ones that break joins with other
datasets keyed by symbol, e.g. LINCS L1000).
"""
from __future__ import annotations

import pandas as pd

from ..io import GENE_COL


def parse_gene_columns(columns) -> pd.DataFrame:
    rows = []
    for c in columns:
        m = GENE_COL.match(str(c))
        if m:
            rows.append({"column": c, "depmap_symbol": m["symbol"], "entrez": m["entrez"]})
    return pd.DataFrame(rows)


def standardize_genes(expr: pd.DataFrame, hgnc: pd.DataFrame | None = None
                      ) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    genes = parse_gene_columns(expr.columns)
    dup = genes["entrez"].duplicated(keep="first")
    genes["duplicate_entrez"] = dup

    if hgnc is not None:
        h = hgnc.dropna(subset=["entrez_id"]).copy()
        h["entrez"] = pd.to_numeric(h["entrez_id"], errors="coerce").astype("Int64").astype(str)
        h = h.drop_duplicates("entrez").set_index("entrez")
        genes["hgnc_symbol"] = genes["entrez"].map(h["symbol"])
        genes["ensembl_gene_id"] = genes["entrez"].map(h["ensembl_gene_id"])
        genes["hgnc_id"] = genes["entrez"].map(h["hgnc_id"])
        genes["symbol_changed"] = genes["hgnc_symbol"].notna() & \
            (genes["hgnc_symbol"] != genes["depmap_symbol"])
    else:
        genes["hgnc_symbol"] = pd.NA
        genes["symbol_changed"] = False

    keep = genes[~dup]
    out = expr[keep["column"].tolist()].copy()
    out.columns = pd.Index(keep["entrez"].tolist(), name="entrez")

    report = {
        "genes": int(len(genes)),
        "duplicate_entrez_dropped": int(dup.sum()),
        "hgnc_table_used": hgnc is not None,
        "no_hgnc_match": int(genes["hgnc_symbol"].isna().sum()) if hgnc is not None else None,
        "symbol_changed_vs_hgnc": int(genes["symbol_changed"].sum()),
    }
    return out, genes, report
