"""Small synthetic GDSC/DepMap-shaped data with known structure and the
pitfalls the real data has: a name-only match, an unmatched line, a duplicate
DRUG_ID for one compound, and extrapolated IC50s."""
import numpy as np
import pandas as pd


def make(n_cells=40, n_drugs=12, seed=0):
    rng = np.random.default_rng(seed)
    cells = [f"C{i:03d}" for i in range(n_cells)]
    drug_eff = rng.normal(0, 2, n_drugs)      # potency differences
    cell_eff = rng.normal(0, 0.7, n_cells)    # general sensitivity
    rows = []
    for ci, c in enumerate(cells):
        for di in range(n_drugs):
            if rng.random() < 0.1:
                continue
            y = drug_eff[di] + cell_eff[ci] + rng.normal(0, 0.3)
            rows.append({"COSMIC_ID": 900000 + ci, "CELL_LINE_NAME": f"line-{ci}",
                         "SANGER_MODEL_ID": f"SIDM{ci:05d}", "DRUG_ID": 1000 + di,
                         "DRUG_NAME": f"Drug{di}", "MIN_CONC": 0.001, "MAX_CONC": 10.0,
                         "LN_IC50": y, "AUC": 0.9})
    gdsc = pd.DataFrame(rows)
    # Same compound under a second DRUG_ID (different supplier)
    dup = gdsc[gdsc["DRUG_ID"] == 1000].copy()
    dup["DRUG_ID"] = 2000
    dup["LN_IC50"] += rng.normal(0, 0.1, len(dup))
    gdsc = pd.concat([gdsc, dup], ignore_index=True)

    model = pd.DataFrame({
        "ModelID": [f"ACH-{i:06d}" for i in range(n_cells)],
        "SangerModelID": [f"SIDM{i:05d}" for i in range(n_cells)],
        "COSMICID": [900000.0 + i for i in range(n_cells)],
        "CellLineName": [f"LINE{i}" for i in range(n_cells)],
        "StrippedCellLineName": [f"LINE{i}" for i in range(n_cells)],
    })
    model.loc[0, ["SangerModelID", "COSMICID"]] = [None, None]  # -> name match
    model = model[model["ModelID"] != f"ACH-{n_cells-1:06d}"]    # -> unmatched

    genes = ["TP53 (7157)", "EGFR (1956)", "KRAS (3845)", "OLDSYM (1017)"]
    expr = pd.DataFrame(rng.normal(size=(len(model), len(genes))), columns=genes)
    expr.insert(0, "", model["ModelID"].to_numpy())
    hgnc = pd.DataFrame({"hgnc_id": ["HGNC:11998", "HGNC:3236", "HGNC:6407", "HGNC:1771"],
                         "symbol": ["TP53", "EGFR", "KRAS", "CDK2"],
                         "entrez_id": [7157, 1956, 3845, 1017],
                         "ensembl_gene_id": ["ENSG00000141510", "ENSG00000146648",
                                             "ENSG00000133703", "ENSG00000123374"],
                         "status": "Approved"})
    return gdsc, model, expr, hgnc
