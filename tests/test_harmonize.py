import pandas as pd

from cellresp.harmonize.assemble import assemble
from cellresp.harmonize.cells import map_cells, norm_name
from cellresp.harmonize.drugs import build_drug_table
from cellresp.harmonize.genes import standardize_genes
from synthetic import make


def test_norm_name():
    assert norm_name("HEL 92.1.7") == norm_name("hel-92-1-7") == "HEL9217"


def test_cell_mapping_priority_and_fallbacks():
    gdsc, model, _, _ = make()
    cells, rep = map_cells(gdsc, model)
    m = cells.set_index("SANGER_MODEL_ID")["match_method"]
    assert m["SIDM00000"] == "name"          # IDs removed on DepMap side
    assert m["SIDM00039"] == "unmatched"     # model missing entirely
    assert m["SIDM00005"] == "sanger_id"
    assert rep["many_to_one_depmap_ids"] == 0


def test_cosmic_float_vs_int():
    gdsc, model, _, _ = make()
    model = model.assign(SangerModelID=None)   # force COSMIC path
    cells, _ = map_cells(gdsc, model)
    assert (cells.set_index("SANGER_MODEL_ID").loc["SIDM00005", "match_method"] == "cosmic_id")


def test_duplicate_drug_ids_pool_under_one_key():
    gdsc, *_ = make()
    pub = {"Drug0": {"status": "ok", "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-N", "n_hits": 1}}
    drugs, rep = build_drug_table(gdsc, pubchem=pub)
    keys = drugs.set_index("DRUG_ID")["drug_key"]
    assert keys[1000] == keys[2000]
    assert rep["unique_drug_keys"] == rep["gdsc_drug_ids"] - 1


def test_gene_standardization():
    _, _, expr, hgnc = make()
    e = expr.set_index("")
    out, genes, rep = standardize_genes(e, hgnc)
    assert list(out.columns) == ["7157", "1956", "3845", "1017"]
    assert rep["symbol_changed_vs_hgnc"] == 1          # OLDSYM -> CDK2


def test_assemble_grain_and_flags():
    gdsc, model, expr, _ = make()
    gdsc.loc[0, "LN_IC50"] = 5.0                         # > ln(10) -> extrapolated
    cells, _ = map_cells(gdsc, model)
    drugs, _ = build_drug_table(gdsc)
    out, rep = assemble(gdsc, cells, drugs, set(expr[""]))
    assert not out.duplicated(["depmap_id", "drug_key"]).any()
    assert out["ic50_extrapolated"].any()
    assert rep["row_flow"]["final_pairs"] == len(out)
