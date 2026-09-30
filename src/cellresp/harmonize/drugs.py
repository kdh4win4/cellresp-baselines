"""Canonical drug identity for GDSC compounds.

GDSC2 has DRUG_IDs that share a name (same compound, different supplier or
screening site). Whether to pool them is a modelling decision, so we keep
DRUG_ID and add drug_key:
  drug_key = InChIKey            if PubChem resolves the name
           = "NAME:<normalized>" otherwise
Pooling on drug_key collapses duplicates; the report lists every collision so
the choice is visible. Only the first 14 characters of the InChIKey (the
connectivity block) are used for the collision check, which also catches
salt / stereo variants.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/{}/JSON"


def norm_drug_name(s) -> str:
    if pd.isna(s):
        return ""
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def _pubchem_lookup(name: str) -> dict:
    for props in ("InChIKey,SMILES", "InChIKey,IsomericSMILES"):
        try:
            r = requests.get(PUG.format(quote(name, safe=""), props), timeout=30)
        except requests.RequestException as exc:
            return {"status": f"error:{type(exc).__name__}"}
        if r.status_code == 404:
            return {"status": "not_found"}
        if r.ok:
            rows = r.json()["PropertyTable"]["Properties"]
            first = rows[0]
            return {"status": "ok", "cid": first.get("CID"),
                    "inchikey": first.get("InChIKey"),
                    "smiles": first.get("SMILES") or first.get("IsomericSMILES"),
                    "n_hits": len(rows)}
    return {"status": f"http_{r.status_code}"}


def resolve_pubchem(names: list[str], cache_path: str | Path, sleep: float = 0.25) -> dict:
    cache_path = Path(cache_path)
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    todo = [n for n in names if n and n not in cache]
    for i, n in enumerate(todo, 1):
        cache[n] = _pubchem_lookup(n)
        time.sleep(sleep)
        if i % 25 == 0:
            cache_path.write_text(json.dumps(cache, indent=1))
            print(f"[pubchem] {i}/{len(todo)}")
    cache_path.write_text(json.dumps(cache, indent=1))
    return cache


def build_drug_table(gdsc: pd.DataFrame, compounds: pd.DataFrame | None = None,
                     pubchem: dict | None = None) -> tuple[pd.DataFrame, dict]:
    drugs = gdsc[["DRUG_ID", "DRUG_NAME"]].drop_duplicates("DRUG_ID").reset_index(drop=True)
    if compounds is not None:
        extra = [c for c in ("SYNONYMS", "TARGET", "TARGET_PATHWAY") if c in compounds.columns]
        drugs = drugs.merge(compounds[["DRUG_ID", *extra]].drop_duplicates("DRUG_ID"),
                            on="DRUG_ID", how="left")
    drugs["name_norm"] = drugs["DRUG_NAME"].map(norm_drug_name)

    pubchem = pubchem or {}
    info = drugs["DRUG_NAME"].map(lambda n: pubchem.get(n, {}))
    drugs["pubchem_status"] = info.map(lambda d: d.get("status", "not_queried"))
    drugs["cid"] = info.map(lambda d: d.get("cid"))
    drugs["inchikey"] = info.map(lambda d: d.get("inchikey"))
    drugs["smiles"] = info.map(lambda d: d.get("smiles"))
    drugs["pubchem_n_hits"] = info.map(lambda d: d.get("n_hits"))

    drugs["drug_key"] = drugs["inchikey"].where(drugs["inchikey"].notna(),
                                                "NAME:" + drugs["name_norm"])

    same_name = drugs.groupby("name_norm")["DRUG_ID"].nunique()
    ik14 = drugs["inchikey"].dropna().str[:14]
    skel = drugs.loc[ik14.index].assign(ik14=ik14).groupby("ik14")["name_norm"].nunique()

    report = {
        "gdsc_drug_ids": int(len(drugs)),
        "unique_drug_keys": int(drugs["drug_key"].nunique()),
        "pubchem_status": drugs["pubchem_status"].value_counts().to_dict(),
        "names_with_multiple_drug_ids": int((same_name > 1).sum()),
        "connectivity_blocks_shared_by_different_names": int((skel > 1).sum()),
        "ambiguous_pubchem_names": int((drugs["pubchem_n_hits"].fillna(0) > 1).sum()),
    }
    return drugs, report
