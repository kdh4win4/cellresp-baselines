# cellresp-baselines

Harmonized GDSC2 × DepMap drug-response table and statistical baselines,
evaluated under three generalization settings. This is the data and evaluation
layer for cell-viability prediction; models are meant to be added on top and
judged against these numbers.

## What it does

1. **Acquire + version** raw sources and record release, URL, size and sha256
   in `data/raw/manifest.json` (optionally logged to W&B as an artifact).
2. **Harmonize identifiers on three axes**
   - *Cell lines*: GDSC → DepMap `ModelID` via Sanger Model ID, then COSMIC ID,
     then normalized name (unique names only). Every name-based match, ID/name
     conflict and many-to-one collision goes to `reports/cells_to_review.csv`.
   - *Drugs*: GDSC name → PubChem InChIKey / SMILES (cached). `drug_key` is the
     InChIKey when resolved, so duplicate `DRUG_ID`s for one compound pool
     together. Collisions go to `reports/drugs_to_review.csv`.
   - *Genes*: DepMap `SYMBOL (ENTREZ)` columns re-keyed to Entrez, with current
     HGNC symbol and Ensembl ID; symbols changed since release are counted.
3. **Assemble** one row per (cell, drug), averaging replicates, and flag IC50s
   extrapolated outside the tested concentration range.
4. **Evaluate baselines** — global mean, drug-conditioned mean,
   cell-conditioned mean, additive (drug mean + cell effect) — under
   `random`, `unseen_drug`, and `unseen_cell` 5-fold splits.

Row counts at every filtering step are in `reports/harmonization.json`.

## Run

```bash
pip install -r requirements.txt && pip install -e .

# 1. Put DepMap Model.csv and the protein-coding TPM log1p expression file in
#    data/raw/ (filenames in configs/sources.yaml). GDSC and HGNC download
#    automatically if the URLs are still valid; verify them first.
python -m cellresp.download --config configs/sources.yaml

# 2. Harmonize
python -m cellresp.pipeline harmonize

# 3. Baselines (add --wandb to log; --measured-only drops extrapolated IC50s)
python -m cellresp.pipeline evaluate
python -m cellresp.pipeline evaluate --measured-only

pytest -q    # synthetic-data tests, no downloads needed
```

## Reading the results

- **Pooled Pearson is not the headline.** Drugs differ in potency by orders of
  magnitude, so a drug-mean baseline scores a high pooled correlation while
  saying nothing about which cell lines respond. The headline metric is
  **per-drug Pearson** (ranking cells within a drug).
- Drug-mean predictions are constant within a drug, so their per-drug
  correlation is undefined; this is reported as a count, not hidden.
- Under `unseen_cell`, every statistical baseline predicts the same value for
  all new cells within a drug, so per-drug correlation is undefined for all of
  them. This is the setting where cell features (expression) have to earn
  their place.
- Under `unseen_drug`, the drug term falls back to the global mean; the
  additive baseline's remaining signal is general cell sensitivity.

## Decisions (and why)

| Decision | Choice | Reason |
|---|---|---|
| Target | `LN_IC50` (AUC also kept) | Standard in GDSC work; AUC is less affected by extrapolation |
| Extrapolated IC50 | flagged, both variants evaluated | Dropping changes the task; keeping trains on curve-fit extrapolations |
| Duplicate DRUG_IDs | pooled on InChIKey, IDs kept in `drug_ids` | Same compound; pooling decision is reversible |
| Cell match by name | allowed only when unique, always reviewed | Most mapping errors happen here |
| Gene key | Entrez | Stable across releases; symbols are not |

## Results

_Fill in after the first real run: data release labels, final cell/drug/pair
counts, and the summary table from `results/baselines_summary_*.csv`._

## Layout

```
configs/sources.yaml       sources, releases, evaluation settings
src/cellresp/download.py   fetch + sha256 manifest
src/cellresp/io.py         raw loaders with column validation
src/cellresp/harmonize/    cells, drugs, genes, assemble
src/cellresp/splits.py     random / unseen_drug / unseen_cell folds + leakage check
src/cellresp/baselines.py  statistical baselines with explicit fallbacks
src/cellresp/metrics.py    pooled and per-drug / per-cell metrics
src/cellresp/pipeline.py   CLI entry points
tests/                     synthetic data reproducing the real pitfalls
```
