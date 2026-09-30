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

Data: GDSC2 fitted dose response `27Oct23` (release 8.5, from Cell Model
Passports); DepMap Public `Model.csv` and
`OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` (DepMap Public 26Q1,
downloaded 2026-09-30). Checksums in `data/raw/manifest.json`.

### Harmonization

| Step | Result |
|---|---|
| GDSC2 rows / cell lines / DRUG_IDs | 242,036 / 969 / 295 |
| Cell lines matched by Sanger Model ID | 967 (1 by unique name, 1 unmatched) |
| ID matches whose names disagree | 29, all correct on review |
| DRUG_IDs → canonical drugs | 295 → 286 (9 names carry >1 DRUG_ID) |
| Drugs resolved in PubChem | 238 of 295 (57 keyed by normalized name) |
| Genes (protein-coding, Entrez-keyed) | 19,215 |
| Cell lines with DepMap expression | 718 of 969 |
| Final (cell, drug) pairs | 177,028 (718 × 286, density 0.86) |
| LN_IC50 above max tested concentration | 71.8% of rows (0.3% below min) |

The 29 name disagreements include `U031`/`UO31` and `786-0`/`786O` (digit
zero vs letter O) and a dozen `NCI-H####`/`H####` prefix differences. Name-only
matching would have silently dropped or mismatched these. The analysis set is
restricted to the 718 cell lines with DepMap expression so that baselines and
future expression-based models are scored on the same cells.

### Baselines (5-fold mean; per-fold sd in `results/`)

| Target | Split | drug_mean pooled r | additive pooled r | additive per-drug r |
|---|---|---|---|---|
| LN_IC50 | random | 0.843 | 0.906 | **0.650** |
| LN_IC50 | unseen_drug | — | 0.335 | 0.651 |
| LN_IC50 | unseen_cell | 0.844 | 0.844 | undefined |
| AUC | random | 0.752 | 0.789 | **0.317** |
| AUC | unseen_drug | — | 0.246 | 0.321 |
| AUC | unseen_cell | 0.753 | 0.753 | undefined |
| LN_IC50, measured only | random | 0.875 | 0.901 | **0.422** |
| LN_IC50, measured only | unseen_drug | — | 0.116 | 0.376 |
| LN_IC50, measured only | unseen_cell | 0.875 | 0.875 | undefined |

### Findings

1. **A model must beat per-drug r = 0.65 (LN_IC50, random split) to show it
   has learned anything beyond general cell sensitivity.** The additive
   baseline uses no expression and no chemistry.
2. **Pooled correlation is dominated by drug potency.** Drug means alone give
   pooled r = 0.84 while saying nothing about which cells respond.
3. **Unseen cell lines are where cell features must earn their place.** Every
   statistical baseline predicts a constant within a drug for a new cell, so
   the effective per-drug baseline there is zero.
4. **Unseen drugs keep within-drug ranking** (per-drug r unchanged at 0.65,
   since cell effects are known) but lose potency calibration (pooled r 0.34,
   large fold-to-fold variance with ~57 drugs per fold).
5. **LN_IC50 likely overstates the baseline.** Per-drug r falls from 0.65 to
   0.42 when extrapolated IC50s are removed, and is 0.32 on AUC. Candidate
   explanations, not yet tested: the joint mixed-effects curve fit pulling
   extrapolated values toward cell-level effects, and AUC saturating near 1
   for resistant cells. Either way, results should be reported on AUC and on
   measured-only IC50 alongside LN_IC50.
6. **The additive correction matters when cells see different drug sets.**
   On measured-only data, cell means give per-drug r = 0.29 versus 0.42 for
   the additive baseline.

Reproduce:

```bash
python -m cellresp.pipeline evaluate --target LN_IC50
python -m cellresp.pipeline evaluate --target AUC
python -m cellresp.pipeline evaluate --target LN_IC50 --measured-only
```

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
