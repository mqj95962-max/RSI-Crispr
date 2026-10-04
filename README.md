# RSI-Crispr

Predicting sgRNA cutting efficiency (`cut.score`) in *E. coli* CRISPR-Cas9
screens, extending the quantum chemical tensor (QCT) feature set of
[Noshay et al. 2023, NAR 51:10147](https://doi.org/10.1093/nar/gkad736).

## Where to start

**`USE THIS!! FINAL_baseline_ecoli_sgRNA_model_5seeds.py`** is the shared
baseline every experiment is compared against: XGBoost, five seeds (41–45),
with feature selection fit on each training fold only so validation folds never
leak in. Trust the five-seed CV scores it reports; the final all-data model it
trains at the end is for inference, not a validation score.

It is written for Google Colab and reads the feature matrix from Google Drive.
Edit `DATA_PATH` near the top if your copy lives elsewhere.

## Layout

Each team member keeps their experiments in their own folder. Most experiment
files are named after what they add on top of the baseline, e.g.
`XGBoost + regional QCT summary feature set + 300 features`.

| folder | what's in it |
|---|---|
| `Jacky/` | gene-essentiality and features A–D experiments, grid search |
| `Joshua/` | QCT-derived feature experiments (regional, PAM-weighted, profile shape, coupling), seed-control tests, SHAP analysis |
| `Runjia/` | swappable baseline (XGBoost / LightGBM / CatBoost / random forest / Ridge), genomic coordinate lookup, CRISPRoff-style ΔG decomposition, flank/folding/read-count extensions; [`sgrna/`](Runjia/sgrna) — the analysis pipeline behind the report (38 modules; needs a data layout not in this repo, see its README); and [`PROGRESS_REPORT.md`](Runjia/PROGRESS_REPORT.md) — the whole project written up for an outside reader (synced copy; edit it at the source) |

## Not versioned

Model training outputs (`*_ecoli/`, `five_seed_champion_ecoli_*/`, `*.pkl`,
`*.npy`) are ignored — rerun the baseline to regenerate them. See
[`.gitignore`](.gitignore).
