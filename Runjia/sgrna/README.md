# `sgrna` — the analysis pipeline behind `Runjia/PROGRESS_REPORT.md`

Every number in that report comes from this package. Committed as one unit
because the modules import each other; individual files will not run on their
own.

## It needs a data layout that is not in this repo

The code resolves paths relative to a working folder three levels up from
`config.py`, and expects:

```
<project>/
├── data/raw/            ecoli_feature_matrix.csv, human_feature_matrix.csv
├── data/reference/      NC_000913.2.fasta, gene table, dinucleotide_shape_scales.tsv
├── data/interim/        caches this code writes (several GB)
├── external_data/       crisprHAL, CRISPRoff, GapR-seq, Hi-C, RegulonDB, C. rodentium genome
├── results/             every .csv the report cites
└── src/sgrna/           this package
```

Those inputs total several GB and are other people's published data, so they
are deliberately not committed. `external_data/README.md` in the working folder
records where each came from. **Clone this and it will import but not run** —
it is here so the analysis is readable and reviewable, not so it is one command
away from reproducing.

`features/dinucleotide_shape_scales.tsv` is the one exception: a small
parameter table that belongs with the code. In the working folder it lives at
`data/reference/`, which is where `i_shape.py` reads it from.

## What is where

| module | what it does |
|---|---|
| `config.py` | all paths and model settings in one place, including `SELECTOR` |
| `genome.py`, `io_utils.py`, `qct.py` | genome access, caching, the recovered quantum tables |
| `decode_v_columns.py`, `featurise.py`, `expand.py` | the decoding work — regenerate 6,169 of the 6,232 published columns from a 20-mer |
| `build_guide_index.py`, `build_features.py`, `build_matrix.py` | locate guides, build feature families, assemble the matrix |
| `features/` | the nine feature families, one module each (`a_flank` is the one that worked) |
| `run_ablation.py` | the ablation harness: grouped CV, permutation controls, the selector |
| `evaluate.py`, `importance.py`, `attribution.py` | model comparison and the interpretability measurements |
| `diagnose.py`, `ceiling2.py` | the redundancy / headroom / ceiling diagnostics |
| `representation.py` | how few of the 6,232 published columns the model needs (427, it turns out) |
| `importance_models.py` | seven model families × three importance methods — why SHAP survives a change of model and split gain does not |
| `make_diagnostic_blocks.py` | the blocks that decomposed the "supercoiling" result |
| `headtohead.py`, `crisprhal_rerun.py` | the comparison against crisprHAL 2, including re-running their model |
| `transfer.py` | cross-enzyme, cross-organism and hybrid experiments |
| `tune.py`, `seqnet.py`, `floors.py` | hyperparameter search, the CNN baseline, permutation floors |
| `notebook.py` | thin wrappers so the Colab notebook stays short |

## Two settings worth knowing before changing anything

**`config.SELECTOR`** picks the model that ranks columns before the champion is
fitted. It defaults to `"xgboost"` because that is what the frozen baseline
used, and switching to `"lightgbm"` — which is 3× faster and ranks just as well
— moves the baseline off the published 0.2937 / 0.5278 figure the project
reproduces as a harness check. Fast sweeps can use it; headline numbers should
not.

**Grouped cross-validation is not optional for positional features.** The
screen puts ~20 guides in every gene, so plain `KFold` lets a model score well
by recognising a locus. Anything positional is reported under `arc5` or a
blocked scheme; see `run_ablation.load_groups`.
