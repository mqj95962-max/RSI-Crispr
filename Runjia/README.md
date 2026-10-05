# Runjia

Two things live here, for two different audiences.

| | what it is | who it's for |
|---|---|---|
| [`PROGRESS_REPORT.md`](PROGRESS_REPORT.md) | the whole project written up — what we found, what we withdrew, and what to do next | **start here**; written to be read by someone outside the project |
| [`baseline_ecoli_sgRNA_model_swappable.py`](baseline_ecoli_sgRNA_model_swappable.py) | the shared baseline with a one-line model switch: XGBoost / LightGBM / CatBoost / random forest / Ridge | anyone wanting to try a different learner in Colab without touching the CV code |
| [`sgrna/`](sgrna) | the 43-module analysis pipeline behind every number in the report | reading and reviewing the analysis |

**`PROGRESS_REPORT.md` is a synced copy.** Edit it at the source (Runjia's
working folder), not here, or the next sync will overwrite your changes.

**`sgrna/` will import but not run from a clone.** It needs a multi-gigabyte
data layout of other people's published data that is deliberately not in this
repo — see [`sgrna/README.md`](sgrna/README.md) for the layout and where each
input came from. It is here so the analysis is readable and reviewable, not so
it is one command away from reproducing.

## Removed on 5 October 2026, and where each one went

Five planning-stage files were deleted to keep this folder readable. All of them
are superseded by measured results in `sgrna/`, and `git log --diff-filter=D` or
commit `42a3ba7^` recovers any of them if needed.

| removed | replaced by | result it produced |
|---|---|---|
| `Baseline Model` | the repo root's `USE THIS!! FINAL_baseline_ecoli_sgRNA_model_5seeds.py` | superseded before any result |
| `SHAP` | `sgrna/importance.py` (`shap_importance`), and Joshua's `shap_baseline_analysis.py` | a Colab snippet that depended on notebook variables and ran nowhere else |
| `crisproff_dg_decomposition.py` | `sgrna/features/b_energy.py` | +0.007 ρ — a near-linear function of the dinucleotide content the matrix already has |
| `proposed_extensions_flanking_folding_readcount.py` | `sgrna/features/a_flank.py`, `c_folding.py`, and the label curation in `build_matrix.py` | flanks **+0.080 ρ**, the project's one large effect; folding +0.007 |
| `genomic_coordinate_lookup.py` | `sgrna/build_guide_index.py` | 13,879 of 13,880 guides located and verified against NC_000913.2 |

Reference genomes (`*.gb`, `*.fasta`) are now git-ignored: download them rather
than versioning them.
