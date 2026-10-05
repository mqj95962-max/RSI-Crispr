# Work-plan status (agent session)

Last updated: agent run on Cloud VM (no team `data/raw` matrix here).

## Done in code (local workspace)

| item | status |
|---|---|
| 3.1 Bootstrap helper | `Runjia/sgrna/evaluate.py` |
| 3.4 `guide_overlap` + `protospacers` | `Runjia/sgrna/transfer.py` |
| 3.2 / 3.3 bootstrap runner | `Runjia/sgrna/workplan.py --bootstrap-transfer` (needs crisprHAL) |
| 4.1–4.3 within-gene + pick-percentile | `workplan.py --within-gene` + OOF in `run_ablation.py` |
| 2A.2 flank high-uniqueness stratum | `workplan.py --flank-stratum` |
| 2B.1 replichore asymmetry | `workplan.py --replichore` |
| 1.4 ceiling from published R² | `results/guo_reliability_ceiling.json` |
| 2A.3 label provenance doc | `CUT_SCORE_PROVENANCE.md` + `--label-provenance` |
| 1.1 table inventory template | `GUO_TABLES_INVENTORY.md` |

## Needs your data (copy into `Runjia/`)

- [ ] `ecoli_feature_matrix.csv` + reference + `external_data/`  
- [ ] Run `python -m sgrna.workplan --all`  
- [ ] Fill `GUO_TABLES_INVENTORY.md` after downloading Guo S3–S10  

## Needs Runjia (one number)

- [ ] **`within_gene_variance_share`** from `results/ceiling.json`

## Not started here (large / wet-lab)

- [ ] 1.3 SRA recount (only if tables lack replicates)  
- [ ] 2A.1 flank vs dCas9 control abundance (needs recount or Guo counts)  
- [ ] 2B.2 ΔrecA label swap (needs Table S8 aligned to matrix IDs)  
- [ ] 2B.3 intergenic orientation (needs Table S6/S1 join)  
- [ ] 3.3 cross-kingdom bootstrap (needs `human_feature_matrix.csv` + loader hook)

## Results files (after you run)

| file | command |
|---|---|
| `guo_reliability_ceiling.json` | `--guo-ceiling` (already created) |
| `within_gene_metrics.json` | `--within-gene` |
| `flank_high_uniqueness_stratum.csv` | `--flank-stratum` |
| `asymmetry_replichore.csv` | `--replichore` |
| `transfer_bootstrap.csv` | `--bootstrap-transfer` |
| `label_provenance.json` | `--label-provenance` |
