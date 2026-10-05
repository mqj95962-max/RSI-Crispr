# How to run the work-plan checks (Jacky)

All agent-written outputs go to **`Jacky/agent-notes/results/`**.  
Code changes live under **`Runjia/sgrna/`** (bootstrap, transfer fix, `workplan.py`).

## What you must have on disk (Runjia’s tree)

Under **`Runjia/`** (same layout as the team repo):

| path | purpose |
|---|---|
| `data/raw/ecoli_feature_matrix.csv` | Noshay matrix (13,880 guides) |
| `data/raw/human_feature_matrix.csv` | cross-kingdom bootstrap (optional) |
| `data/reference/NC_000913.2.fasta` | genome |
| `data/reference/genes_NC_000913.2.tsv` | gene table |
| `data/interim/guide_index.csv` | from `build_guide_index` |
| `external_data/crisprHAL/data/` | WT / Tev screen CSVs |
| `external_data/...` | GapR, crisprHAL, shape scales, etc. |

If anything is missing, Runjia’s `sgrna/README.md` lists the full set.

## Commands (from `Runjia/`)

```bash
cd Runjia
export RSI09_ROOT="$PWD"

# See what is missing
python -m sgrna.workplan --describe

# No wet-lab data — ceiling from Guo Fig 2b/2c (already ran in agent VM)
python -m sgrna.workplan --guo-ceiling

# Label / crisprHAL agreement (needs external crisprHAL)
python -m sgrna.workplan --label-provenance

# Everything that needs the matrix + genome
python -m sgrna.workplan --all
```

Or step by step:

```bash
python -m sgrna.workplan --within-gene
python -m sgrna.workplan --flank-stratum
python -m sgrna.workplan --replichore
python -m sgrna.workplan --bootstrap-transfer
python -m sgrna.workplan --bootstrap-cross-kingdom
```

## One message still for Runjia

Ask for **`within_gene_variance_share`** from `results/ceiling.json`  
(`python -m sgrna.diagnose --what ceiling` on his machine).  
Our runner computes **within-gene ρ on OOF predictions** once you have data; his number is the label-only bound.

## Guo supplementary tables (Thread 1.1 / 2B.2 / 2B.3)

Download from the [NAR article supplementary](https://academic.oup.com/nar/article/51/19/10147/7279034)  
(or PubMed **29982721**): Tables **S3–S10**.  
Inventory notes: `GUO_TABLES_INVENTORY.md`.  
Per-replicate columns are usually **not** published — only geometric-mean scores.

## What the agent already did without your data

- Wrote **`results/guo_reliability_ceiling.json`** (Fig 2b/2c → ceiling bracket).
- Fixed **`transfer.py`**: `protospacers` in `meta.json`, `guide_overlap` computed.
- Added **`evaluate.py`**: bootstrap + within-gene metrics.
- Added **`workplan.py`** runner and **`run_ablation.py`** optional OOF predictions.

## Do not push?

Keep these edits local or on your fork when you are ready. The agent was told not to push to GitHub.
