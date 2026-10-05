# Upload checklist — everything the agent needs from you

Put files under **`Runjia/`** so paths match the team layout. Best: one zip of
`Runjia/data/` + `Runjia/external_data/` (+ optional `Runjia/results/`), or
upload files into those folders in this agent workspace.

After upload, tell me and I will run `python -m sgrna.workplan --describe` then
`--all`.

---

## Priority A — required for the first real runs

These unlock: within-gene ρ, flank high-uniqueness stratum, replichore test,
label provenance, Guo ceiling (already done without these).

### A1. Core E. coli matrix + genome

| # | Path under `Runjia/` | What it is | Approx size |
|---|---|---|---|
| 1 | `data/raw/ecoli_feature_matrix.csv` | Noshay et al. 2023 E. coli matrix (`sgRNAID`, `cut.score`, ~6k feature cols). **Must be read with header row 1** (title row on line 0). | ~200 MB |
| 2 | `data/reference/NC_000913.2.fasta` | *E. coli* K-12 MG1655 genome NC_000913.2 | ~5 MB |
| 3 | `data/reference/genes_NC_000913.2.tsv` | Gene table (name, coords, strand/`ori`, etc.) used by `build_guide_index` and grouped CV | small |
| 4 | `data/reference/dinucleotide_shape_scales.tsv` | Dinucleotide shape scales (also exists in `sgrna/features/` — copy into `data/reference/` too) | tiny |

### A2. Prebuilt interim (strongly preferred — saves hours)

If Runjia already built these, **upload them instead of making me rebuild**.

| # | Path under `Runjia/` | What it is |
|---|---|---|
| 5 | `data/interim/guide_index.csv` | Guide → genomic position, gene, strand |
| 6 | `data/interim/base_X.npy` | Cached float32 base matrix |
| 7 | `data/interim/base_meta.json` | IDs, `y`, feature names for the cache |
| 8 | `data/interim/qct_lookup_tables.json` | Recovered QCT tables (optional if matrix present) |
| 9 | `data/interim/features/a_flank.csv` | Family A flank block |
| 10 | `data/interim/features/g_transcription.csv` | Needed if you want expression-banded asymmetry re-runs |

Without 5–9 I can rebuild from A1, but `build_guide_index` + `build_features --families a_flank` + `build_matrix --cache` must succeed first.

### A3. One number from Runjia (text is fine)

| # | Item | How he gets it |
|---|---|---|
| 11 | `within_gene_variance_share` from `results/ceiling.json` | `python -m sgrna.diagnose --what ceiling` on his machine |

Paste the number (or the whole JSON) in chat — no file required.

---

## Priority B — required for cross-organism / transfer bootstrap

| # | Path under `Runjia/` | What it is |
|---|---|---|
| 12 | `external_data/crisprHAL/data/WT-SpCas9_training_data.csv` | Curated WT screen (headerless: context, score) |
| 13 | `external_data/crisprHAL/data/WT-SpCas9_testing_data.csv` | Same |
| 14 | `external_data/crisprHAL/data/eSpCas9_training_data.csv` | eSpCas9 screen (label provenance / ceiling2) |
| 15 | `external_data/crisprHAL/data/eSpCas9_testing_data.csv` | Same |
| 16 | `external_data/crisprHAL/data/TevSpCas9_training_data.csv` | *C. rodentium* screen |
| 17 | `external_data/crisprHAL/data/TevSpCas9_testing_data.csv` | Same |
| 18 | `external_data/C_rodentium/NC_013716.1.fasta` | *C. rodentium* ICC168 genome (TevSpCas9 context is short; flanks need this) |

Optional but useful if present:

| # | Path | Notes |
|---|---|---|
| 19 | `data/interim/expanded_X.npy` + `expanded_meta.json` | Featurised 33k curated guides |
| 20 | `data/interim/features/a_flank_WT-SpCas9.csv` etc. | Per-screen flank caches |
| 21 | `data/interim/transfer_*/` | Existing transfer caches (`X.npy`, `meta.json`) |

---

## Priority C — Guo supplementary tables (Thread 1 + 2B)

Download from [Guo et al. NAR 2018](https://doi.org/10.1093/nar/gky572) supplementary materials (or PubMed 29982721).

Upload as Excel/CSV/TSV under e.g. `Jacky/agent-notes/guo_tables/` (any clear names):

| # | Table | Why we need it |
|---|---|---|
| 22 | **Table S6** — Cas9 activity scores | Label provenance; check for per-replicate columns |
| 23 | **Table S7** — eSpCas9 scores | Same |
| 24 | **Table S8** — Cas9 ΔrecA scores | Repair / ΔrecA flank test (2B.2) |
| 25 | **Table S10** — high-quality filtered sets | Resistant-loci filter question |
| 26 | **Table S1** (and/or S2) — intergenic library | Both-strand orientation test (2B.3) |
| 27 | **Table S5** — tiling library | Fig 2c independent-library check |

Also useful: the **supplementary PDF** if tables are only in the PDF.

**What I will look for in S6–S8:** columns for replicate 1 / replicate 2 / raw counts vs only a single activity score. That decides whether SRA recount is needed.

---

## Priority D — cross-kingdom (mentor / Part 10)

| # | Path under `Runjia/` | What it is |
|---|---|---|
| 28 | `data/raw/human_feature_matrix.csv` | Noshay human matrix (same column style as E. coli) |

---

## Priority E — not needed for the first work-plan pass

Skip these unless you already have them handy. They support other feature families, not the checklist items we are running first.

| Path | Used by |
|---|---|
| `external_data/GapR-seq/wig/*.wig.gz` (GSE152880) | `d_supercoiling` only |
| `external_data/E_coli_analysis/data/` (Hi-C) | `e_nucleoid` only |
| `external_data/RegulonDB/PromoterSet.tsv` | `g_transcription` rebuild |
| `external_data/precise1k/data/annotation/gene_info.csv` | expression band rebuild |
| `external_data/crisproff/` | `b_energy` |
| `external_data/CRISPRBact/badSeed_public` | bad-seed features |
| Full crisprHAL model code (beyond `data/` CSVs) | parity / re-run (Runjia’s job) |

---

## Suggested upload packs

**Minimum pack (do this first):** items **1–4** + **5–9** if available + **11** (the one number).

**Transfer pack:** + **12–18**.

**Guo pack:** + **22–27**.

**Human pack:** + **28**.

---

## How to upload into this agent

Any of:

1. Attach files / a zip in the next chat message to this Cloud Agent.  
2. Drop files into `/workspace/Runjia/data/...` and `/workspace/Runjia/external_data/...` if you have workspace upload.  
3. Put a zip at `/workspace/uploads/` and tell me the path — I will unpack into `Runjia/`.

### Preferred zip layout

```
RSI09_data.zip
├── data/
│   ├── raw/
│   │   ├── ecoli_feature_matrix.csv
│   │   └── human_feature_matrix.csv          # optional
│   ├── reference/
│   │   ├── NC_000913.2.fasta
│   │   ├── genes_NC_000913.2.tsv
│   │   └── dinucleotide_shape_scales.tsv
│   └── interim/                              # optional but huge time-saver
│       ├── guide_index.csv
│       ├── base_X.npy
│       ├── base_meta.json
│       └── features/
│           └── a_flank.csv
├── external_data/
│   ├── crisprHAL/data/
│   │   ├── WT-SpCas9_training_data.csv
│   │   ├── WT-SpCas9_testing_data.csv
│   │   ├── eSpCas9_training_data.csv
│   │   ├── eSpCas9_testing_data.csv
│   │   ├── TevSpCas9_training_data.csv
│   │   └── TevSpCas9_testing_data.csv
│   └── C_rodentium/
│       └── NC_013716.1.fasta
└── guo_tables/                               # optional folder name
    ├── Table_S6...
    └── ...
```

---

## After you upload — what I will do

1. Unpack / place files under `Runjia/`.  
2. Fix the known `i_shape` path quirk if needed (`dinucleotide_shape_scales.tsv`).  
3. `python -m sgrna.workplan --describe` until everything Priority A shows `ok`.  
4. Run `--within-gene`, `--flank-stratum`, `--replichore`, `--bootstrap-transfer`, `--label-provenance`.  
5. Inventory Guo tables and update Thread 1 / 2B next steps.  
6. Write results into `Jacky/agent-notes/results/` and update `STATUS.md`.

---

## Quick tick list (print / mark off)

```
[ ] 1  ecoli_feature_matrix.csv
[ ] 2  NC_000913.2.fasta
[ ] 3  genes_NC_000913.2.tsv
[ ] 4  dinucleotide_shape_scales.tsv
[ ] 5  guide_index.csv              (preferred)
[ ] 6  base_X.npy                   (preferred)
[ ] 7  base_meta.json               (preferred)
[ ] 8  qct_lookup_tables.json       (optional)
[ ] 9  features/a_flank.csv         (preferred)
[ ] 10 features/g_transcription.csv (optional)
[ ] 11 within_gene_variance_share   (paste from Runjia)
[ ] 12-15 WT + eSpCas9 crisprHAL CSVs
[ ] 16-17 TevSpCas9 crisprHAL CSVs
[ ] 18 NC_013716.1.fasta
[ ] 19-21 transfer/expanded caches  (optional)
[ ] 22-27 Guo Tables S1/S2/S5/S6/S7/S8/S10
[ ] 28 human_feature_matrix.csv     (optional)
```
