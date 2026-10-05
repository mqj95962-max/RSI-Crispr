# RSI09 — derived caches from Runjia, 5 October 2026

Six files. **Everything else on `UPLOAD_CHECKLIST.md` is a public download** —
links are in `READ_ME_FIRST.md`, which also answers item 11 and flags two places
where your work plan corrects the progress report.

Unpack at the root of your working folder so the paths land as shown.

| file | size | what it is | rebuild cost if you lose it |
|---|---:|---|---|
| `data/interim/guide_index.csv` | 28 MB | checklist item 5. Guide → genomic position, gene, strand, PAM, flanks. 13,879 of 13,880 located; every one found with the required NGG in the right place; 13,877 unique in the genome. | items 1–3 + ~1 h (`build_guide_index`) |
| `data/interim/features/a_flank.csv` | 36 MB | item 9. Family A — 348 flanking-DNA columns, the project's one large effect (+0.080 ρ, +0.164 on the curated label). | genome + ~40 min |
| `data/interim/features/g_transcription.csv` | 2.8 MB | item 10. 20 transcription/replication columns, needed for the expression-banded asymmetry re-run. | RegulonDB **and** PRECISE-1K — the worst rebuild in the project |
| `data/interim/qct_lookup_tables.json` | 28 KB | item 8. The recovered quantum-chemical lookup tables: 4 / 16 / 64 / 256 k-mers, no conflicts anywhere. | item 1 + minutes |
| `data/interim/base_meta.json` | 504 KB | item 7. Row IDs, the `cut.score` label, and the 6,232 column names. Useful on its own for checking column names before you rebuild `base_X.npy`. | item 1 + ~1 min |
| `data/interim/v_column_mapping.json` | 600 KB | **not on your list, but include it.** The decoding of the anonymous `V####` columns: 5,853 of 5,887 identified as (position, k-mer) indicators, zero ambiguous, verified against 1.17 M rebuilt values. `featurise`, `importance.block_of` and `representation.column_sets` all read it, and it is the one artefact here that would be genuinely hard to reproduce. | days |

## Deliberately not included

| | why |
|---|---|
| `data/interim/base_X.npy` (330 MB) | `python -m sgrna.build_matrix --cache` regenerates it from checklist item 1 in about a minute of I/O. Not worth moving a third of a gigabyte. |
| `data/interim/expanded_X.npy` (790 MB) | Same, and only needed if you work on the 33k curated arm. |
| `data/raw/*.csv`, genomes, crisprHAL CSVs, Guo tables | Public downloads, and not mine to redistribute. Links in `READ_ME_FIRST.md`. |

## One thing to check first

`i_shape` reads `dinucleotide_shape_scales.tsv` from `data/reference/`, but the
copy in the repo lives at `Runjia/sgrna/features/`. Copy it across before running
anything that touches family I — that is the path quirk you spotted.
