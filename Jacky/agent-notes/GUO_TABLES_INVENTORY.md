# Guo et al. 2018 — supplementary tables (Thread 1.1) — DONE

Paper: *Improved sgRNA design in bacteria via genome-wide activity profiling*,  
*Nucleic Acids Research* 46(14):7052–7069 (2018), doi [10.1093/nar/gky572](https://doi.org/10.1093/nar/gky572).

Source: Google Drive folder shared by Jacky  
(`1Vd7s0CskVTL-tFUDlmmKjZMZFVrnKDJp`). Local copies in `Jacky/agent-notes/guo_tables/`.

**Naming note.** The Drive files are labelled **Data S1–S7** (Excel) plus
`Supplementary Material_180516.pdf` (figures/tables S1–S12). These are the journal
supplementary *datasets*, not the same numbering as the PDF’s Table S1–S10.
Earlier work-plan text that said “Table S6/S8” maps as:

| work-plan name | Drive file / sheet |
|---|---|
| HQ Cas9 scores (≈ Table S10 filter) | **Data S4** sheet `Cas9` (44,163) |
| eSpCas9 HQ | Data S4 `eSpCas9` (45,071) |
| ΔrecA Cas9 HQ | Data S4 `ΔrecA Cas9` (48,112) |
| Full activity scores (Cas9 / eSp / ΔrecA) | **Data S6** sheets |
| Intergenic designs | Data S1 `promoter` + `RBS` |
| Coding / CRISPRi designs | Data S1 `CRISPRi (gene-targeting)` |

## Column inventory

| file | sheets | columns | per-replicate? | used for |
|---|---|---|---|---|
| **Data S1** | promoter (5,560); RBS (4,699); CRISPRi gene-targeting (55,672) | `sgRNAID`, `nucleotide sequence` | no | 2B.3 orientation audit |
| **Data S2** | promoter/RBS entry + statistics | ID + sequence / sgRNA number | no | library design stats |
| **Data S3** | gene-targeting tiling (3,451); negative-control sgRNA | ID + sequence | no | Fig 2c tiling overlap |
| **Data S4** | Cas9; eSpCas9; ΔrecA Cas9 | `sgRNAID`, `score` | **no** | HQ labels; 2B.2 |
| **Data S5** | negative control (2,000) | ID + sequence | no | NC library |
| **Data S6** | README; Cas9; eSpCas9; recA(−) Cas9 (each 64,061) | `sgRNA`, `gene (promoter)`, `Log2_normalized_change`, `Quality` | **no** | WT / ΔrecA labels; Quality flag only |
| **Data S7** | gene / promoter / RBS best guides | locus, best sgRNA, score | no | design recommendations |
| **Supp PDF** | Figs S1–S12; Tables S1–S5 (strains, primers, library sizes, mapping) | — | aggregate R1/R2 **library sizes**, not per-guide | ceiling evidence; 2A.1 blocked |

## Decisive question (Thread 1.1)

**Do S4–S6 include per-replicate scores or raw counts?**  
**No.** Every activity sheet is a single averaged score (`score` or
`Log2_normalized_change`). Data S6’s README states Quality = Good/Bad from a
≥20-read threshold at the **control (dCas9) condition**, but the control counts
themselves are not published.

**Consequence:** Thread **1.3** (SRA recount) is optional for a paper-grade
ceiling; Threads **1.6** (train on one replicate) and **2A.1** (flank vs dCas9
abundance) **cannot** run without SRA. Use Fig 2b/2c + Fig S12 statistics +
`guo_reliability_ceiling.json` instead.

## Label join to Noshay 13,880

Strip `_Cas9` from our `sgRNAID` → Guo `gid`.

| comparison | n overlap | Spearman |
|---|---:|---:|
| our `cut.score` vs Data S4 Cas9 | 13,880 | ≈ **1.000** |
| our `cut.score` vs \|Data S6 Cas9\| | 13,880 | ≈ **1.000** |
| Data S4 Cas9 vs Data S4 ΔrecA | 13,822 | **0.492** |

So Noshay’s published label **is** Guo’s Cas9 HQ score (absolute / Z-score form
as already traced in `CUT_SCORE_PROVENANCE.md`).

## Coding vs intergenic in our matrix

Sequence overlap: **13,879 / 13,880** of the Noshay matrix match Data S1
`CRISPRi (gene-targeting)`. Intergenic promoter/RBS guides exist in Data S1/S4/S6
with both strands, but **are not in the 13,880 matrix**.

## PDF evidence still usable without per-guide replicates

- **Fig S12:** “Read counts of each sgRNA of two biological replicates agree well”
  (Cas9 / eSpCas9 / ΔrecA panels).
- **Tables S4–S5:** R1/R2 library raw sizes and mapping ratios (dCas9-R1/R2 listed
  as separate libraries — counts exist in SRA, not in Excel).
- Main-text **Fig 2b** R² > 0.78 (replicate agreement) and **Fig 2c** R² = 0.771
  (tiling library) remain the ceiling inputs.

Machine-readable inventory: `results/guo_data_S_inventory.json`.
