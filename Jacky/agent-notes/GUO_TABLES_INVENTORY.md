# Guo et al. 2018 — supplementary tables (Thread 1.1)

Paper: *Improved sgRNA design in bacteria via genome-wide activity profiling*,  
*Nucleic Acids Research* 46(14):7052–7069 (2018), doi [10.1093/nar/gky572](https://doi.org/10.1093/nar/gky572).

Download the supplementary PDF/ZIP from the journal page and fill in the column list below.

| table | stated contents | per-replicate columns? | used for |
|---|---|---|---|
| S1 | in silico intergenic sgRNA library | — | 2B.3 orientation |
| S2 | promoter / RBS entries | — | 2B.3 |
| S3 | negative-control sgRNAs | — | — |
| S4 | mapping ratios per library | — | QC |
| S5 | tiling library (3,451 guides) | — | Fig 2c (901 overlap) |
| **S6** | **Cas9 activity scores** | **check** | WT label / ΔrecA compare |
| S7 | eSpCas9 scores | check | ceiling2 / screens |
| S8 | Cas9 (ΔrecA) scores | check | 2B.2 repair test |
| S9 | genome-editing guideline | — | — |
| S10 | high-quality filtered sets | — | resistant-loci filter |

**Decisive question:** Do S6–S8 include **replicate 1 / replicate 2** (or raw counts), or only the **geometric mean** activity score?

- If **only the mean:** Thread 1.3 (SRA recount) is optional; use Fig 2b/2c statistics + `guo_reliability_ceiling.json`.  
- If **replicates exist:** Thread 1.6 (train on one replicate vs mean) becomes runnable without SRA.

Figshare genome-wide map (linked from PubMed):  
https://figshare.com/s/127cecee6f9ea4e814e2 — not the same as per-guide replicate tables.

**Agent status:** inventory template only — **you** download tables and note columns in this file.
