# Genome + gene table — what they are

These are **not Runjia’s private data**. They are public reference files.

## `NC_000913.2.fasta`

- The complete DNA sequence of *E. coli* K-12 MG1655, NCBI accession **NC_000913.2**.
- Length: **4,639,675 bp**.
- Used to: locate each sgRNA on the chromosome, cut flanking sequence for family A, compute mappability/GC, etc.
- **Where from:** NCBI (or Pasteur `badSeed_public/ecoli-COLI-K12.fsa`, same build).
- **Fetched for this workspace from NCBI** on 2026-10-05 → `Runjia/data/reference/NC_000913.2.fasta`.

## `genes_NC_000913.2.tsv`

- Table of ~4,499 genes on that same genome: `name`, `left`, `right`, `ori` (strand).
- Used to: assign each guide to a gene, grouped CV by gene, within-gene ρ, transcription features.
- **Where from:** Pasteur `badSeed_public/genes.txt` (https://gitlab.pasteur.fr/dbikard/badSeed_public).
- Runjia’s note: use **this** gene table with **this** genome build so coordinates match GapR-seq and each other. Mixing NCBI builds can silently shift every locus.
- **Fetched for this workspace** → `Runjia/data/reference/genes_NC_000913.2.tsv`.

You do **not** need to get these from Runjia.
