# What `cut.score` is (provenance chain)

## Chain in this project

1. **Guo et al. 2018** — pooled Cas9 / eSpCas9 depletion screens, ~70k library.  
2. **Noshay et al. 2023** — publish `ecoli_feature_matrix.csv` with column **`cut.score`** (13,880 rows).  
3. **crisprHAL** — read-count-filtered re-derivation, 33,567 guides (Appendix A in progress report).

## Guo activity score (methods summary)

From the Guo *NAR* / bioRxiv methods (verify against PDF):

1. Co-transform sgRNA library with **Cas9** (selective) or **dCas9** (control).  
2. Extract plasmids at exponential phase; count N20 guides in reads (`GCACN20GTTT` motif).  
3. Normalise read counts for library depth (their equation II).  
4. Drop guides with **<20 reads** in the **plasmid** library.  
5. **Two biological replicates** → **geometric mean** per guide.  
6. Activity score = **absolute value of the Z-score** of the depletion metric.

So the label is:

- A **ratio-style** screen (Cas9 arm vs dCas9 arm) — CRISPRi is already in the **denominator**, not a separate assay.  
- Already **denoised** by replicate averaging before anyone else sees it.  
- **Folded** by `|Z|`, which removes sign and compresses tails.

## Implications for this project

- **Thread 1:** Ceiling from **same-enzyme replicate agreement** (Fig 2b) applies to the **published mean label**, not a single replicate. Spearman–Brown: `2r/(1+r)` raises reliability vs one run.  
- **Thread 2A:** Library-composition artefacts **partially cancel** in Cas9/dCas9 ratios; gene-level knockdown lethality in the dCas9 arm can still add **between-gene** structure.  
- **Thread 4:** Check **`within_gene_variance_share`** (label) vs **within-gene OOF ρ** (model) — see `workplan --within-gene`.

## Empirical check (when data present)

```bash
cd Runjia && export RSI09_ROOT=$PWD
python -m sgrna.workplan --label-provenance
```

Writes `Jacky/agent-notes/results/label_provenance.json` including  
`published_vs_wt_spearman` / `published_vs_esp_spearman` from `diagnose.ceiling()`.
