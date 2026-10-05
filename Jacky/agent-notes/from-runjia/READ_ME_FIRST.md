# Reply to Jacky — re: UPLOAD_CHECKLIST.md

**Short version: almost nothing needs to come from me.** 24 of the 28 items are
public downloads, and the links are below. Four are derived caches that only
exist on my machine, and for three of those four it is faster for you to rebuild
than for me to transfer. One is a number, pasted at the bottom.

---

## The one number you asked for (item 11)

From `results/ceiling.json`, already computed:

```json
{
  "n_shared_protospacers": 33567,
  "cross_screen_spearman": 0.80953272564487,
  "implied_ceiling_spearman": 0.8997403656860516,
  "n_of_ours_in_both": 8617,
  "published_vs_wt_spearman": 0.9445582395018755,
  "published_vs_esp_spearman": 0.7417071837026117,
  "genes_with_5plus_guides": 1207,
  "within_gene_variance_share": 0.8451212917093656
}
```

**`within_gene_variance_share = 0.8451`**, over 1,207 genes with ≥5 guides.

---

## Priority A1 — fetch these yourself (items 1–4)

| # | file | where it comes from |
|---|---|---|
| 1 | `ecoli_feature_matrix.csv` (189 MB) | Noshay et al. 2023, *NAR* 51:10147 — **Supplementary Table S2**. https://doi.org/10.1093/nar/gkad736 . Read with `header=1`: row 0 is a title row. |
| 2 | `NC_000913.2.fasta` (4.5 MB) | mine came from `badSeed_public/ecoli-COLI-K12.fsa` (https://gitlab.pasteur.fr/dbikard/badSeed_public), which is the NC_000913.2 build. NCBI also serves it: https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.2 |
| 3 | `genes_NC_000913.2.tsv` (100 KB) | `badSeed_public/genes.txt`, same repo — 4,499 genes on the same build. Columns `name, left, right, ori`. |
| 4 | `dinucleotide_shape_scales.tsv` | **already in the repo**: `Runjia/sgrna/features/dinucleotide_shape_scales.tsv`. Copy it to `data/reference/` — that is the path quirk you noted. |

**Build order matters for 2 and 3.** They must be the *same* assembly as the
GapR-seq wiggles and each other, which is why both come from `badSeed_public`
rather than being mixed from different NCBI downloads. Mixing builds shifts every
coordinate silently.

## Priority B — also public (items 12–18)

| # | file | where |
|---|---|---|
| 12–17 | `WT-SpCas9`, `eSpCas9`, `TevSpCas9` × `_training_data.csv` / `_testing_data.csv` | https://github.com/tbrowne5/crisprHAL , the `data/` folder — all six are there, headerless `context,score`. Clone the whole repo; you will want `models/` later anyway. |
| 18 | `NC_013716.1.fasta` — *C. rodentium* ICC168 | https://www.ncbi.nlm.nih.gov/nuccore/NC_013716.1 |

## Priority C — Guo supplementary (items 22–27)

Guo et al. 2018, *NAR* 46:7052: https://doi.org/10.1093/nar/gky572 — supplementary
tables are on the article page. I do not have these downloaded, so you are not
waiting on me for any of them.

## Priority D — human matrix (item 28)

`human_feature_matrix.csv` (250 MB) is the *H. sapiens* table from the **same
Noshay supplementary** as item 1.

---

## The four items that are genuinely mine (5–10, 19–21)

| # | file | size | my advice |
|---|---|---|---|
| 5 | `guide_index.csv` | 28 MB | **worth transferring.** It is the one artefact with real verification work behind it: 13,879 of 13,880 guides located, PAM checked, 13,877 unique in the genome. Rebuilding needs items 1–3 and about an hour. |
| 9 | `features/a_flank.csv` | 36 MB | **worth transferring.** This is the +0.080 ρ family; rebuilding needs the genome and ~40 min. |
| 10 | `features/g_transcription.csv` | 2.8 MB | **worth transferring** if you want the expression-banded asymmetry re-run — it needs RegulonDB *and* PRECISE-1K, which is the most annoying rebuild in the project. |
| 6, 7 | `base_X.npy` + `base_meta.json` | 330 MB | **rebuild instead.** `python -m sgrna.build_matrix --cache` regenerates it from item 1 in about a minute of I/O. Not worth moving a third of a gigabyte. |
| 8 | `qct_lookup_tables.json` | 28 KB | trivial either way; `qct.learn_lookup_tables` rebuilds it exactly from item 1. |
| 19–21 | `expanded_X.npy` etc. | 790 MB+ | **rebuild**, and only if you actually need the 33k curated arm. |

So the useful pack is **items 5, 9, 10 ≈ 67 MB**, which fits anywhere. Say where
to put it and I will send it.

---

## Two things in your work plan that change *my* report, not yours

These matter more than the file logistics, and I am acting on both.

**1. Your ceiling point is right and mine was too pessimistic.** Part 11 says a
trustworthy ceiling "needs replicates nobody has published". Two corrections from
your note: Guo ran two biological replicates per condition and published the
*agreement* (Fig 2b, R² > 0.78) even though the per-guide values were averaged
away; and the tiling-library comparison (Fig 2c, R² = 0.771, 901 shared guides,
**different library**) is a better instrument than replicates for our purpose,
because it dodges the shared-library objection Part 10 raises against the 0.8095
figure. The Spearman–Brown point is the one I had missed entirely: the label is a
2-replicate geometric mean, so its reliability is `2r/(1+r)` and the ceiling for
predicting *the published score* is higher than for predicting one run. I will
revise Part 10 and credit the correction to you.

**2. The orientation test may be live after all.** `asymmetry.py --audit` found
13,825 of 13,879 guides on one strand and I concluded the test was impossible.
Your reading is better: that is the *coding* sub-library, repurposed from a
CRISPRi design that deliberately binds one strand. If the intergenic sub-library
(10,257 sgRNAs, both strands by design) survives into Table S10, there is a
population to run it on. That makes item 26 (Table S1/S2) the highest-value thing
on your whole list, and it is a download rather than something either of us has
to compute. Worth doing before the ceiling work.

---

## One request back

`UPLOAD_CHECKLIST.md` says you will run `python -m sgrna.workplan --describe`.
**There is no `workplan.py` in `sgrna/`** — 43 modules, and that is not one of
them. If your agent is planning to write it, please keep it in
`Jacky/agent-notes/` rather than inside `Runjia/sgrna/`: every number in
`PROGRESS_REPORT.md` is traceable to that package, so I need its contents to stay
something I can account for. Importing it from your folder is completely fine.
