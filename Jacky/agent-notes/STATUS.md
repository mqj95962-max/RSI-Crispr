# Work-plan status

Last updated after full Priority A–D compute pass on Cloud Agent.

## Done

| item | status | artefact |
|---|---|---|
| 4.0 within_gene_variance_share | **0.845** | `results/ceiling_from_runjia.json` |
| 4.1–4.3 within-gene ρ + pick-percentile | **done** | `results/within_gene_metrics.json` |
| 2A.2 flank high-uniqueness stratum | **done** — Δρ +0.078 → +0.051 | `results/flank_*.csv` |
| 2B.1 replichore asymmetry | **done** — no sign flip | `results/asymmetry_replichore*.csv` |
| 3.1 bootstrap helpers | **done** | `Runjia/sgrna/evaluate.py` |
| 3.2 transfer bootstrap CIs | **done** | `results/transfer_bootstrap.csv` |
| 3.3 cross-kingdom + GC CIs | **done** | `results/cross_kingdom_bootstrap.json` |
| 3.4 guide_overlap fix | **done** | `Runjia/sgrna/transfer.py` |
| 1.4 Guo Fig 2b/2c ceiling | **done** | `results/guo_reliability_ceiling.json` |
| 2A.3 label provenance | **done** | `results/label_provenance.json` + `CUT_SCORE_PROVENANCE.md` |
| 4.4/4.5 reporting notes | **done** | `RESULTS.md`, work plan text |

## Blocked on real Guo 2018 supplementary tables

| item | why |
|---|---|
| 1.1 table column inventory | Drive “Guo” xlsx is the wrong paper |
| 1.3 / 1.6 SRA or replicate ceiling | need per-replicate columns decision |
| 2A.1 flank vs dCas9 abundance | need control-arm counts |
| 2B.2 ΔrecA flank test | need Table S8 |
| 2B.3 intergenic orientation | need Tables S1/S2/S6 |

**Action for Jacky:** download Guo NAR 2018 supplements from
https://doi.org/10.1093/nar/gky572 into `Jacky/agent-notes/guo_tables/`.

## How to re-run

```bash
export RSI09_ROOT=/workspace/Runjia
python3 Jacky/agent-notes/scripts/workplan.py --all
```
