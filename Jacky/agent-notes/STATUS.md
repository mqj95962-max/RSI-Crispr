# Work-plan status

Last updated after Guo Drive tables + ΔrecA/orientation pass.
Full report for Runjia: `FOLLOWUP_REPORT_FOR_RUNJIA.md`.

## Done

| item | status | artefact |
|---|---|---|
| 4.0 within_gene_variance_share | **0.845** | `results/ceiling_from_runjia.json` |
| 4.1–4.3 within-gene ρ + pick-percentile | **done** | `results/within_gene_metrics.json` |
| 2A.2 flank high-uniqueness stratum | **done** — Δρ +0.078 → +0.051 | `results/flank_*.csv` |
| 2B.1 replichore asymmetry | **done** — no sign flip | `results/asymmetry_replichore*.csv` |
| 2B.2 ΔrecA flank ablation | **done** — relative gain does not shrink | `results/reca_flank_*.csv` |
| 2B.3 intergenic orientation **audit** | **done** — both strands exist; 0 in matrix | `results/intergenic_orientation_audit.json` |
| 2B.3 intergenic orientation **ablation** | **done** — asymmetry depends on orientation | `INTERGENIC_ORIENTATION.md`, `results/intergenic_orientation_*.csv` |
| 3.1 bootstrap helpers | **done** | `Runjia/sgrna/evaluate.py` |
| 3.2 transfer bootstrap CIs | **done** | `results/transfer_bootstrap.csv` |
| 3.3 cross-kingdom + GC CIs | **done** | `results/cross_kingdom_bootstrap.json` |
| 3.4 guide_overlap fix | **done** | `Runjia/sgrna/transfer.py` |
| 1.1 Guo Data S column inventory | **done** — no per-guide replicates | `GUO_TABLES_INVENTORY.md` |
| 1.4 Guo Fig 2b/2c ceiling | **done** | `results/guo_reliability_ceiling.json` |
| 2A.3 label provenance | **done** — cut.score = Guo S4 Cas9 | `CUT_SCORE_PROVENANCE.md` |
| Follow-up report for Runjia | **done** | `FOLLOWUP_REPORT_FOR_RUNJIA.md` |

## Still blocked / optional

| item | why |
|---|---|
| 2A.1 flank vs dCas9 abundance | need control-arm counts from SRA |
| 1.3 / 1.6 SRA or replicate training | per-guide replicates not in Excel |
| (none for 2B.3) | ablation complete; optional: bootstrap CI on Δ asymmetry |

## How to re-run

```bash
export RSI09_ROOT=/workspace/Runjia
python3 Jacky/agent-notes/scripts/workplan.py --all
python3 Jacky/agent-notes/scripts/run_reca_and_orientation.py
```
