# Work-plan status

Last updated after Runjia’s 10 Oct reply (`e4a71fc`: orientation accepted;
ceiling → ≈0.94 tiling lead). Reply notes: `RUNJIA_REPLY_2026-10-10.md`.

## Process agreement — `Runjia/sgrna/`

**Any edit under `Runjia/sgrna/` must be called out** in:

1. the PR body (files + one-line intent), **and**
2. the table below.

Runjia’s sync is a straight file copy into his source tree; silent additive
changes on `main` risk being overwritten. Additive / backwards-compatible
edits are welcome — they just need to be visible.

### Pipeline edits touching `Runjia/sgrna/` (log)

| when | files | intent | noted in PR? |
|---|---|---|---|
| follow-up / intergenic work (`7b8266e` era) | `evaluate.py`, `run_ablation.py`, `transfer.py` | bootstrap/within-gene helpers; optional `return_oof`; `protospacers` + real `guide_overlap` (was NaN) | **missed in follow-up text** (“not pushed”) — Runjia found them on main and adopted them. Do not repeat. |

## Done

| item | status | artefact |
|---|---|---|
| 4.0 within_gene_variance_share | **0.845** | `results/ceiling_from_runjia.json` |
| 4.1–4.3 within-gene ρ + pick-percentile | **done** | `results/within_gene_metrics.json` |
| 2A.2 flank high-uniqueness stratum | **done** — Δρ +0.078 → +0.051 | `results/flank_*.csv` |
| 2B.1 replichore asymmetry | **done** — no sign flip | `results/asymmetry_replichore*.csv` |
| 2B.2 ΔrecA flank ablation | **done** — relative gain does not shrink | `results/reca_flank_*.csv` |
| 2B.3 intergenic orientation | **done** — **accepted into Part 5** by Runjia | `INTERGENIC_ORIENTATION.md` |
| 3.1–3.4 transfer CIs + guide_overlap | **done** (helpers adopted by Runjia) | `results/transfer_bootstrap.csv` |
| 1.1 Guo Data S inventory | **done** | `GUO_TABLES_INVENTORY.md` |
| 1.4 ceiling | **settled with Runjia** — working **≈0.94** (Fig 2c); range 0.94–0.97; do **not** lead with 0.968 | `RUNJIA_REPLY_2026-10-10.md` |
| 2A.3 label provenance | **done** — used in his ratio argument against 0.968 | `CUT_SCORE_PROVENANCE.md` |
| Follow-up report | **done**; he replied | `FOLLOWUP_REPORT_FOR_RUNJIA.md` |

## Still blocked / optional

| item | why |
|---|---|
| 2A.1 flank vs dCas9 abundance | need control-arm counts from SRA |
| 1.3 / 1.6 SRA or replicate training | per-guide replicates not in Excel |
| 2B.3 bootstrap CI on Δ asymmetry | optional deepen |
| Better orientation label (TSS/operon) | optional deepen |

## How to re-run

```bash
export RSI09_ROOT=/workspace/Runjia
python3 Jacky/agent-notes/scripts/workplan.py --all
python3 Jacky/agent-notes/scripts/run_reca_and_orientation.py
python3 Jacky/agent-notes/scripts/run_intergenic_orientation.py
```
