# Runjia upstream commit `e563081`

**Commit:** `Runjia: SLICER's hyperparameter search run to convergence -- +0.0131, not +0.0103`
**On:** `smal-boi/RSI-Crispr` main (2026-10-04)
**Status on your fork:** not on `origin/main` yet — fork is still at `3e96c82`

Only file touched: `Runjia/PROGRESS_REPORT.md` (+30 / −10).

## What changed

He re-ran SLICER's hyperparameter search **to convergence** (up to 300 draws per
fold, stop after 80 consecutive non-improving draws) instead of the earlier
fixed 12-draw search.

| | old (12-draw) | new (converged) |
|---|---:|---:|
| default ρ | 0.7078 | 0.7082 |
| tuned ρ | **0.718** | **0.721** |
| tuning gain | **+0.0103 ± 0.0030** (p = 0.0015) | **+0.0131 ± 0.0016** (p = 0.00005) |
| vs crisprHAL margin (+0.0107) | gain ≈ matched the margin | gain **exceeds** the margin |

Per-fold gains now: +0.0139, +0.0111, +0.0149, +0.0117, +0.0136.
627 configs total; each fold found its winner early, then searched 80 more for
nothing. His point in the report: a short search understated the tuning effect
by about a quarter.

## What this gives us

- Confirms the **parity reading more strongly**: tuning alone moves SLICER by
  more than the whole head-to-head margin against untuned crisprHAL. Quoting a
  "win" is even less defensible; parity is the honest claim.
- Updates the ceiling discussion's reference point: "where ρ 0.707–0.721 sits"
  against the 0.90–0.93 bracket (already reflected in `WORKPLAN.md`).

## What this does *not* give us

This commit is **not** the number I asked Runjia for on the within-gene ρ thread.

Still needed from him (already computed on his machine, not in this commit):

- `within_gene_variance_share` from `results/ceiling.json`
  (`diagnose.py --ceiling`)

Still open from the work plan, untouched by this commit:

- same-enzyme ceiling from Guo's Fig 2b / 2c reliability stats
- bootstrap CIs on the cross-organism table
- within-gene ρ on out-of-fold predictions
- the flank artefact / asymmetry mechanism work
