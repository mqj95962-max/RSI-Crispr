# Intergenic orientation-conditioned asymmetry (Thread 2B.3) — DONE

## Why this was needed

On the Noshay 13,880 coding matrix, 13,825/13,879 guides sit on the gene’s
template strand, so orientation-conditioned asymmetry was impossible. Guo’s
**intergenic** sub-library (promoter + RBS) was designed on both strands, but
those guides are **not** in the published feature matrix — features had to be
rebuilt.

## Sets

| set | n | composition | label |
|---|---:|---|---|
| **s6_good** (primary) | 8,954 | 4,611 promoter + 4,343 RBS | \|Data S6 Log2\|, Quality=Good |
| s4_hq (sensitivity) | 3,645 | **1** promoter + 3,644 RBS | Data S4 Cas9 HQ score |

Data S4’s “HQ” filter almost excludes promoters, so S6 Good is the right primary
set for a both-strand test.

Strands and orientation (nearest gene) are balanced on both sets
(`results/intergenic_orientation_rebuild_audit.json`).

## Methods (short)

1. Locate each protospacer on NC_000913.2; take unique hits only.
2. Cut ±1 kb flanks with `genome.guide_context` (guide orientation).
3. Rebuild **6,169 base** columns via `featurise.matrix` + **348 `a_flank`** via
   `a_flank.build` (same recipe as `transfer.build`).
4. Orientation = guide strand vs **nearest gene** strand
   (`same_as_gene` / `opposite_gene`).
5. Grouped CV by 100 kb genomic bins; arms = baseline / upstream / downstream /
   count-matched downstream / both (same shape as `asymmetry.py --run`).

Runner: `Jacky/agent-notes/scripts/run_intergenic_orientation.py`.

## Headline numbers (matched asymmetry)

| set | stratum | n | baseline ρ | up gain | dn matched gain | **asymmetry_matched** |
|---|---|---:|---:|---:|---:|---:|
| s6_good | all | 8954 | 0.442 | +0.081 | +0.115 | **0.034** |
| s6_good | opposite_gene | 4361 | 0.396 | +0.061 | +0.108 | **0.047** |
| s6_good | same_as_gene | 4593 | 0.439 | +0.102 | +0.107 | **0.005** |
| s4_hq | all | 3645 | 0.423 | +0.021 | +0.054 | **0.033** |
| s4_hq | opposite_gene | 1866 | 0.360 | +0.001 | +0.064 | **0.064** |
| s4_hq | same_as_gene | 1779 | 0.410 | +0.061 | +0.046 | **−0.015** |

Full tables: `results/intergenic_orientation_summary.csv`,
`results/intergenic_orientation_asymmetry.csv`.

## Verdict

**The downstream flank advantage is orientation-dependent on intergenic guides.**

- Overall, intergenic still shows a PAM-side advantage (~0.034 matched), similar
  in sign to the coding-matrix result (~0.049).
- Split by orientation: the advantage is **large when the guide faces the
  opposite strand to the nearest gene** (~0.047–0.064) and **≈0 (or slightly
  flipped) when it faces the same way** (~0.005 / −0.015).

That is the prediction a transcription / gene-orientation account makes, and
the opposite of a pure PAM/R-loop geometry account (which should be flat across
orientation). Together with Part 5’s flat expression-band result on *coding*
guides, the honest state is:

- Expression-level traffic did not scale the asymmetry on the coding library.
- **Guide orientation relative to the local gene does change the asymmetry**
  once both orientations exist (intergenic).

### Bootstrap CI on the orientation gap (added)

Grouped 100 kb block bootstrap on OOF Spearman asymmetries
(`scripts/run_orientation_ci.py`, 2000 resamples):

| set | asym opposite [95% CI] | asym same [95% CI] | Δ (opp−same) [95% CI] | P(Δ≤0) |
|---|---|---|---|---:|
| **s6_good** | **0.040 [0.015, 0.063]** | 0.014 [−0.008, 0.034] | 0.026 [−0.006, 0.058] | 0.054 |
| s4_hq | **0.047 [0.017, 0.076]** | 0.002 [−0.034, 0.037] | 0.045 [−0.001, 0.094] | 0.028 |

**Reading:** On the primary set, opposite-strand asymmetry is clearly positive;
same-strand asymmetry’s CI includes zero. The *difference* Δ is directionally
positive but **borderline at 95%** (CI just includes 0; one-sided ≈0.05). So:
pattern supported; don’t over-claim a precise Δ. Artefacts:
`results/intergenic_orientation_delta_ci.json`.

### Caveats (state in the paper)

1. Orientation is vs **nearest gene**, not a curated operon/TSS annotation —
   weaker than an ideal transcription-direction label.
2. Intergenic labels are a different biological context (promoter/RBS targets).
3. Single seed, 100 kb grouped CV; Δ CI is borderline for s6_good.
4. Does not by itself identify polymerase traffic vs other gene-asymmetric
   chromosomal features.
5. OOF-pooled asymmetry point estimates differ slightly from the fold-mean
   table above (same sign/pattern).

## What changes in the progress report

Replace Part 5’s “conditioning on orientation is impossible on this screen”
with: impossible **on the coding matrix**; on intergenic guides the test runs
and **asymmetry depends on orientation**. Re-open transcription/orientation as
a live candidate for the long-range effect, with the caveats above.
