# Follow-up on Runjia’s progress report — Jacky

**Audience:** Runjia (for vetting before it goes into the shared write-up).  
**Date:** 5 October 2026.  
**Scope:** The open threads left after `Runjia/PROGRESS_REPORT.md` Parts 5, 10–12
and the work plan in `Jacky/agent-notes/WORKPLAN.md`. Model tuning and crisprHAL
parity are **out of scope** (yours). SHAP / positions 18–20 remain Joshua’s.

This is a statement of what was measured, what closed, and what is still blocked —
written so you can accept, correct, or reject each claim before anything is
merged into the progress report.

**Artefacts:** `Jacky/agent-notes/results/` (JSON/CSV) and runners under
`Jacky/agent-notes/scripts/`. Small pipeline helpers live locally under
`Runjia/sgrna/` (`evaluate.py`, `transfer.py`, `run_ablation.py`) and have
**not** been pushed unless you ask.

---

## 0 — One-page verdict

| open item (from Parts 5 / 10–12) | what we did | verdict for the paper |
|---|---|---|
| Ceiling needs same-enzyme replicates | Inventoried Guo Data S1–S7; no per-guide replicates; used Fig 2b/2c + Spearman–Brown | Ceiling **~0.94–0.97**; SLICER 0.707–0.721 is **~74%** of the 2-rep-mean ceiling — replace Part 10’s “0.90–0.93 bracket / replicates nobody published” |
| Is global ρ mostly gene ranking? | Within-gene OOF ρ + pick-percentile | **No.** Within-gene ρ ≈ global; flanks still help inside genes |
| Flank effect = sequencing artefact? | High-uniqueness stratum | Survives: Δρ **+0.078 → +0.051** |
| Long-range asymmetry = replication fork? | Per-replichore asymmetry | **No sign flip** — candidate closed |
| Long-range asymmetry = RecA repair? | Flank ablation on ΔrecA labels | Absolute gain shrinks; **relative** gain (Δ/baseline) does **not** — repair account not supported; saturating assay caveat stands |
| Orientation test impossible? | Rebuilt features for 8,954 S6 Good intergenic guides; orientation-conditioned asymmetry | **Asymmetry depends on orientation** (matched Δ ≈ 0.047 opposite vs ≈ 0.005 same) — re-opens transcription/orientation account; see `INTERGENIC_ORIENTATION.md` |
| Cross-organism CIs missing | Grouped bootstrap on test set | *E. coli*→*C. rodentium* **0.701 [0.683, 0.704]**; reverse **0.627 [0.617, 0.636]**; `guide_overlap` fixed (was NaN) |
| Cross-kingdom “≈0” | Bootstrap CIs | Both directions **negative with CIs excluding zero** |
| What is `cut.score`? | Provenance + empirical join | Noshay `cut.score` **=** Guo Data S4 Cas9 (ρ ≈ 1.0); absolute Z of Cas9/dCas9 ratio, 2-rep geometric mean |
| Flank vs dCas9 control abundance (2A.1) | Looked for counts in Data S* | **Blocked** — only Quality flag; needs SRA |

---

## 1 — Guo tables: what we actually have (Thread 1.1)

Drive folder: the shared Guo supplementary set (Data S1–S7 + PDF). Local copies
under `Jacky/agent-notes/guo_tables/`. Full column inventory:
`GUO_TABLES_INVENTORY.md`, `results/guo_data_S_inventory.json`.

**Decisive answer:** every activity sheet is a **single averaged score**. There
are **no** per-guide replicate columns and **no** raw dCas9/Cas9 read counts in
the Excel files.

What *is* published and usable without SRA:

- Data S4: HQ scores for Cas9 (44,163), eSpCas9 (45,071), ΔrecA Cas9 (48,112).
- Data S6: full library Log2 scores + Quality (Good = ≥20 reads at **control /
  dCas9** — threshold only, not the counts).
- Supp PDF Fig S12: replicate-vs-replicate read-count agreement (Cas9 / eSp /
  ΔrecA).
- Supp PDF Tables S4–S5: R1/R2 library sizes and mapping ratios (including
  dCas9-R1/R2 as separate libraries — evidence the counts exist in SRA).

**Label join (strip `_Cas9` from our IDs):**

| comparison | n | Spearman |
|---|---:|---:|
| our `cut.score` ↔ Data S4 Cas9 | 13,880 | ≈ 1.000 |
| our `cut.score` ↔ \|Data S6 Cas9\| | 13,880 | ≈ 1.000 |
| Data S4 Cas9 ↔ Data S4 ΔrecA | 13,822 | 0.492 |

So the Noshay label **is** Guo’s Cas9 HQ score. Provenance write-up:
`CUT_SCORE_PROVENANCE.md`. Empirical numbers also in
`results/label_provenance.json` / `guo_label_join.json`.

**Correction to Part 11’s open line** (“a trustworthy ceiling, which needs
replicates nobody has published”): the *per-guide* values were not released, but
the *agreement statistics* were (Fig 2b, Fig S12). That is enough for a
same-enzyme ceiling without SRA. SRA remains optional (Threads 1.3 / 1.6 / 2A.1).

---

## 2 — Ceiling, done properly (Thread 1 / Part 10)

Using the published Fig 2b / 2c R² figures (verify against the main PDF if you
want a tighter citation; arithmetic in `results/guo_reliability_ceiling.json`):

| assumption | ceiling |
|---|---:|
| One replicate (√√0.78) | 0.940 |
| **2-replicate mean (Spearman–Brown)** | **0.968** |
| Independent tiling library (√0.771) | 0.937 |

SLICER tuned **0.721 / 0.968 ≈ 74%** of the way to the mean-label ceiling.

**Suggested replacement text for Part 10’s ceiling subsection:**

- Drop the two-enzyme 0.90–0.93 bracket as the *primary* ceiling (keep
  `ceiling2.py` as the negative result that killed that proxy).
- Quote the same-enzyme figure: ~0.94 (one run) / ~0.97 (published average).
- Caveats to keep: (i) attenuation is derived for Pearson under additive noise —
  √Spearman is an approximation; (ii) shared library composition can inflate
  Fig 2b; (iii) the tiling comparison does not share the library and lands in
  the same place.

Thread **1.6** (train on one replicate vs mean) and a bootstrap CI on the tiling
901-guide set remain nice-to-haves; neither blocks the revised sentence.

---

## 3 — Metric discipline: within-gene ρ (Thread 4 / Part 2)

Label `within_gene_variance_share` = **0.845** (your number, reproduced).

| metric | baseline | +`a_flank` |
|---|---:|---:|
| Global OOF ρ | 0.528 | 0.606 |
| **Within-gene mean ρ** (980 genes ≥10 guides) | **0.510** | **0.571** |
| Pick-percentile (library draw) | 68.7th | — |
| Pick-percentile (**within gene**) | **72.4th** | — |

Δ within-gene from flanks = **+0.061**.

**Verdict for Part 2:** global ρ is *not* mostly gene-ranking. The
practitioner-facing claim (pick a better guide *inside* a gene) survives, and
the flank effect is not an artefact of between-gene structure alone.

Artefact: `results/within_gene_metrics.json`.

---

## 4 — Flank effect controls (Thread 2A / Part 5)

### 4.1 High-uniqueness stratum (2A.2) — runnable, done

| stratum | n | baseline ρ | +flank ρ | Δρ |
|---|---:|---:|---:|---:|
| all guides | 13,880 | 0.529 | 0.607 | **+0.078** |
| high uniqueness (top quartile mean uniq) | 4,949 | 0.490 | 0.540 | **+0.051** |

Survives the mappability control. Smaller than on all guides, but still large
and positive — not the same class of artefact as the overturned supercoiling
result. Artefacts: `results/flank_*.csv`.

### 4.2 Flank vs dCas9 control abundance (2A.1) — blocked

Data S6 only publishes a Quality flag derived from the ≥20-read dCas9 threshold.
Per-guide control-arm abundances are not in Data S1–S7. **Cannot run without
SRA recount** (PRJNA450978 and co-accessions — still verify).

What we *can* already say for the paper without 2A.1: the cut score is a
Cas9/dCas9 **ratio**, so pure library-composition bias largely cancels (defence
against the supercoiling-class artefact). Gene-level knockdown lethality in the
dCas9 arm can still add **between-gene** structure — which is why Thread 4’s
within-gene numbers matter.

---

## 5 — Long-range asymmetry candidates (Thread 2B / Part 5)

### 5.1 Replication-fork / replichore (2B.1) — closed

| replichore | up gain | down gain | asymmetry (matched) |
|---|---:|---:|---:|
| left | +0.023 | +0.069 | **0.049** |
| right | +0.023 | +0.072 | **0.048** |

No sign flip, no material magnitude change. Long-range downstream asymmetry is
**not** explained by replication-fork co-orientation. Artefacts:
`results/asymmetry_replichore*.csv`.

### 5.2 Locus-dependent repair / ΔrecA (2B.2) — run; not supportive of repair

Flank ablation under gene-grouped CV on Guo ΔrecA labels joined to the same
13,880 rows (`results/reca_flank_summary.csv`):

| label | baseline ρ | +flank ρ | Δρ | Δ / baseline |
|---|---:|---:|---:|---:|
| published `cut.score` | 0.529 | 0.607 | +0.078 | 0.148 |
| Guo S4 Cas9 | 0.531 | 0.607 | +0.076 | 0.144 |
| **Guo S4 ΔrecA** | **0.227** | **0.267** | **+0.040** | **0.175** |
| Guo S6 \|Cas9\| | 0.531 | 0.607 | +0.076 | 0.143 |
| Guo S6 \|ΔrecA\| | 0.229 | 0.271 | +0.042 | 0.181 |

**Reading (please check):**

1. Absolute flank gain shrinks in ΔrecA (~half). That is **expected** from the
   saturating assay you and Guo already noted (their ΔrecA model ρ 0.328 vs
   Cas9 0.542; our baseline falls from ~0.53 to ~0.23).
2. Normalising by each label’s own baseline (as the work plan required), the
   **relative** flank contribution does **not** shrink — it is slightly larger
   in ΔrecA (0.175–0.181 vs 0.144–0.148).
3. So this is **not** evidence that the long-range flank effect is primarily
   RecA-dependent repair. Report as: candidate not confirmed; absolute shrink
   explained by compressed dynamic range; relative ratio uninformative for a
   positive repair claim.

### 5.3 Intergenic orientation (2B.3) — features rebuilt; asymmetry depends on orientation

**Update (done):** rebuilt base + `a_flank` for intergenic guides and ran the
orientation-conditioned asymmetry test. Full write-up:
`INTERGENIC_ORIENTATION.md`.

Primary set = Data S6 Quality=Good promoter+RBS (**8,954**). Data S4 HQ
intergenic is almost only RBS (1 promoter / 3,644 RBS) — used as sensitivity.

| set | stratum | n | asymmetry_matched |
|---|---|---:|---:|
| s6_good | all | 8954 | 0.034 |
| s6_good | **opposite_gene** | 4361 | **0.047** |
| s6_good | **same_as_gene** | 4593 | **0.005** |
| s4_hq | opposite_gene | 1866 | **0.064** |
| s4_hq | same_as_gene | 1779 | **−0.015** |

**Verdict for Part 5:** Soften “orientation impossible on this screen” →
impossible on the *coding* matrix. On intergenic guides the PAM-side advantage
**depends on orientation vs the nearest gene** (large when opposite, ~0 when
same). That re-opens a transcription/orientation account, with caveats
(nearest-gene proxy; different label context; single seed).

Artefacts: `results/intergenic_orientation_*.csv`,
`intergenic_orientation_rebuild_audit.json`.

---

## 6 — Transfer intervals + cross-kingdom (Thread 3 / Parts 5 & 10)

### 6.1 Cross-organism bootstrap CIs

Grouped 100 kb bootstrap on the **test** set (fixes Part 12 item 2). Also fixed
`guide_overlap`, which was always NaN (`transfer.py` now carries protospacers in
meta).

| train → test | ρ | 95% CI | guide overlap |
|---|---:|---|---:|
| *E. coli* WT → *C. rodentium* Tev | **0.701** | [0.683, 0.704] | ≈ 0 |
| *C. rodentium* Tev → *E. coli* WT | **0.627** | [0.617, 0.636] | ≈ 0 |

Cross-organism transfer is sequence generalisation, not label leakage.
Artefact: `results/transfer_bootstrap.csv`.

### 6.2 Cross-kingdom (shared 6,216 columns)

| direction | ρ | 95% CI |
|---|---:|---|
| *E. coli* → human | **−0.048** | [−0.062, −0.033] |
| human → *E. coli* | **−0.056** | [−0.073, −0.038] |
| *E. coli* → *E. coli* (80/20) | 0.513 | [0.485, 0.541] |
| human → human (80/20) | 0.380 | [0.350, 0.411] |

Both cross-kingdom CIs **exclude zero** — worse than chance, not just “weak.”
(Your Part 10 quoted −0.017 for human→*E. coli*; our shared-column rerun lands
at −0.056. Worth reconciling which split/seed you used — direction and sign
agree.)

### 6.3 Spacer GC quadratic turning points (CIs new)

| | linear ρ | quadratic R² | turning point | CI | sits at |
|---|---:|---:|---:|---|---|
| *E. coli* | −0.201 | 0.045 | 0.290 | [0.220, 0.334] | **1st** pct |
| human | +0.017 | 0.0042 | **0.571** | [0.555, 0.588] | **60th** pct |

Matches Part 10’s point estimates; turning-point CIs are new.
`results/gc_quadratic_proper.json`.

---

## 7 — What I would change in the progress report (for your edit)

Suggested edits only — do **not** treat as already applied.

1. **Part 10 ceiling.** Replace the two-enzyme bracket as the headline ceiling
   with same-enzyme Fig 2b/2c + Spearman–Brown (~0.94 / ~0.97). Keep
   `ceiling2` as the reason the two-enzyme proxy fails.
2. **Part 11 open.** Strike “replicates nobody has published” as an absolute
   blocker; note agreement stats exist; SRA only needed for per-guide replicate
   training / dCas9 abundance.
3. **Part 5 §3.** Add: (i) high-uniq flank survival +0.051; (ii) replichore
   closed; (iii) ΔrecA relative flank gain does not shrink; (iv) intergenic
   both-strand population exists but outside the 13,880 matrix.
4. **Part 2 / metric.** Add within-gene ρ ≈ 0.51 / 0.57 and within-gene
   pick-percentile 72.4th.
5. **Part 5 / Part 12 transfer table.** Attach bootstrap CIs; note
   `guide_overlap ≈ 0` for cross-organism cells.
6. **Part 10 cross-kingdom.** Prefer CI-backed numbers; reconcile the
   human→*E. coli* point estimate if you want a single quoted figure.
7. **Appendix A / methods.** Point at `CUT_SCORE_PROVENANCE.md` for what
   `cut.score` is.

---

## 8 — Still open / blocked (honest list)

| item | status | unblock |
|---|---|---|
| crisprHAL parity / tuning | **yours** | — |
| 2A.1 flank vs dCas9 abundance | blocked | SRA recount |
| 1.3 / 1.6 per-replicate ceiling demo | blocked / optional | SRA |
| 2B.3 orientation ablation | **done** — asymmetry orientation-dependent | optional: CI on Δ asymmetry |
| Distant-bacterium screen | wet-lab | — |
| Physical mechanism of 250–500 nt flank gradient | open | new hypothesis + test |

---

## 9 — How to reproduce

```bash
export RSI09_ROOT=/workspace/Runjia   # or your checkout of Runjia/
python3 Jacky/agent-notes/scripts/workplan.py --all
python3 Jacky/agent-notes/scripts/run_reca_and_orientation.py
```

Guo Excel files belong in `Jacky/agent-notes/guo_tables/` (Data S1–S7 + PDF).

---

## 10 — Asks for you (vet checklist)

Please mark each as OK / change / reject:

- [ ] Ceiling sentence → ~0.94 / ~0.97, SLICER ~74% of mean-label ceiling  
- [ ] Within-gene ρ table into Part 2  
- [ ] High-uniq flank +0.051 into Part 5 controls  
- [ ] Replichore candidate **closed**  
- [ ] ΔrecA: “not supportive of repair” (with saturation caveat) — wording OK?  
- [ ] Soften “orientation impossible” → coding-matrix only  
- [ ] Transfer CIs + guide_overlap fix into the cross-organism table  
- [ ] Cross-kingdom CIs excluding zero  
- [ ] `cut.score` = Guo S4 Cas9 (ρ≈1) into methods  
- [ ] Leave 2A.1 / SRA as optional follow-up, not a paper blocker  

Happy to rewrite any section into progress-report prose once you’ve marked these up.
