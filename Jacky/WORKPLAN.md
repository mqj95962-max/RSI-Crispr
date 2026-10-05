# Jacky — work plan

Scope: the three threads left open after `Runjia/PROGRESS_REPORT.md`, plus a metric
item. **Model tuning and the crisprHAL parity result are Runjia's and are not in
here.** Positions 18–20 SHAP is Joshua's — `Joshua/shap_baseline_analysis.py`
already does the hard part.

Every number quoted from the progress report is cited by part. Numbers taken from
the source papers are marked **verify** where I have read them from a text
rendering rather than the published table.

---

## Three things that change the earlier plan

Worth reading before the threads, because two of them cancel work I previously
suggested and one re-opens a question the report closed.

**1. Read-count filtering is already done, and it is not where the headroom is.**
The curated crisprHAL label (33,567 guides) *is* the read-count-filtered
re-derivation of the same Guo screen (Appendix A), and Part 10 closes the label
as a source of remaining headroom: "not the label being dirty (cleaning reveals
the flank effect rather than raising the ceiling)". So wiring
`filter_by_control_reads` and sweeping a threshold reproduces a decision already
taken upstream. **Do not spend time on it.** The live question is the *ceiling*,
which is a different quantity: not "are the labels clean" but "how repeatable is
the assay".

**2. The replicates exist, and the source paper publishes reliability numbers.**
Guo et al. 2018 (*Improved sgRNA design in bacteria via genome-wide activity
profiling*, NAR; bioRxiv 272377) ran **two biological replicates for each of five
conditions** (Cas9, eSpCas9, Cas9 ΔrecA, dCas9, eSpdCas9 — ten libraries), and
averaged the two replicates as a geometric mean before publishing. It also reports:

| comparison | figure | statistic | what it is |
|---|---|---|---|
| replicate vs replicate | Fig 2b | **R² > 0.78** | same library, same enzyme, independent transformations |
| genome-wide vs an independent tiling library, 901 shared guides | Fig 2c | **R² = 0.771** | **different library**, independent experiment |
| screen vs individual colony counting, 15 guides | Fig 2d | R² = 0.840 | orthogonal assay |

(All **verify** against the published figures.) Part 11 lists "a trustworthy
ceiling, which needs replicates nobody has published" as open. That is too
pessimistic — the replicates were run and averaged away, and the *tiling-library*
comparison is better than replicates for our purpose because it does not share
the library, which is the objection Part 10 raises against the 0.8095 figure
("the shared library pulls the estimate the other way again").

**3. The orientation test may not be impossible after all.** `asymmetry.py --audit`
correctly finds no orientation variation — 13,825 of 13,879 guides on one strand.
But that is the *coding* sub-library, which was repurposed from a CRISPRi library
deliberately designed to bind one strand only. Guo's paper says the **intergenic
sub-library (10,257 sgRNAs targeting 3,142 promoters and 4,174 RBSs) was designed
to target either strand**, explicitly in contrast to the earlier library. If those
guides survived into Guo's high-quality Cas9 set (44,163 guides, Table S10), the
transcription-orientation test has a population to run on. It would be run on
intergenic guides with orientation defined relative to the adjacent operon, not
relative to a host gene, so it is a weaker test than the one originally
designed — but it is not an impossible one.

---

## Thread 1 — The ceiling, done properly

### What it entails

The report's ceiling is a bracket, **0.90–0.93**, from the attenuation argument
applied to two screens agreeing at ρ 0.8095 (Part 10). `ceiling2.py` already
showed the assumption behind it is false: the disagreement between WT-SpCas9 and
eSpCas9 is itself predictable at ρ 0.546 against a control of 0.708, so it is
systematic enzyme biology, not independent noise. That is a good negative result
and it leaves us without a ceiling.

This thread replaces the two-enzyme proxy with a same-enzyme repeatability
estimate, which is what the attenuation argument actually requires.

### Why it matters

It sets the denominator for every claim in the paper. At a ceiling of 0.90, ρ
0.707 is 79% of the way there and the honest summary is "near the limit of the
assay". If the replicate-based ceiling is nearer 0.94–0.97 — which is what the
arithmetic below suggests — then there is materially more headroom, "the features
are exhausted" (Part 10) becomes a statement about *these* features rather than
about the problem, and the paper's framing changes. It is also the one open item
other than parity that can move a headline sentence.

### Steps

- [ ] **1.1 Get the supplementary tables.** Download Guo et al. Tables S3–S10.
      Record, per table, exactly which columns exist. The decisive question is
      whether any table carries **per-replicate** scores or read counts, or only
      the geometric-mean activity score. Cheap, and it determines whether 1.3 is
      needed at all.
- [ ] **1.2 Read the reliability figures off the paper properly.** Confirm the Fig
      2b and 2c statistics: what quantity is correlated (raw counts, normalised
      counts, or final activity score), on how many guides, and whether R² is a
      squared Pearson. This matters: a correlation computed on read counts across
      the whole library — including thousands of inactive guides and a huge
      dynamic range — will overstate the repeatability of the filtered activity
      score we actually train on. Report it as an upper bound until checked.
- [ ] **1.3 If per-replicate numbers are not published, recount from raw reads.**
      The authors describe their pipeline precisely enough to reimplement: merge
      paired-end reads with FLASH, find the `GCACN20GTTT` 28-mer and its reverse
      complement, discard reads with mutations in either 4 bp flank, map the N20
      back to the in-silico library, normalise for library depth, drop guides with
      <20 reads in the plasmid library. Ten libraries. Verify the SRA accession
      first — `ceiling2.py` cites PRJNA450978 "and others" and I could not confirm
      it, so treat the accession as unchecked.
- [ ] **1.4 Compute the ceiling, and say which correlation it is in.** Note a
      methodological nit in the current derivation: the attenuation correction is
      derived for Pearson correlation under additive noise, and Part 10 applies √
      to a *Spearman* 0.8095. Either do it in Pearson and convert, or state the
      approximation. Then:
      - reliability of a single replicate = replicate-replicate correlation;
      - reliability of the published 2-replicate mean, by Spearman–Brown,
        `2r / (1 + r)`;
      - ceiling = √reliability, separately for "predicting one replicate" and
        "predicting the published average", because they are different targets.

      Illustrative arithmetic on the published figures, to be redone on real
      numbers: r = √0.78 = 0.883 gives a single-replicate ceiling of 0.940 and a
      2-replicate-mean ceiling of √(2·0.883/1.883) = 0.969. The independent tiling
      comparison gives √0.771 → 0.937. Both land above the current 0.90–0.93
      bracket, from two routes that do not share the library.
- [ ] **1.5 Put an interval on it.** The tiling comparison is 901 guides; bootstrap
      it. A ceiling quoted without a CI is the same mistake as the cross-organism
      table (Thread 3).
- [ ] **1.6 Test the "average over replicates" claim.** Part 10 asserts that
      predicting the mean of several screens is an easier target and the apparent
      ceiling rises. With per-replicate scores this stops being an assertion:
      train the frozen baseline against one replicate, then against the mean, and
      report both. Cheapest demonstration in the whole thread.

### Deliverable

One short section replacing Part 10's ceiling subsection: a same-enzyme
repeatability estimate with a CI, the ceiling for both targets, and a statement of
where ρ 0.707–0.718 sits against it. Plus a one-line correction to Part 11's
"replicates nobody has published".

### Pitfalls

- Guo's high-quality sets **exclude guides in genes resistant to Cas9 editing**.
  If our labels inherit that exclusion, the positional variance is censored, which
  affects both the ceiling and the flank effect size. Check whether the 13,880 and
  33,567 sets are downstream of that filter.
- Replicate agreement is *not* an upper bound on predictability if the two
  replicates share a systematic bias — the same hole `ceiling2.py` found in the
  two-enzyme estimate. Independent transformations of the same plasmid library
  still share library composition. Say so.

---

## Thread 2 — The long-range downstream asymmetry

### What it entails

Part 5 is the report's strongest biology and it is half-finished. Flank features
are worth +0.080 ρ on the original rows and +0.164 on the curated set, larger than
everything else tested combined; the signal is a smooth gradient rather than a
motif (a CNN on raw flanking DNA gets +0.019, a factor of 4.2 less); it peaks at
250–500 nt; and downstream beats upstream by about 3.4×, with the side rather than
the column budget responsible (+0.0790 count-matched against +0.0782).

Short range has a mechanism — PAM-first engagement and R-loop propagation away
from the PAM. **Long range has none.** Transcription was the obvious candidate and
is not supported: the asymmetry is flat across expression bands (0.0486 / 0.0458 /
0.0479). There are two things to do, and the first one is more important than
finding a mechanism.

### Why it matters

The flank effect is the headline result. Part 7 found that our largest
non-sequence gain, supercoiling, was **97% read depth** — +0.047 total of which
only +0.004 was twisting-specific — and Appendix A records that mappability+GC is
worth +0.038 alone and **−0.000 on top of `a_flank`**. That collapse is usually
read as "the flank family already contains it", and it does. But redundancy is
symmetric and does not identify the causal direction: it is equally consistent
with part of the flank gain being the same sequencing artefact that inflated
supercoiling. Having already overturned one of our own positive results this way,
the headline result should face the same control before it goes in a paper.

### Steps

**2A — the artefact control on the flank effect (do this first)**

- [ ] **2A.1 Does `a_flank` predict control-library abundance?** Guo's activity
      score is a depletion of the Cas9 arm relative to a **dCas9 control arm**, so
      binding-only effects are already divided out — good. What is not divided out
      is composition bias in the control library itself (synthesis, PCR, GC-biased
      amplification), which is exactly what the mappability block was built to
      detect. With the control read counts from Thread 1.3, fit `a_flank` against
      the control abundance. A strong fit means part of +0.080 is assay, not
      biology. **This is the one place where Threads 1 and 2 share a
      prerequisite**, which is the argument for doing the recount.
- [ ] **2A.2 Stratified control, no new data needed.** Restrict to guides in
      high-uniqueness windows using the existing `eng.mapgc.uniq.w*` columns from
      `make_diagnostic_blocks.py`, and re-run the flank ablation inside that
      stratum under grouped CV. If the gain survives on unambiguously mappable
      DNA, the artefact account is largely dead and the result is stronger for
      having been attacked. Runnable today.

**2B — candidate mechanisms for the 250–500 nt scale**

- [ ] **2B.1 Replication-fork co-orientation.** The sharpest available test and
      the machinery exists: `genome.py` already has `replication_landmarks()`,
      `replichore()` and `replication_fork_direction()`, and
      `make_diagnostic_blocks.py` already builds `eng.dose.right_replichore`.
      Replichore membership flips the relationship between genomic coordinate and
      fork direction, so if the asymmetry is set by replication it should **change
      sign or magnitude between the two replichores**; if it is set by the Cas9
      complex or by gene-level composition it should not. Measure the
      upstream/downstream gain separately per replichore, the same shape of
      experiment as `asymmetry.py --run` with replichore in place of expression
      band. A sign flip would be the mechanism; no difference closes another
      candidate, which is also publishable given Part 6's framing.
      Note the prior: Part 7 found an analytic replication gradient worth +0.013
      against abundance's +0.045, so replication timing is already known to be
      part of the story but not most of it.
- [ ] **2B.2 Locus-dependent repair.** Guo observed that regions resistant to Cas9
      in the WT background **become vulnerable in ΔrecA**, and concluded that
      endogenous repair mitigates DSB lethality in a locus-dependent way. Repair
      capacity is a property of chromosomal neighbourhood, which is the right
      length scale. Test: rebuild the flank ablation on the Cas9 ΔrecA label
      (Table S8) and compare the long-range component against the WT Cas9 label
      (Table S6). If the long-range flank gain shrinks when repair is knocked out,
      the effect is about repair rather than cutting.
      **Caveat to state up front:** Guo report ρ 0.328 for their own ΔrecA model
      against 0.542 for Cas9, because selection pressure saturated the assay. A
      smaller gain in ΔrecA is therefore expected from compressed dynamic range
      alone. Normalise by each label's own baseline ρ and report the ratio, not
      the raw Δ — and if the dynamic range is too compressed, report the test as
      inconclusive rather than as evidence.
- [ ] **2B.3 Re-open orientation on the intergenic sub-library.** Per correction 3
      above: check whether both-strand intergenic guides exist in Guo's Table S6,
      and if so extend `asymmetry.py --audit` to report orientation variation on
      that subset. If there is variation, run the orientation-conditioned
      asymmetry test the module currently declares impossible. Smaller n, weaker
      definition of orientation, and still the only direct test of the
      transcription account.

### Deliverable

A replacement for the last paragraph of Part 5 §3. Best case, a mechanism. Worst
case, two more candidates closed with controls and a headline result that has
survived the same artefact test that killed the supercoiling claim. Both are
results under the paper's framing.

### Pitfalls

- Grouped CV is not optional here. `run_ablation.py --group` exists for this
  reason and its help text says so: roughly 20 guides per gene means ungrouped
  folds let the model memorise loci, and every one of these tests is positional.
- Replichore and expression and GC are all correlated with genome position. Any
  positive result needs the count-matched control that Part 5 §3 already uses.

---

## Thread 3 — Cross-species and cross-organism, with error bars

### What it entails

Two separate claims currently carry no uncertainty.

**The 2×2 organism transfer** (Part 5 §5) — *E. coli* WT-SpCas9 against
*C. rodentium* TevSpCas9:

| | tested on *E. coli* | tested on *C. rodentium* |
|---|---:|---:|
| trained on *E. coli* | 0.708 | 0.700 |
| trained on *C. rodentium* | 0.628 | 0.764 |

Each cell is one fit on a full source screen predicting a full target screen.
Part 2 flags it as having no interval over the test set and Part 12 calls it
"the weakest-supported table in the report". Three model seeds exist per cell but
all share one test set, so their spread understates the uncertainty.

**The cross-kingdom result** (Part 10) — *E. coli* → human −0.048 and human →
*E. coli* −0.017, with the mechanism inverting: GC → cut score is −0.201 in
*E. coli* and +0.017 in human, and the human relationship is U-shaped with an
interior optimum at GC 0.571 and a quadratic fit explaining 22× more variance
than the linear one. Also no intervals, and this one is load-bearing for the
mentor's cross-species question and for the benchmark paper's own Pearson 0.0656
cross-prediction.

### Why it matters

These are the numbers a reviewer will attack first, and the fix is mechanical.
A negative ρ of −0.048 is indistinguishable from zero without a CI, and "useless
in both directions" is a much safer sentence when the interval is printed next to
it. Part 12 estimates the bootstrap at about an hour of work.

### Steps

- [ ] **3.1 Write one bootstrap helper.** There is no bootstrap anywhere in
      `Runjia/sgrna/` — I checked. One function taking `y_true`, `y_pred`, a metric
      and `n_boot`, resampling test rows with replacement, returning the
      percentile interval. It belongs next to `top_decile_precision` and
      `expected_pick_percentile` in `evaluate.py` so everything else can use it.
- [ ] **3.2 Bootstrap the 2×2.** Resample the *test* screen, not the training one;
      the uncertainty being quantified is "would this transfer number hold on
      another sample of target guides". Resample by **gene or genomic block, not
      by row** — neighbouring guides are not independent, which is the same reason
      grouped CV exists, and row bootstrap will produce intervals that are too
      narrow. Report seed spread and test-set interval separately; they measure
      different things.
- [ ] **3.3 Bootstrap the cross-kingdom numbers**, including the GC quadratic fit:
      an interval on the turning point (GC 0.571) is more informative than the
      p-value, because the whole claim is that the optimum is *interior* to the
      human data and below the 1st percentile in *E. coli*.
- [ ] **3.4 Fix `guide_overlap` while in `transfer.py`.** It is recorded in every
      output row and is always `NaN`: `transfer_model()` reads a `protospacers`
      key out of `meta.json`, and `build()` never writes one — it writes `y`,
      `names`, `n_base`, `dataset`, `left`, `genome_length`. The docstring says
      overlap is "measured and reported rather than removed, because for the two
      *E. coli* screens it is near-total", and that caveat is exactly what
      distinguishes a *label*-transfer test from a *sequence*-generalisation test.
      Write `protospacers` in `build()` and let the existing intersection compute.
      Then recheck whether any conclusion in Part 5 §5 depended on reading that
      column as if it were populated.
- [ ] **3.5 Keep the honest caveats attached.** Both new screens are
      *Enterobacteriaceae*, the *C. rodentium* screen covers 4.3% of one
      chromosome, and no genome-wide Cas9-cutting screen exists outside this
      family (Part 11). Error bars make the table defensible; they do not make it
      cross-species.

### Deliverable

Both tables reprinted with intervals, the `guide_overlap` column populated, and a
reusable bootstrap in `evaluate.py`.

---

## Thread 4 — Metric discipline, and the within-gene ρ check

### Where I agree with the report

Part 2's position is right and is **not** "ρ instead of R²" — it is report both,
and the report says why: a LightGBM with a robust loss scored **R² 0.192 and
ρ 0.568**, same model and data, because the loss shrank predictions toward the
middle and destroyed the scale while preserving the order. ρ is the right headline
because the decision is ranking candidates, and because the labels are a depletion
proxy on an arbitrary scale — min-max normalised by Noshay et al. across species,
and running at MSE ≈ 79.5 in *E. coli* against ≈ 0.0224 in human, two numbers with
no common unit. ρ is invariant to that; R² and MSE are not.

R² has to stay for two reasons. Any sentence of the form "we can predict
efficiency" needs it. And it is the only column that connects to the benchmark
paper, which reports Pearson and R² and **never reports Spearman** — the reason
re-running their iRF to produce ρ 0.4785 was necessary at all (Part 3).

### The gap

**Global ρ is probably not the ρ a practitioner experiences.** With ~20 guides per
gene, a ρ pooled across a whole fold gets credit for ordering *genes* — gene A is
more cuttable than gene B — but a practitioner has already chosen the target and is
ranking guides inside it. Between-gene ordering is free credit on a decision
nobody makes. This bears directly on the headline: the flank effect is long-range
and positional, so some of +0.080/+0.164 may be gene-level signal that vanishes
once the gene is fixed. Grouped CV tests generalisation to unseen genes; it does
not convert the metric into a within-gene number. Part 2's pick-percentile
simulation has the same exposure — it draws ten candidates from the library, not
ten candidates for one gene.

### Steps

- [ ] **4.1 Within-gene ρ.** Take the existing out-of-fold predictions, group by
      `gene_name` (already in the guide index — `build_guide_index.py` writes it,
      and `run_ablation.py`'s `gene` group option uses it), compute ρ within each
      gene with ≥10 guides, and report the mean with its spread. Near 0.70 means
      the headline is safe and now defensible against the obvious question. A
      large drop is a more interesting finding than anything left in the QCT
      direction.
- [ ] **4.2 Re-run the pick-percentile simulation within gene**, drawing the ten
      candidates from one gene instead of from the library. This is the number that
      belongs in an abstract.
- [ ] **4.3 Repeat 4.1 for the flank ablation specifically**, not just the
      champion: is Δρ from `a_flank` still positive within gene? If the flank gain
      is mostly between-gene, the paper's central claim needs restating as "which
      loci cut well" rather than "which guide to pick" — a true statement either
      way, but a different one.
- [ ] **4.4 Standardise the reporting block.** Every results table carries ρ, R²
      and Pearson together — all three are already logged per fold by
      `FINAL_baseline_ecoli_sgRNA_model_5seeds.py` and `run_ablation.py`, so this
      costs nothing — plus one decision metric. `evaluate.py` already has
      `top_decile_precision`, `expected_pick_percentile` and `calibration_slope`.
- [ ] **4.5 Carry Part 2's warning into the text.** Ratios of Δρ are not ratios of
      information: 0.53 → 0.61 is +0.09 on the ρ² scale while 0.82 → 0.90 is +0.18.
      The "4.2× better than a CNN" comparison is valid because both arms sit at the
      same baseline, and the paper should say that where the claim appears rather
      than in a methods note.

---

## Priority

1. **2A.2** — stratified mappability control on the flank effect. Runnable today,
   needs no downloads, and defends the headline result.
2. **3.1–3.4** — bootstrap helper, the 2×2 and cross-kingdom intervals, and the
   `guide_overlap` fix. Mechanical, about an hour each, removes the weakest table.
3. **4.1–4.3** — within-gene ρ. Cheap, uses artefacts that already exist, and can
   change how the central claim is phrased.
4. **1.1–1.2** — read the supplementary tables and the reliability figures. Decides
   whether the recount is needed.
5. **2B.1** — replichore test. Machinery exists; sharp prediction.
6. **1.3–1.6** — the recount, if 1.1 shows per-replicate numbers are unpublished.
   Largest piece of work here, and it unlocks both the real ceiling and 2A.1.
7. **2B.2, 2B.3** — ΔrecA and intergenic orientation. Both need Guo's tables.

## Not doing

No more regional QCT summaries, coupling terms, PAM-weighted variants or
profile-shape features. The evidence that this line is finished: `regional_qct`
was selected **zero times in all 25 seed-folds** of the latest run, mean R² moved
0.2923 → ~0.298, and Part 9 shows sixteen model classes span only 0.116 ρ while
300 → 1,200 features buys +0.001. Part 6's rule explains why in advance, and it
was correct 8 times out of 8.

Also not doing: `filter_by_control_reads` and the read-count threshold sweep
(correction 1), model tuning or the crisprHAL comparison (Runjia), and the
positions 18–20 SHAP analysis (Joshua).
