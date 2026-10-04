# SLICER — what determines whether a CRISPR guide cuts, in *E. coli*

**S**imple **L**ightweight **I**nterpretable **C**RISPR **E**fficiency **R**anker.

RSI09. State of the project as of 4 October 2026. Written to be read without a
background in molecular biology or machine learning: Part 1 defines every
technical term, Part 2 explains how to read the numbers, and each result says
what it means in plain words before giving the figure.

This document states what we currently think and the evidence for it. It
deliberately contains no record of what we used to think — fourteen claims have
been withdrawn along the way, and those live in a separate working file
(`local/MISCONCEPTIONS.md`, not in the repo) so that this one reads as a
statement rather than a diary. `FINDINGS_LOG.md` is the chronological record and
`../results/RESULTS.md` has the full protocol behind every number.

---

## Contents

- [Part 0 — Why this problem is worth working on](#part-0--why-this-problem-is-worth-working-on)
- [Part 1 — The words you'll need](#part-1--the-words-youll-need)
- [Part 2 — How to read the numbers](#part-2--how-to-read-the-numbers)
- [Part 3 — Where SLICER stands](#part-3--where-slicer-stands)
- [Part 4 — What made it possible: decoding the published dataset](#part-4--what-made-it-possible-decoding-the-published-dataset)
- [Part 5 — The one large effect, and what kind of thing it is](#part-5--the-one-large-effect-and-what-kind-of-thing-it-is)
- [Part 6 — The rule: when a feature cannot help](#part-6--the-rule-when-a-feature-cannot-help)
- [Part 7 — Two effects that are not what they look like](#part-7--two-effects-that-are-not-what-they-look-like)
- [Part 8 — Controls](#part-8--controls)
- [Part 9 — The model, its settings, and what it can and cannot explain](#part-9--the-model-its-settings-and-what-it-can-and-cannot-explain)
- [Part 10 — How much room is left](#part-10--how-much-room-is-left)
- [Part 11 — What is new, what is imported, what is open](#part-11--what-is-new-what-is-imported-what-is-open)
- [Part 12 — Next steps and how to frame the write-up](#part-12--next-steps-and-how-to-frame-the-write-up)
- [Appendix A — Every feature set, grouped, with its data source](#appendix-a--every-feature-set-grouped-with-its-data-source)
- [Appendix B — Papers referred to](#appendix-b--papers-referred-to)

---

## Part 0 — Why this problem is worth working on

### What CRISPR is used for

CRISPR-Cas9 cuts DNA at a chosen place. Once cut, the cell's repair machinery
either breaks the gene (useful for finding out what it does) or pastes in a
replacement (useful for changing it). Four uses, roughly by how established they
are:

1. **Research** — by far the largest. If you want to know what a gene does, you
   break it and see what changes.
2. **Medicine** — Casgevy was approved in late 2023 for sickle-cell disease and
   beta-thalassaemia, the first approved CRISPR therapy, and has since been
   extended to younger patients.
3. **Agriculture and industry** — crop traits, and **engineering bacteria** to
   manufacture insulin, fragrances, biofuels and drug precursors. A large
   existing industry, and it runs on *E. coli*.
4. **Antimicrobials and diagnostics.**

In every one of these somebody has to choose a guide. Choosing badly wastes an
experiment; at scale it wastes a screen of tens of thousands.

### Is a bacteria-specific result worth having?

We can answer this with our own measurement rather than an opinion. **A model
trained on bacteria does not transfer to human cells at all** — in both
directions the correlation is slightly *negative* (−0.048 and −0.017, Part 10).
So nothing here should be sold as relevant to human gene therapy.

What it is relevant to: **bacterial engineering**, where guide choice is a daily
practical problem; **mechanism**, because the reason the human transfer fails is
itself a finding we can quantify; and **method**, because the most transferable
thing here — the rule in Part 6 for telling in advance whether a feature can
help — is not about bacteria at all.

The split is worth stating explicitly in the paper: **the predictions are
bacterial, the methods and reasoning are general.**

---

## Part 1 — The words you'll need

### The biology

**DNA** — a string in a four-letter alphabet: A, C, G, T, normally in a double
helix where A pairs with T and C with G.

**GC content** — the fraction of letters that are G or C. G–C pairs are held by
three hydrogen bonds and A–T pairs by two, so **GC-rich DNA is harder to pull
apart**. This single fact explains a surprising amount of what follows.

**CRISPR-Cas9** — molecular scissors. **Cas9** is a protein that cuts DNA but
does not know where.

**Guide RNA (sgRNA)**, or **guide** — the address label, carrying a 20-letter
sequence. Cas9 cuts wherever the DNA matches it.

**Protospacer** — the 20 letters of DNA being targeted. **Target** is used
interchangeably.

**PAM** — a three-letter signal in the DNA immediately after the target, which
must read **NGG** (any letter, G, G) or Cas9 will not cut. Note *where* it is:
immediately 3′ of the protospacer, i.e. on the **downstream** side. Different Cas
enzymes need different PAMs, which matters in Part 5.

**Upstream / downstream** — conventionally, 5′ and 3′ along the strand the guide
matches. The PAM is downstream.

**Flank / flanking DNA** — the DNA *surrounding* the target. The guide does not
match it and Cas9 does not read it, so the expectation is that it is irrelevant.
The central finding here is that it is not.

**R-loop** — the structure Cas9 forms once engaged: the guide RNA paired with one
DNA strand, the other DNA strand displaced. It forms starting at the PAM and
propagates away from it.

***E. coli*** — a gut bacterium, the standard workhorse. ***C. rodentium*** — a
mouse pathogen and close relative; both are in the family *Enterobacteriaceae*.

**Nucleoid-associated proteins** — bacteria pack DNA with proteins such as HU and
H-NS. This is **not** the same as the **histones** and **nucleosomes** of humans,
animals and plants. A nucleosome is a spool of DNA wrapped around eight histone
proteins, and it physically blocks Cas9. Bacteria have no histones. This matters
in Part 10.

**Supercoiling** — DNA can be over- or under-twisted. Over-twisted DNA is harder
to open.

### The experiment the data comes from

**A screen** — an experiment testing tens of thousands of guides at once. In
bacteria, cutting the chromosome usually kills the cell, so you grow a mixed
population, sequence what survives, and a guide whose cells **disappeared** cut
well.

**Cut score** — the efficiency number a screen produces.

**The catch, which runs through everything:** the screen measures *survival*, not
cutting. Survival mixes cutting with DNA repair, with growth rate, and with
sequencing noise.

**Sequencing read depth** — how many times a given stretch of DNA was read.
Explained properly in Part 7, because mistaking it for biology cost us a result.

**Data cleaning / curation** — here specifically: discarding guides whose score
came from too few sequencing reads to be reliable. A guide read 5 times has a far
noisier score than one read 500 times. The dataset we compare against discards
**45%** of guides this way. It does not change any remaining score.

**CRISPRi** — a related technique using a *disabled* Cas9 that binds without
cutting, switching a gene off. A different measured quantity, so CRISPRi datasets
cannot substitute here. Worth knowing because most bacterial screens outside
*E. coli* are CRISPRi (Part 5).

### The modelling

**Feature** — one number describing a guide. A model never sees the guide itself,
only its features.

**Feature set / family** — a group of features from one idea and one data source.
Nine of them; Appendix A.

**Model** — shown thousands of guides with their measured scores, predicts the
score of one it has not seen.

**Cross-validation** — hide part of the data, train on the rest, test on the
hidden part, repeat. **Grouped** cross-validation is the stricter version used
for anything positional: hide a whole contiguous chunk of chromosome, because the
screen puts ~20 guides in every gene and a model told *where* a guide sits can
otherwise score well by memorising "guides around here do about this well".

**Ablation** — add one feature set, re-measure, see whether the model improved.

**Permutation control** — re-run a feature set with its rows **shuffled**: the
features are intact but no longer describe the right guide. Whatever that scores
is the set's noise floor.

**Hyperparameters** — a model's dials, set before it sees data: how many trees,
how deep, how fast it learns. **Tuning** means searching for good values.

**Gradient boosting (LightGBM, XGBoost)** — many small decision trees, each
correcting the previous one's errors. Fast, CPU-only, named features.

**Neural network (CNN, BiGRU)** — a more flexible family learning from raw data.
A **CNN** scans for short recurring patterns. Usually needs a GPU and is hard to
interpret. The best published competitor is one.

---

## Part 2 — How to read the numbers

This part exists because the main metric is easy to over-interpret, and we have
done so at least twice.

### Spearman ρ, and what it is not

**Spearman ρ** measures whether the model gets the **ordering** right: rank the
guides by prediction, rank them by true score, and correlate the ranks. 0 is
random, 1 is perfect. It is the right headline because nobody needs a guide's
exact efficiency — they have candidates and want the best one.

But ρ is a correlation, and **differences in ρ are not differences in an amount
of anything.** Three specific errors to avoid:

**"4× the increase in ρ means 4× better" — no.** A ratio of two Δρ values is not
a ratio of information, accuracy, or usefulness. If you want a quantity that
behaves additively, ρ² is closer (it is roughly a share of variance), and
−½·log(1−ρ²) is closer still (it is in bits). On the ρ² scale, going 0.53 → 0.61
is worth +0.09 while 0.82 → 0.90 is worth +0.18 — **the same Δρ is worth twice as
much higher up.** So when this report says the flank features are worth 4.2× what
a CNN extracts from the same DNA, that is a statement about two measurements
taken at the same baseline, and it does not generalise to comparisons at
different baselines.

**"We are 70% of the way to the ceiling" — only on a stated scale.** 0.707/0.90
is 79% of the way on the ρ scale. On the scale that matters for a decision it is
different, and the honest thing is to say which you mean. Simulating the decision
directly (take ten candidate guides, let a model of given ρ choose, record where
its pick really falls):

| model ρ | pick lands at |
|---:|---|
| 0 (blind) | 50th percentile |
| 0.53 | 73rd |
| 0.61 | 77th |
| **0.71** | **81st** |
| 0.90 | 88th |
| 1.00 | 91st |

Two things to take from this. **The whole range from useless to perfect spans
only the 50th to the 91st percentile** — a perfect model does not hand you a
perfect guide, because ten random candidates may not contain a great one. And
**most of the available benefit arrives early**: the first 0.53 of ρ buys 23
percentile points, the next 0.47 buys 18. So on the decision scale we are
(81−50)/(88−50) = **82%** of the way to the ceiling, which happens to be close to
the 79% the ρ scale gives — coincidence at these values, not a rule.

**ρ versus R².** R² is the share of *variance* of the actual values explained; ρ
is agreement on *ranks*. R² therefore cares about predicting the right scale
(calibration) and is sensitive to outliers; ρ ignores both. They can disagree
violently. Our own clearest case: a LightGBM with a robust loss function scored
**R² 0.192 and ρ 0.568** — the same model, on the same data, looking broken on
one metric and respectable on the other, because the robust loss shrinks
predictions towards the middle and destroys the scale while preserving the order.
**Report both.** For choosing a guide, ρ is the relevant one; for claiming you
can predict efficiency, R² is.

### Statistical significance, and the right comparison

A single fold's ρ is less precise than it looks. With a validation fold of
n = 6,713, the standard error of ρ is about 1/√(n−3) = **0.0122**, so a single
fold's 0.707 carries a 95% interval of roughly **[0.683, 0.731]** — ±0.024, wider
than most effects in this project.

**That is the wrong interval to use for comparing models**, though, and using it
would make almost nothing here significant. Two models evaluated on the *same*
folds share all the fold-to-fold difficulty, so the comparison should be
**paired**: take the difference per fold and test that. Our head-to-head against
crisprHAL 2:

| | |
|---|---|
| per-fold differences | +0.0089, +0.0091, +0.0142, +0.0090, +0.0125 |
| mean | **+0.0107** |
| standard error of the mean | **0.0011** |
| paired *t* | p = 0.0006 |
| Wilcoxon signed-rank | p = 0.0625 |

**The paired standard error is eleven times smaller than the single-fold
interval**, which is the whole reason paired designs are used. Note also the two
*p*-values disagree by two orders of magnitude: with n = 5 the signed-rank test
cannot go below 0.0625 however large the effect, so it is the *t*-test that has
power here and the sign test that bounds how much a rank-based reading can claim.
Both are quoted throughout for that reason.

**And statistical significance is not importance.** The +0.0107 above is highly
significant and, as Part 3 explains, still not a defensible claim of superiority,
because an effect the same size is available from tuning either model.

### Do all the numbers have error estimates?

Honest audit:

| kind of number | replication | error estimate |
|---|---|---|
| headline ablations and head-to-head | 5–25 folds, 1–5 seeds | **yes** — fold sd, and paired tests against the relevant arm |
| transfer within a screen | 5 folds | **yes** — fold sd and paired p |
| cross-organism model transfer | **one fit per cell** | **no** — single number, no interval |
| correlations of one column with the label | n = 13,880–59,489 | **implicit** — SE ≈ 1/√n ≈ 0.004–0.008, not usually quoted |
| correlations *between* columns | same | **not quoted** |
| the redundancy and cross-redundancy R² | 5-fold out-of-fold | partial — a distribution over columns is reported, not an interval per column |
| the ceiling | — | **now a range**, Part 10 |
| timings | 5 repetitions, pinned threads | **yes** — IQR |

**The weakest cell is cross-organism transfer**: each train→test figure is a
single fit, so the 92%/89% retentions have no interval and should not be compared
with each other. Fixing it means repeating with several seeds and bootstrapping
the test set — cheap, not yet done, and listed in Part 12.

### How correlation between columns is calculated

Used in several places (features correlating 0.90–0.96 with GC; read-depth tracks
correlating −0.70 with mappability; attribution in Part 9), so worth stating
once. Unless a result says otherwise it is the **Pearson** correlation: subtract
each column's mean, divide by its standard deviation, take the average product
over guides. For a pair of columns *x* and *y* over *n* guides that is
Σ(xᵢ−x̄)(yᵢ−ȳ) / (n·sₓ·s_y). Where the relationship is clearly non-linear, or a
column is an indicator, **Spearman** (the same formula on ranks) is used instead
and is named. Missing values are mean-filled before correlating, which biases
estimates slightly towards zero for sparse columns.

---

## Part 3 — Where SLICER stands

| | what it is | Spearman ρ | guides |
|---|---|---:|---|
| Noshay et al. 2023 | the paper SLICER inherits its features from | 0.502 (Pearson) | 40,468 |
| crisprHAL 2 (2026) | the best published bacterial model | **0.697 ± 0.007** | 33,495 curated |
| **SLICER** | this project | **0.707 ± 0.008** | 33,567 curated |
| **SLICER, tuned** | after searching its settings | **0.718** | 33,567 curated |
| apparent ceiling | Part 10 | 0.90–0.93 | — |

### The head-to-head, and why it is parity

Earlier numbers in this project were measured on different guides from the model
they were compared against. That was fixed in two steps. Because of Part 4 the
published feature set can be computed for *anyone's* guides, so SLICER was run on
crisprHAL 2's own 33,567 curated guides with their label: **ρ 0.7078**. Then
rather than trust their published figure, their actual model — their
architecture, hyperparameters, 48 epochs, code imported unchanged — was re-run on
our exact folds: **0.6971 ± 0.0067**. As a check that the harness is faithful
rather than flattering, it was also run on their own split and scored **0.6940**
against their published 0.695.

So SLICER is ahead by **+0.0107**, in 5 of 5 folds, paired p = 0.0006.

**That is not a claim of superiority, because neither model was tuned.** Searching
SLICER's own hyperparameters — inside the training folds only, so the choice never
sees test data — gains **+0.0103 ± 0.0030** (5/5 folds, p = 0.0015). **Tuning
alone moves our number by as much as the entire margin.** Tuning theirs costs
about a day of CPU and was not done.

**The defensible claim is parity**: two models built on completely different
principles landing within a hundredth of each other on identical data, with the
difference between them no larger than the difference between tuned and untuned
versions of either. The tuned 0.718 may be quoted as what it is — tuned against
untuned.

### What they do differ in: cost

Measured on one machine, CPU to CPU, threads pinned, five repetitions:

| | SLICER | crisprHAL 2 |
|---|---|---|
| one fold | **17.1 s** (IQR 0.7) | ~21 min |
| five folds | **1.4 min** | **105 min** |
| peak memory | 4.03 GB | 4.55 GB |
| features | named quantities you can look up | learned, inside a network |

**About 74× less time for equal accuracy, at essentially the same memory.** The
saving is time, not footprint — SLICER holds a 33,567 × 6,517 table in memory and
is not light on RAM.

> **Could a different model fix the footprint? Mostly no — measured, and not
> adopted** (`results/memory_comparison.csv`). Peak RSS, threads pinned,
> separate processes: XGBoost selector 3.35 GB, LightGBM selector 2.93 GB, so
> **the model is worth ~13%**. The obvious next move backfires: building the
> dense matrix and *then* converting it to a sparse one costs *more*
> (3.51–3.64 GB), because both copies exist at once. Reading the rows in chunks
> and never materialising the dense array reaches **1.98 GB at identical
> accuracy (ρ 0.7011)** — a 41% saving, but a change to the data path that
> touches every caller, so it was left unimplemented. The larger lever is not
> memory engineering at all: Part 4 shows **427 of the 6,232 columns are
> statistically indistinguishable from all of them**, and 427 columns is 23 MB
> where 6,232 is 350 MB.

Two caveats. **crisprHAL 2 as published is GPU-trained**, so this shows their
architecture needs ~74× more *CPU* time, not that a laptop beats a GPU. And
**most of our 17 s is not the model** — 14 s is a feature-selection step whose
only output is a ranking. Swapping its XGBoost for LightGBM makes the fold 6.6 s
(≈190×) with accuracy indistinguishable (+0.0008), but it has not been adopted,
because every number here was produced with the current selector.

---

## Part 4 — What made it possible: decoding the published dataset

### The problem

The project began from a 2023 paper describing 13,880 guides with 6,232 numbers
each, mostly quantum chemistry. Two things made it unusable by anyone else: **it
never recorded the DNA sequence of any guide**, and **5,887 of the 6,232 columns
were named `V1`, `V2`, `V3`…** with no description. So the method could only ever
be applied to guides its authors had already processed.

### What we found

**The sequence is hidden in the numbers.** One column records the electron count
of the DNA letter at each position, and across all 13,880 guides it takes exactly
four values — 42, 48, 50, 56 — the valence-electron counts of C, T, A and G. That
column *is* the sequence in another alphabet.

**It checks out against the real genome.** 13,879 of 13,880 decoded; every one
found in the *E. coli* genome with the required NGG PAM in the right place;
13,877 appearing exactly once. Independently, a standard method for locating
where bacterial DNA copying starts put it at position 3,923,620 against a
textbook value near 3,923,800.

**The quantum tables are recoverable exactly.** Each quantum column is a fixed
lookup from a short run of letters to a number. All 4 single letters, 16 pairs, 64
triples and 256 quadruples appear **with no contradictions anywhere** — so the
parameterisation is recovered, not approximated, and can be applied to DNA the
paper never touched. That is what Part 5's feature set does.

**The anonymous columns were identified.** 5,853 of the 5,887 are "is there letter
X at position Y" indicators, none ambiguous, verified against **1.17 million
rebuilt values at 100% agreement**.

### How few of the 6,232 columns does the model actually need?

Decoding says the matrix is a re-spelling of 20 letters. That is an argument
about *information*, and it does not by itself mean the model can be handed
less: a tree can only split on what it is given, and a redundant encoding might
be the shape that makes the signal reachable. Nobody had tested it, so we did
(`results/representation*.csv`). Same rows, same cross-validation, same model;
only the columns change. Seed-to-seed spread is 0.001–0.003, so a difference
above ~0.006 is real.

| what the model is given | columns | ρ |
|---|---:|---:|
| everything published | 6,232 | 0.5272 |
| the 20 letters, one indicator per letter per position | **74** | **0.5081** |
| …plus letter-pair indicators | 364 | 0.5189 |
| only the quantum-chemistry columns | 316 | 0.5154 |
| only the 5,853 anonymous indicators | 5,853 | 0.5119 |
| only the 63 hand-named columns | 63 | 0.3044 |

**74 columns recover 96.4% of what 6,232 do**, and those 74 are a complete
encoding — every position carries at least three of its four letter indicators,
so the fourth is implied and the guide's sequence is recoverable exactly. That
is an 84× reduction for 0.019 Spearman.

Two readings of the table are worth stating plainly. **Giving the model the
5,853 anonymous indicators scores *worse* than giving it 364 columns** — more
columns, worse model, because the extra triple- and quadruple-letter indicators
are restatements that dilute feature selection. And **422 columns beat all
6,232**: the 74 letter indicators plus Part 5's flanking-DNA family reach 0.5462,
+0.019 ahead of the entire published matrix. One family describing DNA *outside*
the guide outweighs the whole published feature set.

**Where the last 0.019 comes from — after two wrong guesses.** It cannot be
information, since nothing in 6,232 columns is absent from the 74. The first
explanation tried was that the redundant columns pre-compute combinations our
trees are too shallow to build, which predicts that deeper trees close the gap.
They do the opposite, at every encoding: 0.5081 → 0.4953 → 0.4739 → 0.4557 as
depth goes 4 → 6 → 8 → 12. Supplying explicit triple-letter indicators hurts
too. The gap is not about combinations.

What accounts for it is the **63 hand-named columns**, which are the one part of
the matrix that is *not* derived from the guide's 20 letters — they include how
far along its gene the guide sits, and letter counts that disagree with a direct
count by up to 13. Added to the 364-column encoding they give **427 columns at
ρ 0.5249, statistically indistinguishable from all 6,232** (paired across five
seeds, mean difference 0.0023, p = 0.11).

So the matrix is **74 columns of sequence, 63 columns of something that is not
sequence, and about 6,095 columns of restatement.** A 427-column replacement
exists; it is not the default only because the 6,232-column form is what the
reproduction check in Part 3 anchors to. One result we cannot explain and are not
claiming to: adding the 63 columns to the bare 74 makes things *worse* (0.5008),
and the sign flips only once letter-pair indicators are present.

### Why it matters

**99.0% of the published dataset — 6,169 of 6,232 columns — now regenerates from
a guide's 20 letters alone.** A resource that worked only for its authors works
for anyone, on any bacterial CRISPR screen. The head-to-head in Part 3, the
cross-organism tests in Part 5 and the human boundary in Part 10 all depend on it.

---

## Part 5 — The one large effect, and what kind of thing it is

One argument in six steps: the DNA around a cut site predicts cutting; the signal
is a smooth gradient rather than a pattern; it is stronger downstream for reasons
we can partly name; it belongs to the DNA rather than the enzyme; it crosses
between organisms; and it has a mechanism.

### 1. The surrounding DNA predicts cutting, and it is the largest effect found

Features describing the flanks are worth **+0.080 Spearman** on the original data
and **+0.164** on the cleaned version — larger than every other idea tested here
combined.

The features are deliberately simple: for windows of 50, 250, 500 and 1,000
letters each side, the GC fraction, purine fraction, longest run of one letter,
and averages of the recovered quantum tables; plus the identity of each of the 10
letters immediately either side.

### 2. The signal is a gradient, not a pattern — and that is testable

If the effect were a short recurring motif, a convolutional network should find
it better than hand-computed averages, because that is what CNNs are for. We
built one with the same architecture family as the best published competitor and
gave it raw DNA.

On the **same rows and label**, ±100 letters of raw flanking sequence is worth
**+0.019 ρ** to the network; the hand-computed windowed composition is worth
**+0.080** — a factor of **4.2** from the same DNA. (Both at the same baseline;
see Part 2 on why that ratio does not generalise to other baselines.)

So there is no motif to find. The signal is closer to "how GC-rich are the next
500 letters" — a smooth average that a windowed mean computes exactly and a
5-letter pattern detector must approximate badly. **Using the right instrument for
the shape of the signal is the whole of SLICER's advantage.**

The effect **peaks at 250–500 letters** and has faded by 1,000, so there is no
point looking further out.

### 3. Why downstream matters more than upstream

Downstream is worth about **3.4×** upstream. Two separate explanations are needed,
because the effect has two length scales.

**For the nearest ~30 letters, enzyme geometry is a plausible cause.** Cas9's
engagement with DNA is strongly asymmetric. It recognises the **PAM first**, and
the PAM sits immediately downstream of the target. Only after PAM binding does it
unwind the duplex, and the R-loop then propagates **away** from the PAM. The
displaced strand exits on the PAM side. So the bases just downstream are where
Cas9 first makes contact and where the displaced strand is threaded — and the
published bacterial model independently found downstream context helping out to
about +8 letters while upstream did not help at all. **This is a mechanistic
expectation consistent with our data, not something we tested.**

**For the 250–500-letter scale it cannot be enzyme geometry.** Cas9 contacts
nothing that far away. A compositional gradient over hundreds of bases has to be
a property of the DNA's local state or of the assay. One concrete hypothesis we
have *not* tested: genes have a direction, so "downstream of the PAM" is
systematically related to the guide's orientation relative to the gene it sits in
— which would make the asymmetry a transcription effect rather than a structural
one. Testing it means repeating the upstream/downstream split while conditioning
on whether the protospacer lies on the coding or template strand, which is cheap
and listed in Part 12.

**Until that is done, the asymmetry is a solid observation with a partial
explanation**, and should be written that way.

### 4. It is a property of the DNA, not of this particular enzyme

Three usable screens:

| screen | organism | enzyme | guides | base | + flanks | gain |
|---|---|---|---:|---:|---:|---:|
| WT-SpCas9 | *E. coli* | SpCas9 | 33,567 | 0.544 | **0.707** | +0.163 |
| eSpCas9 | *E. coli* | eSpCas9 | 59,489 | 0.685 | **0.789** | +0.104 |
| TevSpCas9 | *C. rodentium* | TevSpCas9 | 25,210 | 0.704 | **0.764** | +0.060 |

All under strict grouped cross-validation, all winning 5 of 5 folds (p = 7×10⁻⁶
and 3×10⁻⁵ for the two new ones).

A fourth screen, **TevSaCas9**, was excluded. Its enzyme needs a different PAM
(NNGRRT rather than NGG), so only 45% of its sites carry the NGG our pipeline
assumes and the rest would be silently mis-positioned. **The pipeline would have
run to completion with no error and returned plausible results measuring the wrong
positions** — more dangerous than a crash, and the reason for leaving it out
rather than caveating it.

**A limitation to state rather than let a reader find.** eSpCas9 is WT-SpCas9 with
three point mutations, and the two screens use the **same guide library** — 100%
identical sequences. So "two enzymes agree" is weaker than it sounds. The
genuinely different enzyme, TevSpCas9, is also in a different organism, so that
arm confounds the two. **The cross-enzyme evidence is real but narrow.**

### 5. A model trained in one organism ranks another organism's guides

Set up as a 2×2 using the only pairing the data allows — *E. coli* WT-SpCas9 and
*C. rodentium* TevSpCas9. Each cell is a single fit on the full source screen
predicting the full target screen, so these figures have **no error bars** (Part
2) and should not be compared with each other:

| | tested on *E. coli* | tested on *C. rodentium* |
|---|---:|---:|
| **trained on *E. coli*** | **0.708** | 0.700 |
| **trained on *C. rodentium*** | 0.628 | **0.764** |

Reading each off-diagonal against the diagonal below it: *E. coli* →
*C. rodentium* keeps **92%** of a locally-trained model, and the reverse **89%**.
The guides share no sequence and the organisms no chromosome.

**The flank features are the portable part.** Without them the same table reads
0.626 and 0.478 — the inherited representation transfers poorly, and nearly all
retention comes from the flanks. Crossing out of *C. rodentium*, 88 long-range
columns beat 260 short-range ones (+0.107 against +0.058).

**Why the enzyme necessarily differs between cells.** No screen uses the same
nuclease in both organisms. That this is tolerable is measurable rather than
assumed: within *E. coli*, training on WT-SpCas9 and testing on eSpCas9 retains
**87%** — the same band as the cross-organism 89–92% — so changing enzyme and
changing organism cost about the same, and neither dominates. The WT→eSp arm is
**not** a generalisation test, though (those screens share all their guides, so
every test sequence was in training); it bounds enzyme sensitivity only.

### 6. How a borrowed model sees something the local data cannot

This looks paradoxical and is worth spelling out, because it was our own first
objection. The *C. rodentium* screen covers **229 kb — 4.3% of one chromosome** at
110 guides per kb, which is the authors' stated design. Inside it, flank GC
correlates **−0.022** with the label, against +0.159 and +0.138 in the two
*E. coli* screens. Yet an *E. coli*-trained model transfers the effect in
successfully.

The resolution is that **each guide still has its flanks; what the screen lacks is
variation between guides.** Every *C. rodentium* guide's ±1 kb context was read
off that organism's reference genome, so the features exist and are perfectly
well-defined for all 25,210 of them. But because all the guides sit within one
229 kb window, the *range* of long-range composition across the screen is narrow,
and a relationship cannot be *estimated* from a predictor that barely varies.
Applying a relationship learned elsewhere needs no variation at all — only a
value per guide.

So the limitation is **statistical, not informational**: the sequence shown to the
model is not too short, it is too uniform. **Discovery and verification have
different data requirements**, and that distinction is worth a sentence in the
methods.

### 7. The mechanism: which step is the bottleneck

GC-rich *targets* cut worse (−0.20 to −0.15), and a guide binding its target
*more* strongly also cuts worse (−0.21).

The second is the informative one. Before the guide can pair with its target, the
helix has to be **pried apart**. If grabbing on were the slow step, stronger
binding would help; because it *hurts*, the slow step must be the prying apart —
and GC-rich DNA, with its extra hydrogen bond per pair, is harder to pry apart.

Cutting here is **strand-invasion limited, not hybridisation limited.** This
replicates in all three bacterial screens and is *strongest* in *C. rodentium*
(−0.379), so unlike the flank gradient it does not depend on a screen's genomic
span.

---

## Part 6 — The rule: when a feature cannot help

Six of nine feature sets improved the model by less than 0.01. That is the most
transferable result here.

### The rule

Every one of the 20 target positions in the published dataset carries indicators
for at least three of its four possible letters — and if a position is not three
of them it must be the fourth. So **the dataset already encodes the guide's
sequence completely.**

Hence: **a feature computable from information the model already has adds no new
information.** It can only restate it in a shape the model may find easier.

Note the phrasing: *what the model already has*, not *the guide's 20 letters*.
The second is a special case, and the difference decides the hardest case below.

### Tested as a prediction, not an explanation

Explaining six failures with a rule invented after seeing them is weak. So here
is the rule applied **as a forecast**, with the reasoning fixed before the result
column is read. Two questions: can this be computed from what the model already
holds? If not, is the outside information actually about the locus?

| feature set | computable from what we already have? | **predicted** | **measured Δρ** | right? |
|---|---|---|---:|---|
| `b_energy` — binding energy | yes, a sum over letter pairs | nothing new | +0.007 | ✓ |
| `d_mechanics` — duplex stability | yes, also over letter pairs | nothing new | +0.020 | ✓ |
| `c_folding` — RNA folding | yes, the guide folds alone | nothing new | +0.007 | ✓ |
| `f_methylation` — methylation motifs | yes, a text search | nothing new | +0.003 | ✓ |
| `i_shape` — DNA bendability | **largely yes** — median R² 0.52 from flank composition | little on top of flanks | **+0.001** on top of flanks | ✓ |
| `a_flank` — flanking composition | **no** | genuinely new ⇒ should help | **+0.080** | ✓ |
| `d_supercoiling` / `e_nucleoid` / `g_transcription` | **no**, each needs a measurement | new ⇒ should help | +0.045 / +0.016 / +0.015 | ✓, but see Part 7 |

**Eight of eight, once the rule is stated correctly.**

**`i_shape` is the case that sharpened it.** DNA bendability is outside the
20-mer, so the narrow version of the rule predicted it would help — and alone it
does (+0.036). On top of the flank features it adds +0.001 against a seed spread
of ±0.003. Asking the right question explains it: regressing the 152 shape
columns on the 348 flank-composition columns gives a **median R² of 0.52**, half
of them more than half explained, by a *linear* model — and the reverse direction
gives only **0.08**. So composition largely determines shape and not vice versa:
shape is a lossy re-description of something already present.

Two things are worth separating, because only the first is about redundancy.
Median 0.52 means roughly half of a typical shape column is *not* linearly
explained by composition. That residual half is simply uninformative about
cutting — independently confirmed by the residual test, where the best of all 152
shape features correlates with the model's errors at 0.028, below the 0.035 noise
level of features it already uses, **including** the one component that letter
composition provably cannot express (a phased bend sum at the helical repeat, at
0.025). So: much of shape *is* composition, and the part that isn't doesn't
matter.

### What the rule does not cover

**`h_offtarget` is outside its scope.** Off-target burden needs genome-wide
information, so the rule cannot rule it out — but the rule says when a feature
**cannot** help, never that new information **must** help. Its measured gain is
+0.0003 against its own permutation floor of +0.0004, so it is **a null result**.

Is there a biological reason to have expected an effect? A weak one. A guide with
near-matches elsewhere might cut at those sites too, and extra cuts kill the cell,
which would inflate apparent depletion; and a guide targeting a repeated sequence
cuts multiple copies. Both are real mechanisms, but *E. coli* has few exact
repeats and the screen's guides were largely chosen to be unique, so there was
little room for either. **So this is close to confirming that something expected
to have no effect has no effect** — worth one line as a documented null, not a
test of the rule.

### Why `d_mechanics` is the exception worth reporting

Among the redundant group, `d_mechanics` gains +0.020 — three times the others and
153× its own floor. A physical reparameterisation of information already present
*can* help a tree, because it says *where* along the guide the duplex is weak
rather than how GC-rich it is overall.

So "a better-shaped version of the same information" has now been tested twice and
gone both ways: it won for `d_mechanics` (+0.020 over raw GC) and lost for
`i_shape`. Reporting both is more useful than quoting only the success.

### The practical upshot

Run the check **before** building a feature. It takes seconds and would have saved
four of nine feature sets from being built. Measured recoverabilities: binding
energy **98%** from the guide alone by a straight line; duplex stability **85%**;
RNA folding 30%; methylation motifs 14%; bendability **52% from the flank
features**; everything that helped, **≈0%**.

---

## Part 7 — Two effects that are not what they look like

### 1. The "supercoiling" effect is not supercoiling

Our largest non-sequence gain (+0.047 R²) came from a dataset measuring DNA
twisting. Then we ran the source experiment's own **negative control** — the same
measurement with the biological part deliberately removed, which should contain
nothing. **It predicted cutting just as well.**

| what the model was given | ΔR² |
|---|---:|
| the whole feature set | +0.047 |
| read depth only, twisting removed | **+0.045** |
| **only** the twisting-specific part | **+0.004** |

97% of the effect is read depth. **The write-up must not say "supercoiling
predicts sgRNA efficiency."**

### 2. What read depth is, and why it correlated at all

**Read depth** is how many sequencing reads came from a given stretch of DNA. To
measure a screen you sequence a pool of cells and count how often each guide's
barcode appears; separately, experiments like the one above sequence the
chromosome itself and count reads per position. Depth varies for reasons that have
nothing to do with cutting: GC-rich fragments amplify less efficiently in the PCR
step, and reads from repeated sequence cannot be assigned to one location.

**So why should it predict anything?** Not because depth causes cutting — it
cannot. The chain is a confound: **GC content and repetitiveness affect both the
measured read depth and the measured cut score**, the second because the cut score
is itself derived from counting sequencing reads. Two quantities sharing an
upstream technical cause correlate without either causing the other.

The evidence is direct. The read-depth tracks correlate with local GC and with
25-letter uniqueness at up to **−0.70**, and **8 columns computed from the
reference genome alone reproduce the entire gain** (+0.038 against +0.036).

> **What those 8 columns are.** Two quantities at four window sizes (100, 500,
> 2,000 and 10,000 letters centred on the guide): the fraction of positions whose
> surrounding 25-letter sequence occurs **exactly once** in the genome — 97.35% of
> the genome passes, and the failures are repeats — and the window's GC fraction.
> 2 × 4 = 8. No sequencing data at all.

**And then it collapses too.** Stacked on the flank features, the genome-derived
version adds **−0.0002** (8 of 15 folds, p = 0.84). The sequence-intrinsic
explanation of read depth *is* windowed flank composition under another name.

What survives is small and real: the *measured* depth still adds **+0.0061 ρ,
winning 15 of 15 folds** (p = 2.5×10⁻⁵), and that part is not reproducible from the
genome. About a tenth of what the standalone +0.045 implied.

### 3. Four feature sets are one variable

Four of the nine trace to the same quantity. `d_supercoiling` (read depth),
`e_nucleoid` (3D packing from Hi-C), `g_transcription` (gene activity) and the
genome-derived mappability block have leading features correlating 0.77–0.96 with
each other, because Hi-C contacts and ChIP coverage are *both* read counts
inheriting the same GC and mappability bias. Individually +0.047, +0.016, +0.015 —
three findings. **Together: +0.051.** And all of it then collapses into the flanks.

Appendix A groups all nine this way, which is how they should appear in the paper:
**four distinct ideas, not nine.**

---

## Part 8 — Controls

Two different kinds, and the distinction matters.

### Permutation controls — run on every family with a gain

| feature set | real Δρ | its own floor | ratio | verdict |
|---|---:|---:|---:|---|
| `a_flank` | +0.0795 | −0.0009 | 88× | clear |
| `d_mechanics` | +0.0201 | +0.0001 | 153× | clear |
| `b_energy` | +0.0074 | −0.0001 | 92× | clear |
| `c_folding` | +0.0074 | +0.0010 | 7.7× | clear |
| `f_methylation` | +0.0034 | +0.0011 | 3.0× | **marginal** |
| `h_offtarget` | +0.0003 | +0.0004 | 0.8× | **a null result** |

Each floor is that family's own — it depends on how many columns the family
contributes and how many survive selection by luck — and every delta is paired
fold by fold so both arms see identical data.

One sobering detail: **20 of `a_flank`'s *shuffled* columns were still chosen by
the model as "important"**. A feature being selected proves nothing.

### Source-experiment controls — available for exactly one data source

A permutation control asks "better than noise?". It cannot ask "is my biological
interpretation right?" — that needs a control built into the original experiment.

| data source | mock/negative control? |
|---|---|
| GapR-seq (supercoiling) | **yes** — untagged, no-antibody and rifampicin arms. Used, and it overturned the result |
| Hi-C (`e_nucleoid`) | no mock arm published. Checked only indirectly, via its 0.77–0.96 correlation with tracks the untagged arm had already discredited — weaker evidence, and described as such |
| RegulonDB / PRECISE-1K | not that kind of data — annotation, not a measurement with a control |
| reference genome (`a_flank`, `f_methylation`, `h_offtarget`, `i_shape`, mappability) | not applicable — no experiment to control |

---

## Part 9 — The model, its settings, and what it can and cannot explain

### Model choice barely matters

Sixteen model classes on identical data and features:

| | ρ | | ρ |
|---|---:|---|---:|
| stacked combination | **0.611** | PLS | 0.575 |
| LightGBM | 0.609 | small neural net | 0.572 |
| XGBoost | 0.608 | extra trees | 0.556 |
| CatBoost | 0.601 | random forest | 0.550 |
| hist gradient boosting | 0.600 | LightGBM, robust loss | 0.500 |
| support vector machine | 0.586 | nearest neighbours | 0.495 |
| ridge regression (a straight line) | 0.577 | | |

The whole range is 0.116, and on the decision scale about 3 percentile points.
**A straight line gets within 0.03 of the best**, so the effects mostly add up
rather than interacting. **Nearest neighbours fails**, so guides with similar
features do not have similar scores — consistent with many small independent
effects. **Combining models wins +0.002 at 48× the cost**, with the random forest
given a weight of −0.003.

### Settings: tuning, and what happens with more features

Both belong together, because they are the same question — how much does the
model's configuration matter, given the features are fixed?

**Tuning.** Twelve random draws over capacity (`num_leaves`, `max_depth`,
`min_child_samples`), learning rate, `n_estimators`, subsampling and `lambda_l2`.
The search runs **inside each training fold** on an inner split, and the chosen
configuration touches the validation fold exactly once — picking by validation
score would report the best of twelve draws on the test set, which is not a
held-out number. Feature selection is done once per outer fold and shared across
draws, since it is 14 s of the 17 s fold and does not depend on these settings.

| | ρ |
|---|---:|
| default settings | 0.7078 |
| tuned | **0.7181** |
| gain | **+0.0103 ± 0.0030, 5/5 folds, p = 0.0015** |

Per fold: +0.0087, +0.0125, +0.0143, +0.0077, +0.0081.

**More features.** The model is capped at the top 300. Raising it:

| features allowed | ρ |
|---:|---:|
| 300 | 0.6068 |
| 600 | 0.6079 |
| 1,200 | 0.6082 |

**+0.001 for four times as many.** Beyond a few hundred the extra columns
correlate with ones already in and produce no new splits worth making.

Put together: **the configuration is worth about +0.010 and the feature budget
about +0.001, against +0.164 for the right feature set.** Settings are a
rounding error next to features — which is the whole argument for spending effort
where this project spent it.

### And the model is not leaving information unused

If a model were failing to exploit a feature it holds, that feature would still
correlate with the model's remaining errors. Of 6,480 features, **none** exceeds
0.10 (the largest is 0.035). It has squeezed its features dry, which is why model
choice costs so little.

### How a gradient-boosted tree uses features

Many small decision trees in sequence. Each asks a few yes/no questions ("is
downstream GC above 54%?"), each splitting guides into two groups whose average
scores differ as much as possible, and each new tree is fitted to the *errors* of
those before it. Each tree uses only a few features, and a feature that never
produces a good split is never used.

### Why LightGBM, and what the importance disagreement means

LightGBM is not chosen for accuracy — it ties XGBoost at 0.609 — but because its
explanation is reproducible:

| | LightGBM | XGBoost |
|---|---:|---:|
| do its two importance methods agree? | **+0.44** | **−0.02** |
| same features chosen on a different split? | 0.72 | 0.49 |
| features for half the importance | 151 | 330 |

That −0.02 looks alarming, so we tested what it means. Letting each model choose
its own 300 features and comparing the top 50: **14 of 50 shared**, overall
importance rankings correlating +0.69, held-out scores 0.627 and 0.622 — tied. So
they really do choose differently and it really costs nothing.

**It is not that they found different biology.** For each column only one model
picked, the closest counterpart in the other's set has median |correlation|
**0.456**, against **0.069** for randomly chosen columns.

**But substitution is not one-to-one either** — only 22% of unshared picks have a
counterpart above 0.7. A tree does not need a single substitute for a dropped
column; it can rebuild the same function from several weakly-correlated ones. So
pairwise correlation only puts a floor under replaceability.

**What the test did show: the two models prefer different *kinds* of column, and
the preference is algorithmic.**

| kind of column | XGBoost | LightGBM |
|---|---:|---:|
| target: position/letter indicators (binary) | **49%** | 18% |
| target: quantum descriptors (continuous) | 19% | **33%** |
| flanking-DNA windows (continuous) | 9% | **25%** |
| flanking DNA, nearest 10 letters | 11% | 17% |
| other published columns | 12% | 6% |

Almost the whole disagreement is one axis: **XGBoost leans on binary indicators,
LightGBM on continuous columns** — consistent with LightGBM's histogram binning
and leaf-wise growth making continuous features cheap to split on repeatedly. So
a large part of what an importance ranking reflects is **the splitting algorithm's
affinity for a column's data type**, not how much the column matters.

> ⚠️ **This axis is an artefact of the measurement, not of the models.** The
> table above uses split gain. Measured with SHAP instead, XGBoost's share on the
> binary indicators falls from 57% to 31% and lands on LightGBM's 31% exactly —
> see "Does this generalise beyond two boosting libraries?" below, which also
> narrows the conclusion drawn from it.

That also explains the −0.02. Gain credits whichever interchangeable column a
tree used first; permutation asks what breaks when that column is destroyed, and
answers "barely" if the others rebuild it. The two diverge most where importance
is spread thinly across many interchangeable columns — XGBoost's situation
exactly (330 columns, mostly binary slivers) against LightGBM's (151, continuous,
individually harder to replace).

**The decisive check says neither has found the true features.** Dropping each
model's top 20 and refitting costs barely more than dropping 20 at random (+0.008
and +0.003). **Method agreement measures how concentrated and irreplaceable an
attribution is, not whether it is right.** LightGBM's ranking is *reproducible*,
which is the reason to use it; reproducible is not correct.

### Does this generalise beyond two boosting libraries? Partly — and SHAP changes the answer

The section above rests on two models of the same family and on one way of
measuring importance. The previous version of this report flagged both gaps and
said the test had not been run. It has now been
(`src/sgrna/importance_models.py`, `results/importance_model_*.csv`): seven model
families — XGBoost, LightGBM, CatBoost, random forest, extra trees, ridge,
elastic net — and three ways of measuring importance, on identical folds. Pairs
are matched so the three methods are compared on exactly the same model pairs.

| how importance is measured | agreement between models, per column | per block of related columns |
|---|---:|---:|
| what the model reports by default | **0.050** | 0.531 |
| destroy one column and see what breaks | 0.258 | 0.581 |
| **SHAP** | **0.496** | **0.800** |

**The first finding holds and widens.** Default importance agrees between model
families at 0.050 per column against 0.531 per block, including across bagged
forests, which fit independently rather than on residuals, and penalised linear
models, which have no selection step at all. So the column-level instability is a
property of this matrix, not of boosting.

**The second finding reverses a prediction we recorded in advance.** The
expectation was that SHAP would fail the same way: Shapley values divide credit
for one prediction among its inputs, so two columns carrying identical
information should split it. **That is not what happens.** XGBoost and LightGBM
rank each other's columns at **−0.07 by default importance and +0.87 by SHAP**
(per block, 0.59 → 0.98). And the data-type axis above collapses:

| share of importance on the binary indicators | default | SHAP |
|---|---:|---:|
| XGBoost | **57%** | **31%** |
| LightGBM | 30% | 31% |
| CatBoost | 27% | 25% |
| random forest | 16% | 24% |
| extra trees | 46% | 47% |

Default gain counts how much a split improved the fit at the moment it was made,
which over-credits a binary column used in many shallow nodes; SHAP measures the
effect a column has on the output. These seven models fit nearly the same
function — they span 0.08 Spearman — so once the *function* is measured rather
than the *fitting procedure*, they largely agree. Under SHAP all five tree models
pick out the same two leading blocks: the quantum columns of the nearest 10
flanking letters (≈14%) and the guide's letter-pair quantum columns (≈13%), while
the 3,383 middle-position indicators carry ≈10% between them.

**So the conclusion narrows rather than disappears.** "Do not read mechanism off
an importance plot" should be "do not read it off a **gain** plot". SHAP's 0.496
per column is substantial agreement, not unanimity, so a claim still belongs at
block level, and the three kinds of evidence below are still what the biology
here rests on. But a SHAP ranking is reproducible across model classes in a way a
gain ranking is not, which makes it a usable result rather than only a caveat.

### What to trust instead

Three kinds of evidence survived every model, and they are what the biological
claims here rest on:

1. **Group-level ablations defined by hypothesis before fitting**, against the
   family's own shuffled control.
2. **The sign and size of a single named quantity** — target GC −0.201, flank GC
   +0.159. Anyone can recompute these in one line; no model involved.
3. **Transfer** — does the relationship hold in a different screen, enzyme or
   organism? A spurious attribution does not survive it.

None of the biology in Part 5 came from an importance ranking. That should be
stated as policy.

### A consequence for the paper this project extends

Noshay et al. read their biological conclusion — quantum-chemical properties at
the 3′ end of the guide — off the importance ranking of an iterative random
forest on this exact matrix. The allocation of importance *between kinds of
column* here swings from 49% to 18% depending on the algorithm.

This does **not** show their conclusion is wrong; LightGBM puts its largest share
(33%) on the quantum descriptors, which is consistent with it, and under SHAP all
five tree models agree that quantum columns lead. It shows that **a gain-based
importance ranking on this matrix cannot establish the claim on its own** — the
49%-to-18% swing is largely an artefact of how gain is computed, and a method
without that artefact was available. That is a precise, constructive criticism,
and it now comes with the alternative attached.

---

## Part 10 — How much room is left

### The ceiling, recomputed

The standard way to bound a model is the **attenuation argument**: if two
measurements of the same quantity each equal a shared signal plus *independent
noise*, their correlation is the reliability, and no model can correlate better
than its square root with a single observed measurement. Two screens of the same
guide library agree at **ρ 0.8095**, giving √0.8095 ≈ **0.90**.

**Everything rests on "independent noise", and we tested it.** Published read
counts would settle it directly, and they do not exist — the underlying data is
raw reads in SRA (PRJNA450978 and others), so replicate-level scores would mean
re-running the authors' counting pipeline. So instead: **systematic effects are
predictable, noise is not.** Fit a model to the *difference* between the two
screens' scores for the same guide.

| | ρ |
|---|---:|
| predicting the cut score itself (control) | 0.708 |
| **predicting the disagreement between the two screens** | **0.546** |

**The disagreement is highly predictable.** That is not noise — it is a
systematic, sequence-dependent difference between what WT-SpCas9 and eSpCas9 do,
which is exactly what you would expect from two enzymes engineered to differ in
specificity.

| assumption about the 0.19 disagreement | reliability | ceiling |
|---|---:|---:|
| all of it is noise (the original assumption) | 0.810 | **0.90** |
| the predictable ~30% is biology | 0.866 | **0.93** |

So the usable figure is a **range, 0.90–0.93**. But the more important conclusion
is that **these two screens cannot establish a ceiling properly**: they differ in
a way that is neither shared signal nor independent noise, and the shared library
pulls the estimate the other way again. A trustworthy ceiling needs true
replicates — the same library, the same enzyme, two independent experiments —
which nobody has published for a bacterial Cas9 cutting screen.

**SLICER at 0.707–0.718 is therefore somewhere around 80% of the way to a limit
we can only bracket.** That is enough to say neither model is near it, and not
enough to quote a precise remaining headroom.

### Can the ceiling be raised?

It is a property of the *measurement*, so raising it means measuring differently.

1. **Average over replicates.** The ceiling applies to predicting a *single* noisy
   observation. Predicting the mean of three independent screens is an easier
   target and the apparent ceiling rises. Cheapest route, no new technique.
2. **Measure cutting rather than survival.** The current label is depletion from
   a growing population, mixing cutting with repair and growth. A direct readout
   would remove whole categories of noise.
3. **Publish counts.** Had the original screens released per-replicate counts,
   none of the above reasoning would have been necessary.

### Why the learning curve saturates

Four times the training data buys **+0.008 ρ**, flat between 20,000 and 26,000
guides. Three reasons, and together they make it expected rather than surprising:

- **Precision has run out.** The signal is largely additive, and estimating a few
  hundred additive coefficients from 13,000 examples is already comfortable.
- **Complexity cannot be bought.** More rows support a richer function only if
  the extra structure exists in the features — and nothing correlates with the
  model's errors above 0.035.
- **Label noise does not shrink.** More rows average out noise in the *fitted
  parameters*, not in the test labels you are scored against. That part of the gap
  is fixed by the assay.

### Where the remaining room is not

Each was a live hypothesis that got closed: not more rows (+0.008), not model
class (0.116 across sixteen), not a neural network on raw sequence (ties on the
guide, loses on the flanks), not more features (+0.001), not DNA shape (+0.001),
not chromosome position (collapses into the flanks), and not the label being
dirty (cleaning reveals the flank effect rather than raising the ceiling).

### Where everything here stops: human cells

The source paper also published a human dataset sharing 6,216 columns with the
bacterial one, so a model crosses with **no change of representation** — a failure
cannot be blamed on mismatched features.

| | ρ |
|---|---:|
| *E. coli* model on *E. coli* | 0.531 |
| human model on human | 0.404 |
| ***E. coli* model on human** | **−0.048** |
| **human model on *E. coli*** | **−0.017** |

**Useless in both directions**, slightly worse than guessing. And the mechanism
inverts, which is why it is negative rather than merely weak:

| | *E. coli* | human |
|---|---:|---:|
| GC content → cut score | **−0.201** | **+0.017** |
| melting temperature → cut score | −0.201 | +0.017 |

**The GC penalty does not generalise to humans.** In *E. coli* it is the project's
strongest single mechanism; in human data it is inert.

**Do histones explain this?** Partly, and the literature supports the mechanism
though we have not tested it. **Nucleosomes physically block Cas9**, and bacteria
have no histones — so in human cells a major determinant is whether the target is
*accessible*, a variable absent from our data and not inferable from sequence.
The published human GC effect is also **non-monotonic** rather than absent: very
high and very low GC both work less well, with roughly 40–60% the useful range,
and a *linear* association is not significant. Our +0.017 is a linear
correlation, so it is consistent with a real U-shape being invisible to the
measure used. And accessibility dominates in a way it cannot in bacteria: targets
in promoter regions, which are kept open, cut better than intergenic ones.

Fair summary: **in human cells accessibility is a first-order determinant and
target GC at best a weak non-linear one; in bacteria there are no nucleosomes and
GC is a strong monotonic one.** Bacteria are not chromatin-free — HU and H-NS,
which is what `e_nucleoid` was about — but that packaging is not nucleosomal.

---

## Part 11 — What is new, what is imported, what is open

### Imported, and openly so

`b_energy` is the CRISPRoff authors' energy model; `c_folding` is ViennaRNA;
`f_methylation` a text search; `h_offtarget` a genome scan; `i_shape` uses
published dinucleotide scales; the position sets read published GapR-seq, Hi-C and
RegulonDB/PRECISE-1K data. **Using them is not a contribution.** What is, is that
**seven were measured against controls and found unable to help**, with Part 6's
rule explaining why in a way that generalises.

### New

1. **The published dataset decoded and made portable** — 99% of 6,232 columns
   regenerate from any 20-letter guide, verified at 100% on 1.17 million values.
   Everything else depends on it.
2. **A rule for when a feature cannot help, with a seconds-long test**, applied as
   a forecast and correct 8 times out of 8 once stated correctly — with the
   sharpening of its statement being part of the result.
3. **A control that overturned our own positive result**, of a kind not standard
   in this field; and a second effect traced to a measurement artefact.
4. **Interpretability measured and then diagnosed** — the disagreement between
   importance methods is a consequence of redundancy plus each algorithm's
   preference for binary or continuous columns, which means an importance ranking
   on this matrix cannot carry a biological claim, including the inherited
   paper's.
5. **The flank effect characterised** — a gradient rather than a motif, peaking at
   250–500 letters, downstream-weighted, transferring across enzyme and organism
   at ~90% and across kingdoms at 0%.
6. **The ceiling argument tested rather than assumed**, and found to rest on a
   false independence assumption — with the systematic difference between two
   Cas9 variants quantified as a by-product.
7. **Measured boundaries**: where the method stops, what the data cannot answer,
   and why.

### Open

- **Cross-species.** Both new screens are *Enterobacteriaceae*, and the
  *C. rodentium* one covers 4.3% of one chromosome. **No genome-wide Cas9-cutting
  screen exists outside this family** — what exists in *Bacillus*,
  *Mycobacterium* and so on is CRISPRi, a different quantity. The competition is
  in the same position: their two non-*E. coli* validations are a 236 kb fragment
  and 296 guides on a plasmid inside *E. coli*. Closing this needs wet-lab work.
- **The long-range mechanism.** We know the effect is compositional and
  downstream-weighted; we do not know what it physically is.
- **A trustworthy ceiling**, which needs replicates nobody has published.

---

## Part 12 — Next steps and how to frame the write-up

### Next, in order of value for effort

1. **Tune both models, or publish the parity result.** Our tuning gain is the size
   of the whole margin, so either crisprHAL 2 gets the same search — about a day
   of CPU — or the paper claims parity and says why. *Parity is a perfectly good
   result; an unsupported lead is not.*
2. **Test whether the upstream/downstream asymmetry is a transcription effect.**
   Repeat the side-by-side split while conditioning on whether the protospacer
   lies on the coding or template strand of its gene. Cheap, and it would convert
   a solid observation with a partial explanation into a mechanism.
3. **Put error bars on the cross-organism transfers.** Each cell of the 2×2 is a
   single fit. Several seeds plus a bootstrap over the test set is an afternoon,
   and it is currently the weakest-supported table in the report.
4. **Confirm the human GC relationship is U-shaped** in the human matrix we
   already hold. Cheap, and it upgrades "the mechanism does not transfer" to "the
   mechanism is replaced by a different one".
5. **Find or generate a genome-wide screen in a distant bacterium.** The one open
   question analysis cannot close.
6. **Adopt SHAP as the reported importance method.** It agrees across seven model
   families at 0.50 per column and 0.80 per block where gain agrees at 0.05, and
   it removes a 26-point artefact from XGBoost's attribution (Part 9). The code
   exists; what remains is regenerating the figures that currently show gain.

### A note on tooling, deliberately parked

A guide-design tool or web UI is **not** a deliverable of this project. The
findings are the contribution, and the field already has several predictors, so a
tool without the findings would add little. Everything needed to build one later
exists — `featurise.matrix()` turns any 20-mer into the published representation,
and the model runs in under a second — but it is a separate piece of work and
should not compete with finishing the paper.

### How to frame it

**Not** "a better guide predictor". At parity that framing invites the one
comparison this project loses — against a full-time lab with a GPU.

Frame it as: **what determines whether a CRISPR guide works in *E. coli*, and how
to tell in advance whether a proposed explanation can possibly help.** Then every
result is load-bearing:

- a published dataset decoded, verified and made usable by anyone;
- a rule predicting which feature ideas cannot work, applied as a forecast and
  correct 8 of 8, with its statement sharpened by the hardest case;
- our own positive result overturned by the source experiment's own control, and a
  second traced to a sequencing artefact;
- the first *measured and then diagnosed* interpretability analysis in this
  literature, with a specific consequence for the inherited paper — including
  which importance method survives a change of model family and which does not;
- a 6,232-column published feature set reduced to 427 columns with no measurable
  loss, and decomposed into sequence, non-sequence and restatement;
- the flank effect characterised as a gradient and shown to transfer across enzyme
  and organism at ~90%, across kingdoms at 0%;
- a ceiling argument tested and found to rest on a false assumption;
- and, as a by-product, parity with the state of the art at about 1/74th of the
  CPU time.

### Honest caveats to carry into the paper

- **Parity, not a win.** The +0.0107 margin is inside the tuning effect.
- **The speed advantage is CPU-to-CPU.** crisprHAL 2 as published is GPU-trained.
- **Cross-species is untested**, and the competition's position is no better,
  which does not make ours good.
- **Cross-enzyme evidence is narrow** — three point mutations on an identical
  library.
- **The ceiling is a bracket, not a number**, and the assumption behind the usual
  derivation is false here.
- **Everything is bacterial.** Measured, not hedged: 0% transfer to human cells.
- **`h_offtarget` is a null result** and `f_methylation` marginal, by their own
  floors.
- **Cross-organism transfer figures have no error bars yet.**
- **The downstream/upstream asymmetry has a mechanism for its short-range half
  only.**

---

## Appendix A — Every feature set, grouped, with its data source

Nine sets, **four distinct ideas**. This is how they should appear in the paper.

### Group 1 — the guide's own sequence (redundant by Part 6's rule)

| set | what it measures | how computed | data source | Δρ |
|---|---|---|---|---:|
| `b_energy` | guide–DNA binding energy, by component | the CRISPRoff energy model run on the guide | [CRISPRoff](https://github.com/RTH-tools/crisproff) repo, imported unchanged | +0.007 |
| `c_folding` | whether the guide RNA folds on itself | ViennaRNA folding of spacer, and spacer+scaffold | ViennaRNA library | +0.007 |
| `d_mechanics` | how easily the duplex opens, position by position | nearest-neighbour thermodynamics | SantaLucia & Hicks (2004) parameter tables | +0.020 |
| `f_methylation` | Dam/Dcm motifs at the PAM and seed | text search for `GATC`, `CCWGG` | reference genome | +0.003 |

### Group 2 — the surrounding DNA (the one that worked)

| set | what it measures | how computed | data source | Δρ |
|---|---|---|---|---:|
| `a_flank` | flank composition at four scales | GC, purine fraction, longest letter run and averages of the recovered quantum tables over windows of 50/250/500/1000 letters each side; plus letter identity at the 10 nearest positions | reference genome + the quantum tables recovered in Part 4 | **+0.080** |
| `i_shape` | physical bendability and stiffness of the flanks | nine dinucleotide-step scales averaged over windows, plus a phased bend sum at the 10.5-letter helical repeat | published dinucleotide shape scales — two needed correcting: one scale was corrupt (excluded), two were not strand-symmetric as published (symmetrised) | +0.036 alone, **+0.001** on top of `a_flank` |

### Group 3 — chromosome position (one variable; collapses into Group 2)

| set | what it measures | how computed | data source | Δρ |
|---|---|---|---|---:|
| `d_supercoiling` | DNA twisting, nominally | read density of GapR ChIP tracks in windows, plus ratios against controls | GEO **GSE152880**, including its untagged and rifampicin control arms | +0.045, of which +0.004 is twisting-specific |
| `e_nucleoid` | 3D crowding, packaging-protein dependence | contact counts from a Hi-C matrix | Lioy et al. 2018 Hi-C matrices | +0.016 |
| `g_transcription` | distance to a promoter, strand, expression | promoters re-located by matching the 80-letter sequence each entry ships with, avoiding a coordinate-system mismatch | RegulonDB + PRECISE-1K | +0.015 |
| mappability+GC | the sequencing artefact behind all of the above | 25-letter uniqueness and GC fraction at four windows — **the 8 columns of Part 7** | reference genome only | +0.038 alone, **−0.000** on top of `a_flank` |

### Group 4 — the rest of the genome

| set | what it measures | how computed | data source | Δρ |
|---|---|---|---|---:|
| `h_offtarget` | copy number and near-match burden | genome-wide scan for sequences within a few letters of the target | reference genome | +0.000 — **a null result by its own floor** |

**Labels:** the published `cut.score` from Guo et al.'s 2018 *E. coli* depletion
screen (13,880 guides, via Noshay et al.'s matrix), and the read-count-filtered
re-derivation shipped with crisprHAL (33,567 guides). **Genomes:** *E. coli*
NC_000913.2 and *C. rodentium* ICC168 NC_013716.1 (5,346,659 letters, 54.7% GC).

---

## Appendix B — Papers referred to

- [Noshay et al. — *Quantum biological insights into CRISPR-Cas9 sgRNA efficiency*, NAR 2023](https://academic.oup.com/nar/article/51/19/10147/7279034) — the matrix this project decodes
- crisprHAL 2 — *Better data for better predictions*, PeerJ 2026 · [PeerJ](https://peerj.com/articles/20706/) · [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12903899/)
- [crisprHAL 1 — *A generalizable Cas9/sgRNA prediction model*, Nat Commun 2023](https://www.nature.com/articles/s41467-023-41143-7) — the cross-species claims belong here
- DeepCC9 — *An interpretable deep learning framework*, Bioinformatics 2026 · [Oxford Academic](https://academic.oup.com/bioinformatics/article/42/7/btag483/8723703) · [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13384063/)
- [Liu et al. — *Sequence features associated with cleavage efficiency*, Sci Rep 2016](https://www.nature.com/articles/srep19675) — the human GC U-shape and chromatin accessibility
- [*Nucleosomes impede Cas9 access to DNA*, eLife 2016](https://elifesciences.org/articles/12677) and [*Nucleosomes inhibit Cas9 cleavage in vivo*, PNAS 2018](https://www.pnas.org/content/115/38/9351)
- [*Improved prediction of bacterial CRISPRi guide efficiency*, Genome Biology 2023](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03153-y) — the closest adjacent problem
- [FDA — approval of the first CRISPR therapy](https://www.fda.gov/news-events/press-announcements/fda-approves-first-gene-therapies-treat-patients-sickle-cell-disease)
