# SLICER — what determines whether a CRISPR guide cuts, in *E. coli*

**S**imple **L**ightweight **I**nterpretable **C**RISPR **E**fficiency **R**anker.

RSI09. State of the project as of 4 October 2026. Each result says what it means
before giving the figure; Parts 1 and 2 cover only the points where a reader who
knows the background still gets *this* dataset wrong.

This document states what we currently think and the evidence for it. It
deliberately contains no record of what we used to think — sixteen claims have
been withdrawn along the way, and those live in a separate working file
(`local/MISCONCEPTIONS.md`, not in the repo) so that this one reads as a
statement rather than a diary. `FINDINGS_LOG.md` is the chronological record and
`../results/RESULTS.md` has the full protocol behind every number.

---

## Contents

- [Part 0 — Why this problem is worth working on](#part-0--why-this-problem-is-worth-working-on)
- [Part 1 — Six things about this problem that are easy to get wrong](#part-1--six-things-about-this-problem-that-are-easy-to-get-wrong)
- [Part 2 — What the headline number does and does not mean](#part-2--what-the-headline-number-does-and-does-not-mean)
- [Part 3 — What SLICER is, and where it stands](#part-3--what-slicer-is-and-where-it-stands)
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

## Part 1 — Six things about this problem that are easy to get wrong

Not a glossary. These are the points where a reader who knows the general
background still reaches the wrong conclusion about *this* dataset, and each one
is load-bearing somewhere below.

**The screen measures survival, not cutting.** Cutting the chromosome usually
kills a bacterium, so the experiment grows a mixed population, sequences what
survives, and scores a guide by how far its cells *disappeared*. Survival mixes
cutting efficiency with DNA repair, growth rate and sequencing noise. Every
"efficiency" number in this report is that composite, and Part 7 contains a
result we lost to forgetting it.

**The flanking DNA is not supposed to matter.** The guide matches 20 letters;
Cas9 reads those plus a three-letter PAM immediately downstream of them. It does
not read the surrounding DNA. The central finding here is that the surrounding
DNA predicts cutting anyway, which is why it took controls rather than an
ablation to believe.

**"Downstream" means PAM-side, not gene-side.** Upstream and downstream are
defined along the strand the guide matches, so the PAM is downstream by
construction. The asymmetry in Part 5 is about the PAM side, and Part 5 also
shows why this screen cannot tell that apart from the gene's own direction.

**GC content is the mechanism, not a covariate.** G–C pairs are held by three
hydrogen bonds and A–T by two, so GC-rich DNA is harder to pull apart. Cas9 has
to open the duplex to engage, which makes "how GC-rich is the DNA around here" a
physical cause rather than a proxy, and explains why four separate feature
families turned out to be one variable (Part 7).

**Bacteria have no nucleosomes.** Humans, animals and plants wrap DNA around
histone spools that physically block Cas9, and in human screens chromatin
accessibility is the dominant predictor. Bacteria pack DNA with proteins like HU
and H-NS instead, which is a different kind of obstacle. This is the reason the
cross-kingdom transfer in Part 10 is zero rather than merely weak.

**~20 guides sit in every gene, so plain cross-validation leaks.** A model told
roughly *where* a guide is can score well by recognising a locus. Anything
positional in this report is measured under grouped cross-validation — hide a
whole contiguous stretch of chromosome — and where that changes the answer, it
is said so.

One dataset caveat that belongs with these: **CRISPRi screens cannot substitute
for cutting screens.** CRISPRi uses a disabled Cas9 that binds without cutting,
so it measures a different quantity. Most bacterial screens outside *E. coli*
are CRISPRi, which is why Part 5 ends with an open question rather than a
cross-species result.

---

## Part 2 — What the headline number does and does not mean

Three facts about the metrics that changed a conclusion in this project, and an
audit of which numbers carry error bars.

### A difference in ρ is not a difference in an amount of anything

Spearman ρ measures whether the ordering is right, which is the right headline
because nobody needs a guide's exact efficiency — they have candidates and want
the best one. But **ratios of Δρ are not ratios of information.** On the ρ²
scale, 0.53 → 0.61 is worth +0.09 while 0.82 → 0.90 is worth +0.18: the same Δρ
is worth twice as much higher up. So when this report says the flank features
extract 4.2× what a CNN gets from the same DNA, that compares two measurements
taken **at the same baseline** and does not generalise to other baselines.

### What ρ buys on the decision a practitioner actually makes

Simulating it directly — take ten candidate guides, let a model of given ρ
choose, record where its pick really falls:

| model ρ | pick lands at |
|---:|---|
| 0 (blind) | 50th percentile |
| 0.53 | 73rd |
| 0.61 | 77th |
| **0.71** | **81st** |
| 0.90 | 88th |
| 1.00 | 91st |

**The whole range from useless to perfect spans the 50th to the 91st
percentile** — a perfect model does not hand you a perfect guide, because ten
random candidates may not contain a great one. And **most of the benefit arrives
early**: the first 0.53 of ρ buys 23 percentile points, the next 0.47 buys 18.
On this scale we are 82% of the way to the ceiling, against 79% on the ρ scale;
the agreement is a coincidence at these values, not a rule.

### Report ρ and R² together, because they can disagree violently

Our own clearest case: a LightGBM with a robust loss scored **R² 0.192 and
ρ 0.568** — same model, same data, broken on one metric and respectable on the
other, because the robust loss shrinks predictions toward the middle and
destroys the scale while preserving the order. For choosing a guide ρ is the
relevant one; for claiming you can predict efficiency, R² is.

### Which numbers have error bars

Models are compared on the *same* folds, so comparisons are paired: for the
head-to-head the paired standard error is **0.0011** against a single fold's
**0.0122**, eleven times smaller, which is the whole reason paired designs are
used. With five folds the Wilcoxon signed-rank test cannot go below p = 0.0625
however large the effect, so both it and the paired *t* are quoted — the *t* has
the power, the sign test bounds what a rank-based reading can claim.

| kind of number | replication | error estimate |
|---|---|---|
| headline ablations and head-to-head | 5–25 folds, 1–5 seeds | **yes** — fold sd, paired tests |
| transfer within a screen | 5 folds | **yes** — fold sd and paired p |
| cross-organism model transfer | 3 model seeds, one test set | **partial** — no interval over the test set |
| correlations with the label | n = 13,880–59,489 | implicit, SE ≈ 0.004–0.008 |
| redundancy and cross-redundancy R² | 5-fold out-of-fold | a distribution over columns, not an interval per column |
| the ceiling | — | **a range**, Part 10 |
| timings and peak memory | 5 repetitions, pinned threads | **yes** — IQR |

**And statistical significance is not importance.** The head-to-head margin is
highly significant at every stage — +0.0111 untuned, +0.0161 with both models
tuned, 5 of 5 folds either way — and it is still small enough that what makes it
quotable is the *control*, not the *p*-value: it only became a claim once both
models had been searched, because before that an effect of the same size was
available to whichever one got tuned.

---

## Part 3 — What SLICER is, and where it stands

### The architecture, and why each piece is what it is

Five stages, and every one of them is a choice that was tested rather than
inherited.

**1. Input — named columns describing the guide and its surroundings.** The
6,232 published columns (Part 4 shows these are a re-spelling of the guide's 20
letters) plus 348 columns describing the flanking DNA in windows of 50 to 1,000
letters. 6,580 columns on the original screen, 6,517 on the curated one. *Why
named features rather than raw sequence:* Part 5 shows the signal is a smooth
compositional gradient, and a windowed mean computes a gradient exactly where a
pattern detector has to approximate it — on the same DNA, hand-computed windows
are worth 4.2× what a convolutional network extracts.

**2. Imputation — column medians, from the training fold only.** Mundane but
load-bearing: medians taken over all rows would leak the validation fold's
distribution into training.

**3. Feature selection — rank by importance, keep the top 300.** *Why select at
all:* most columns are restatements, and a tree spends its depth budget on
whatever is offered. *Why 300:* raising the cap to 1,200 is worth +0.001. *Why
XGBoost's gain to do the ranking,* when LightGBM ranks as well 5× faster: this
is the one place the project prefers continuity over speed, because these are
the settings under which the published 0.2937 / 0.5278 baseline reproduces, and
that reproduction is the check that tells a reader the harness is faithful. The
fast alternative is one config setting away (`config.SELECTOR`) and is used for
sweeps.

**4. The model — gradient-boosted trees: 400 of them, depth 4, learning rate
0.05, 80% row and column subsampling.**

> **Which library, and why there are two answers.** This matters for reading the
> rest of the report, so it is stated once here. The **baseline configuration**
> on the 13,880 published rows uses **XGBoost**, because those are the settings
> under which Noshay et al.'s published figure reproduces, and that reproduction
> is the harness check. The **head-to-head configuration** on the 33,567 curated
> rows uses **LightGBM** as the predictor, with XGBoost still doing the column
> selection. So ρ 0.5278 is XGBoost and ρ 0.7078 is LightGBM, and the tuning
> gain of +0.0131 that exceeds the whole crisprHAL margin is a LightGBM result.
> The two libraries are within 0.001 of each other everywhere they have both
> been run (Part 9), so nothing substantive turns on the choice — but a reader
> tracking a specific number should know which one produced it.

*Why boosting:* of sixteen model classes
on identical data the whole field spans 0.116 ρ, and boosting is at the top of
it. Two of the failures are informative rather than embarrassing — a plain ridge
regression gets within 0.03, so **the effects mostly add up rather than
interacting**, and nearest-neighbours fails outright, so guides with similar
features do *not* have similar scores. That combination is the signature of many
small independent effects, which is what boosting with shallow trees is for.
*Why depth 4:* deeper is monotonically worse at every representation tested
(Part 4).

**5. Evaluation — five-fold cross-validation, grouped wherever a feature is
positional.** The screen puts ~20 guides in every gene, so plain `KFold` lets a
model score well by recognising a locus rather than reading a guide.

**What is deliberately absent.** The default configuration is *not* tuned —
tuning is worth +0.010 and is reported separately, because the untuned
configuration is the one that anchors to the published baseline. The matrix is
*not* stored in its reduced 427-column form, which would cut memory 5× (Part 4),
for the same reason. Both are available; neither is the default while the anchor
is worth more than the saving.

**Cost, so the shape of the thing is clear:** one fold is 17 s of CPU and 3.3 GB,
of which 12 s is the selection step whose only output is a ranking. No GPU, no
accelerator, nothing that is not a laptop.

### Where it stands

**Everything in this table was run here, on the same rows and the same folds.**
No row is a figure copied from a paper, which is why they can be compared at all.

| | what it is | Spearman ρ | CPU per fold | peak RAM |
|---|---|---:|---:|---:|
| **SLICER, tuned** | this project, settings searched to convergence | **0.7212 ± 0.0016** | 17 s | 3.3 GB |
| **SLICER** | this project, default settings | **0.7082** | 17 s | 3.3 GB |
| **crisprHAL 2, tuned** | the best published bacterial model, searched | **0.7051 ± 0.0092** | 69 min | 4.6 GB |
| **crisprHAL 2** | as its authors published it | **0.6971 ± 0.0067** | 21 min | 4.6 GB |
| our CNN on raw sequence | the gradient-vs-motif test, Part 5 | 0.5165 | 2.3 min | — |
| Noshay et al.'s iRF | the model whose features SLICER inherits | 0.4785 | 6.2 min | 2.2 GB |
| apparent ceiling | Part 10 | ≈0.94 (range 0.94–0.97) | — | — |

Four things that table says at a glance. **The two serious models are within
0.016 of each other** and everything else is far behind. **The gap between
tuned and untuned (0.008–0.013) is comparable to the gap between the models**,
which is why tuning both was necessary before any comparison could be believed.
**Replacing Noshay's iRF with gradient boosting on their own matrix is worth
+0.049** — more than anything else in the table except the flanking-DNA features
themselves. And **the cost column spans three orders of magnitude** at almost
constant accuracy between the top two.

Two honest caveats on the cost column: it is **CPU-to-CPU**, and crisprHAL 2 as
published is GPU-trained, so this is a statement about their architecture's CPU
cost and not about beating a GPU. And their tuned configuration is wider and runs
longer than their default, which is why tuning moves their cost from 21 to 69
minutes while ours does not move at all.

### Their model, actually run

The first row used to be a Pearson correlation copied from their Table 1, which
is not comparable with the rest of the column. Their model has now been
reimplemented from the methods section and run here — iterative Random Forest,
1,000 trees per forest, ten iterations, five-fold cross-validation on the full
6,232-column matrix (`src/sgrna/noshay_replicate.py`):

| | Pearson | R² | Spearman |
|---|---:|---:|---:|
| their Table 1 | 0.5019 | 0.2491 | not reported |
| **our run of their algorithm** | **0.4881 ± 0.0119** | **0.2345 ± 0.0102** | **0.4785** |

**Within 0.014 Pearson of their published figure on a third of their rows** —
their Table 1 used 40,468 sgRNAs, and the supplementary matrix they released
covers 13,880, so some shortfall was expected. That is close enough to treat the
reimplementation as faithful, and it supplies the Spearman their paper never
reported, which is what makes the comparison column above honest.

**The like-for-like model comparison this allows is the useful part.** On the
*same* 13,880 rows, the *same* 6,232 columns and the same folds, their iRF
reaches ρ 0.479 and our gradient-boosted baseline reaches **0.527** — so
**+0.049 from changing the model alone**, before any new feature does anything.
The flank family then adds +0.080 on top. Worth stating in that order, because it
separates what we contributed from what a more modern learner contributes for
free. (One asymmetry: iRF uses all 6,232 columns by design, where our pipeline
selects 300. That is part of what is being compared, not a confound to remove.)

**And the cost gap is the same story as crisprHAL's.** Their iRF takes **6.2 min
per fold** on ten cores against our 17 s, at 2.2 GB against 3.3 GB — 22× the
time for 0.049 less ρ. Their paper ran on Oak Ridge's CADES cluster; the point is
not that they were wasteful but that the iterative reweighting buys nothing here,
which is consistent with Part 4: there is no deep interaction structure to
amplify, because the matrix is a re-spelling of 20 letters.

### The head-to-head, with both models tuned

Earlier numbers in this project were measured on different guides from the model
they were compared against. That was fixed in two steps. Because of Part 4 the
published feature set can be computed for *anyone's* guides, so SLICER was run on
crisprHAL 2's own 33,567 curated guides with their label. Then rather than trust
their published figure, their actual model — their architecture, their
hyperparameters, 48 epochs, code imported unchanged — was re-run on our exact
folds. As a check that the harness is faithful rather than flattering, it was
also run on their own shipped split and scored **0.6940** against their published
0.695.

That left one objection, and it was the right one: **neither model was tuned, so
the margin might only reflect whose defaults suited this data.** Both have now
been searched, each inside the training folds so no configuration ever sees a
validation fold.

| | default | tuned | tuning gain | folds improved |
|---|---:|---:|---:|---|
| crisprHAL 2 | 0.6971 ± 0.0067 | **0.7051 ± 0.0092** | **+0.0080** | 5/5, p = 0.018 |
| **SLICER** | 0.7082 | **0.7212 ± 0.0016** | **+0.0131** | 5/5, p = 0.0001 |

| comparison, paired by fold | mean | folds | paired *t* |
|---|---:|---:|---:|
| untuned | +0.0111 | 5/5 | p = 0.0002 |
| **both tuned** | **+0.0161** | **5/5** | **p = 0.0003** |

**Tuning both sides widens the margin rather than closing it.** The objection is
answered: the lead is not an artefact of default settings. Both models improve
under search, and ours improves more.

**What the claim can now be, and what still limits it.** The previous honest
statement was *parity*. It can now be **a small, consistent lead: +0.016 on
identical rows, labels and folds, in 5 of 5 folds** — with one asymmetry stated
openly, because it is large. SLICER's search ran **627 configurations** across
the five folds; crisprHAL's ran **10**, on a single inner split, because one of
its fits costs about a thousand times one of ours. Theirs is much the weaker
search, so +0.016 is an upper bound on the true tuned gap.

**How much that asymmetry can be hiding is itself measurable.** SLICER's search
trace records best-so-far at every draw, so we can ask what a 10-draw budget
would have found for *us*:

| draws | share of the total tuning gain captured |
|---:|---:|
| 5 | 47% |
| **10** | **66%** |
| 20 | 74% |
| 40 | 100% |

A ten-draw search recovers about two thirds of what a converged one finds. If
crisprHAL's ten draws behaved similarly, its converged gain would be near
**+0.012** and its tuned score near **0.709** — leaving SLICER ahead by about
**0.012** rather than 0.016. To erase the margin entirely their tuning would have
to find **+0.024**, three times what it found and well outside that pattern.

So: **the margin is small, it survives tuning both sides, and the correction for
their smaller search shrinks it by about a quarter without removing it.** That is
the strongest statement the evidence supports, and it should be written in exactly
those terms rather than as a win.

### What they do differ in: cost

Measured on one machine, CPU to CPU, threads pinned, five repetitions:

| | SLICER | crisprHAL 2 |
|---|---|---|
| one fold | **17.1 s** (IQR 0.7) | ~21 min |
| five folds | **1.4 min** | **105 min** |
| peak memory, one fold | 3.33 GB (3.23–3.68) | 4.55 GB |
| features | named quantities you can look up | learned, inside a network |

**About 74× less time for equal accuracy, at broadly similar memory.** The saving
is time, not footprint — SLICER holds a 33,567 × 6,517 table in memory and is not
light on RAM.

> **Where the footprint actually comes from — measured, and one lever worth
> having** (`src/sgrna/bench.py`, `results/memory_comparison*.csv`). Peak RSS,
> threads pinned, every repetition in a fresh process, three repetitions each,
> because a single shot is not reliable here:
>
> | configuration | columns | peak (median, range) | selection | ρ |
> |---|---:|---:|---:|---:|
> | imports only — the floor | — | 0.17 GB | — | — |
> | all columns + XGBoost (current) | 6,517 | **3.33** (3.23–3.68) | 12.0 s | 0.6939 |
> | all columns + LightGBM | 6,517 | 3.05 (3.05–3.05) | 2.4 s | 0.6945 |
> | reduced columns, subset after loading | 712 | 2.83 | 1.4 s | 0.6823 |
> | **reduced columns, stored reduced** | 712 | **0.65** | 1.4 s | 0.6823 |
>
> **The model is not the lever.** The 0.3 GB between the two selectors is the
> same size as XGBoost's own run-to-run spread, so the two are indistinguishable
> on memory; what LightGBM buys is 5× on selection time, which is a separate
> argument (Part 9).
>
> **Subsetting columns after loading is not the lever either** — 14%, because
> peak RSS is already set by the time the subset is taken. **Storing the reduced
> matrix is**: reading a 712-column store instead of a 6,517-column one is a
> **5× reduction, 3.33 GB → 0.65 GB**, for −0.0116 ρ, which is the same cost the
> reduced representation shows in Part 4. Not adopted as the default, because the
> full-matrix form is what the reproduction check above anchors to; available if
> footprint ever becomes the binding constraint.
>
> **Protocol, because it decides the answer:** peak memory must be measured one
> configuration per process, repeated. `ru_maxrss` is a high-water mark that
> never falls, so two arms sharing a process share a number; and run-to-run
> spread for one arm here is ~0.3 GB, wider than most differences worth
> reporting.

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

**Where the last 0.019 comes from.** It cannot be information, since nothing in
6,232 columns is absent from the 74. The obvious candidate is that the extra
columns pre-compute letter combinations our depth-4 trees are too shallow to
build themselves — which would mean deeper trees close the gap. **They do the
opposite, monotonically, at every encoding:**

| | depth 4 | 6 | 8 | 12 |
|---|---:|---:|---:|---:|
| 74 columns | **0.5081** | 0.4953 | 0.4739 | 0.4557 |
| 364 columns | **0.5189** | 0.5120 | 0.4988 | 0.4858 |
| 1,510 columns | **0.5163** | 0.5081 | 0.4972 | 0.4781 |

Supplying explicit triple-letter indicators hurts too (0.5163 against 0.5189 for
pairs alone). So depth is not the missing ingredient, and the gap is not about
letter combinations.

What accounts for it is the **63 hand-named columns**, which are the one part of
the matrix that is *not* derived from the guide's 20 letters — they include how
far along its gene the guide sits, and letter counts that disagree with a direct
count by up to 13. Added to the 364-column encoding they give **427 columns at
ρ 0.5249, statistically indistinguishable from all 6,232** (paired across five
seeds, mean difference 0.0023, p = 0.11).

So the matrix is **74 columns of sequence, 63 columns of something that is not
sequence, and about 6,095 columns of restatement.** A 427-column replacement
exists; it is not the default only because the 6,232-column form is what the
reproduction check in Part 3 anchors to.

### What the reduced column set is worth in practice

Three measured consequences, all at the same accuracy cost of about 0.01 ρ
(`results/memory_comparison.csv`, `results/net_cost.csv`):

| | all columns | reduced | change |
|---|---:|---:|---|
| peak memory, if the reduced matrix is what gets stored | 3.33 GB | **0.65 GB** | **5× less** |
| feature-selection step | 12.0 s | **1.4 s** | **8× faster** |
| a dense neural net on the feature table: training | 6.1 s | **1.4 s** | **4.4× faster** |
| that net's parameters | 1,685,121 | **199,041** | 8.5× fewer |

The net result is the one worth separating out, because it is easy to
over-generalise. **A network that reads the feature table does get much cheaper**
— almost all of a dense net's parameters live in its first layer, so cutting the
input width cuts the model nearly proportionally, and accuracy is unchanged
(ρ 0.6401 against 0.6413). **But the sequence CNN gets nothing**: `seqnet.py` and
crisprHAL 2 both read raw flanking DNA and never touch the feature table, so
their cost is epochs × convolutions over sequence length. For the architecture
this project actually competes against, a reduced column set is not a speed-up at
all.

The memory saving has the same shape as the column-count saving and is available
only under the same condition: **it has to be the reduced matrix that gets
stored.** Subsetting columns after the full matrix is loaded saves 14%, because
peak memory is set by the load. One result we cannot explain and are not
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
a property of the DNA's local state or of the assay.

**The obvious candidate is transcription, and the evidence is now split.** Genes have
a direction, so "downstream of the PAM" is systematically related to the guide's
orientation inside its gene, and RNA polymerase unwinds DNA ahead of itself. Two
tests follow from that, and the first one cannot be run at all:

- **Conditioning on orientation is impossible on *this matrix*.** Checked two
  independent ways — joining the guide index to the reference gene table, and
  reading `g_transcription`'s separately-derived `template_strand` column — the
  published matrix targets the gene's template strand for **13,825 of 13,879
  guides, with 4 exceptions.** There is no variation to condition on. (It also
  means `eng.txn.gene.template_strand` is a near-constant column, which is part
  of why that whole family could only manage +0.015.) An earlier version of this
  report concluded from that the test was impossible *full stop*. It is not —
  see below.
- **The asymmetry does not scale with transcription.** What the screen *can*
  support is splitting guides by their gene's expression and measuring each
  side's gain within each band (`results/asymmetry_expression.csv`):

| expression band | baseline | upstream gain | downstream gain | downstream, count-matched | asymmetry |
|---|---:|---:|---:|---:|---:|
| low | 0.4733 | +0.0296 | +0.0782 | +0.0790 | 0.0486 |
| mid | 0.4802 | +0.0177 | +0.0635 | +0.0604 | 0.0458 |
| high | 0.4836 | +0.0177 | +0.0656 | +0.0539 | 0.0479 |

**Flat.** If polymerase traffic were driving the asymmetry it should grow with
expression; it does not move (0.0486 / 0.0458 / 0.0479), and on the count-matched
comparison it if anything declines. One seed and ~4,600 guides per band, so this
rules out a gradient rather than proving strict independence — but the
transcription account predicted a gradient and there is none.

The count-matched column is its own control and worth keeping. Family A is not
symmetric — 194 downstream columns against 144 upstream, because the windows and
the PAM-side k-mers are not mirror images — so a third of the downstream
advantage could have been column count. Restricting downstream to a random 144
leaves the asymmetry essentially unchanged (+0.0790 against +0.0782 in the low
band), so **it is the side, not the budget.**

#### The orientation test does run — on guides that are not in the published matrix

The coding sub-library binds one strand because it was **repurposed from a
CRISPRi design**, where binding the non-template strand is the point. Guo's
*intergenic* sub-library was built differently: its Methods say the sgRNAs "were
designed to target either of the two DNA strands in this new library, in contrast
to the sgRNAs in our previous library, which bind only the nontemplate strand."
That library is 10,257 guides — 5,559 across 3,146 promoters and 4,698 across
4,140 RBSs — and it is **not in the published feature matrix**, so the features
had to be rebuilt from the genome.

Jacky did that: 8,954 intergenic guides located on NC_000913.2, flanks cut with
the same recipe family A uses, orientation defined against the nearest gene,
grouped CV by 100 kb bins, arms the same shape as `asymmetry.py --run`
(`Jacky/agent-notes/INTERGENIC_ORIENTATION.md`).

| stratum | n | upstream gain | downstream, count-matched | asymmetry |
|---|---:|---:|---:|---:|
| all intergenic | 8,954 | +0.081 | +0.115 | **0.034** |
| **guide opposite the nearest gene** | 4,361 | +0.061 | +0.108 | **0.047** |
| **guide same way as the nearest gene** | 4,593 | +0.102 | +0.107 | **0.005** |

**The asymmetry depends on orientation.** It is large when the guide faces
against the neighbouring gene and ≈0 when it faces with it. A second, stricter
set (Data S4 high-quality, n = 3,645) shows the same pattern more sharply:
0.064 against −0.015.

That is what a transcription account predicts and **the opposite of what pure
PAM/R-loop geometry predicts**, which should be flat across orientation. So the
candidate is re-opened rather than closed.

**How to hold the two results together.** On the *coding* library, the asymmetry
did not scale with how much the gene is transcribed. On *intergenic* guides, it
does depend on which way the guide faces relative to the local gene. Those are
not contradictory — expression level and orientation are different variables, and
only the second is what a "DNA is asymmetric around a transcription unit" account
actually requires. The simplest reading is that the asymmetry tracks the
*geometry* of the local transcription unit rather than its *traffic*.

Four caveats belong with it in the paper: orientation is defined against the
**nearest gene**, not a curated operon or TSS annotation; the intergenic guides
target promoters and RBSs, a different biological context from coding knockouts;
it is one seed with no bootstrap interval on the orientation *difference* yet; and
it does not distinguish polymerase traffic from any other gene-asymmetric
chromosomal feature.

**Where that leaves it:** short range has a mechanistic account from PAM-first
engagement. Long range has a measured asymmetry, one ruled-out explanation
(replication — no sign flip between replichores, 0.049 against 0.048), and one
live candidate that now has direct supporting evidence rather than none.

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

The search is run **to convergence rather than to a fixed count**: up to 300
draws per fold, stopping when 80 consecutive draws fail to beat the best, with
every draw written to `tuning_trace.csv` so the plateau can be inspected instead
of taken on trust.

| | ρ |
|---|---:|
| default settings | 0.7082 |
| tuned | **0.7212** |
| gain | **+0.0131 ± 0.0016, 5/5 folds, paired p = 0.00005** |

Per fold: +0.0139, +0.0111, +0.0149, +0.0117, +0.0136. And the plateau is real —
each fold found its winner early and then searched 80 more draws for nothing:

| fold | draws run | best found at draw | Δρ |
|---:|---:|---:|---:|
| 1 | 115 | 35 | +0.0139 |
| 2 | 157 | 77 | +0.0111 |
| 3 | 141 | 61 | +0.0149 |
| 4 | 117 | 37 | +0.0117 |
| 5 | 97 | 17 | +0.0136 |

627 configurations in total. An earlier 12-draw version of this search found
+0.0103 (`tuning_trials12.csv`), so **a short search understates the tuning
effect by about a quarter** — which matters for how the head-to-head is read,
because the number being compared against crisprHAL's margin is the converged
one.

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

**And model choice costs little because there is nothing left for a better model
to find.** Of the 6,480 columns the champion is given, not one still correlates
with its out-of-fold errors above 0.031. Part 10 is the measurement, including
why that test can be trusted and whether it survives the 427-column reduction.

### Which model explains itself, and what the disagreement was really about

The two boosters tie on accuracy at 0.609, so the choice between them was made
on explanation quality — and by the library's default importance LightGBM wins
clearly:

| | LightGBM | XGBoost |
|---|---:|---:|
| do its two importance methods agree? | **+0.44** | **−0.02** |
| same features chosen on a different split? | 0.72 | 0.49 |
| features for half the importance | 151 | 330 |

**Read through SHAP instead, that advantage disappears**: stability 0.74 against
0.75, method agreement +0.40 against +0.39, and **43 against 44** columns
carrying half the importance. XGBoost's bad row was the measurement, not the
model — and fixing the measurement is worth far more than swapping the model
(XGBoost's own stability goes 0.49 → 0.75 and its load-bearing columns 330 → 44
without anything about the fit changing).

So **interpretability no longer argues for either library**, which settles the
question the other way round from how it was posed: the baseline configuration
keeps XGBoost because that is what reproduces the published figure (Part 3), and
the head-to-head configuration keeps LightGBM because that is what every number
in it was produced with. Neither choice is now defended on explanation quality,
and neither needs to be.

The rest of this section is about what the original −0.02 meant, because the
answer is a result in its own right.

That −0.02 looks alarming, so we tested what it means. Letting each model choose
its own 300 features **by gain** and comparing the top 50: **14 of 50 shared**,
overall gain rankings correlating +0.69, held-out scores 0.627 and 0.622 — tied.
So they really do choose differently and it really costs nothing.

**It is not that they found different biology.** For each column only one model
picked, the closest counterpart in the other's set has median |correlation|
**0.456**, against **0.069** for randomly chosen columns.

**But substitution is not one-to-one either** — only 22% of unshared picks have a
counterpart above 0.7. A tree does not need a single substitute for a dropped
column; it can rebuild the same function from several weakly-correlated ones. So
pairwise correlation only puts a floor under replaceability.

**And the disagreement is about the measurement, not the models.** Divide each
model's importance by what the column describes, under the method each library
reports by default and then under SHAP (`results/attribution_summary*.csv`):

| kind of column | gain: XGB | gain: LGBM | **SHAP: XGB** | **SHAP: LGBM** |
|---|---:|---:|---:|---:|
| target: position/letter indicators (binary) | **49%** | 18% | **21%** | **21%** |
| target: quantum descriptors (continuous) | 19% | **33%** | 29% | 29% |
| flanking DNA, nearest 10 letters | 11% | 17% | 22% | 22% |
| other published columns | 12% | 6% | 18% | 19% |
| flanking-DNA windows, downstream | 4% | **14%** | 5% | 5% |
| flanking-DNA windows, upstream | 4% | **12%** | 4% | 4% |

Under gain the two models look like they disagree about biology, and the
disagreement is almost entirely one axis — XGBoost leaning on the binary
indicators, LightGBM on the continuous columns, consistent with LightGBM's
histogram binning making continuous features cheap to split on repeatedly.
**Under SHAP the disagreement is gone**: the largest gap on any kind of column
is **0.3 percentage points**, against 31 points under gain. The same
regeneration takes the shared columns in the two models' top 50 from **14 to 41**.

So the "different kinds of column" finding was real about gain and false about
the models. Gain counts how much a split improved the fit at the moment it was
made, which over-rewards a binary column used in many shallow nodes; SHAP
measures the effect a column has on the output. **What an importance ranking
reflects is partly the splitting algorithm's affinity for a column's data type —
and that part is removable.**

That also explains the −0.02. Gain credits whichever interchangeable column a
tree used first; permutation asks what breaks when that column is destroyed, and
answers "barely" if the others rebuild it. The two diverge most where importance
is spread thinly across many interchangeable columns — XGBoost's situation
exactly (330 columns, mostly binary slivers) against LightGBM's (151, continuous,
individually harder to replace).

### The check that looked decisive, and what it was actually measuring

The obvious test of whether a ranking means anything is to remove what it points
at and see whether the model suffers. Done on the full matrix, it says no:
dropping each model's top 20 columns and refitting costs **−0.0024**, against
**+0.0015** for dropping 20 at random — indistinguishable. For a long time this
report drew the strong conclusion from that, namely that neither model had found
"the true features".

**That conclusion was wrong, and the reason is the subject of Part 4.** With
6,232 columns of restatement available, removing the 20 columns a model leaned on
leaves many other ways to say the same thing. The test cannot tell an unfaithful
ranking apart from a redundant matrix. Separating them needs a representation
where columns are *not* interchangeable — which Part 4 now provides. Re-run on
the 775-column reduced set (`faithfulness_groups*.csv`, 3 folds, both boosters,
SHAP-ranked):

| columns dropped and refitted | full, 6,580 cols | reduced, 775 cols |
|---|---:|---:|
| top 20 | −0.0024 | **−0.0150** |
| 20 at random | +0.0015 | −0.0006 |
| top 50 | −0.0102 | **−0.0453** |
| 50 at random | −0.0008 | −0.0015 |
| top 100 | −0.0226 | **−0.0816** |
| 100 at random | +0.0016 | −0.0054 |

**On a representation without restatement, the ranking is load-bearing.** The
random control sits at ≈0 in both, so the top-versus-random gap is the statistic,
and it is 3–6× larger on the reduced set.

**Normalising matters here, and the obvious normalisation is the wrong one.**
Twenty columns is 0.3% of the full matrix and 2.6% of the reduced one, so a
column-count comparison looks unfair. But matching on column *fraction* is worse,
not better: dropping the top 170 of 6,580 — the same 2.6% — costs **−0.027**
against +0.001 at random, but it removes **80% of the model's attributed
importance** rather than 35%, which is a far larger ask than 20 of 775. The
right normalisation is the share of importance removed, and on that scale the two
are closely matched (`faithfulness_shap_share.csv`):

| | share of SHAP in the top 20 | top 50 | top 100 |
|---|---:|---:|---:|
| full, 6,580 columns | 34.6% | 53.1% | 68.6% |
| reduced, 775 columns | 37.5% | 57.5% | 73.8% |

So at an almost identical share of attributed importance removed, the full matrix
loses 0.002 and the reduced one 0.015. **Same ranking, same share of the
explanation taken away, six times the damage — because in one case the
information is still reachable elsewhere and in the other it is not.**

Two things follow. **Faithfulness cannot be measured on a redundant
representation**; a "the top features don't matter" result there is a statement
about the matrix, not about the ranking. And **the rehabilitated reading of
method agreement** is narrower than before: it measures how concentrated and
individually irreplaceable an attribution is, which on this matrix is mostly a
property of the encoding. It still does not measure whether an attribution is
*biologically* right — for that, the three kinds of evidence below are what the
claims rest on.

### Does this generalise beyond two boosting libraries? Partly — and SHAP changes the answer

The section above rests on two models of the same family. The wider test
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

**The second finding is that the method matters more than the model.** XGBoost
and LightGBM rank each other's columns at **−0.07 by default importance and
+0.87 by SHAP** (per block, 0.59 → 0.98), and the data-type axis above
collapses:

| share of importance on the binary indicators | default | SHAP |
|---|---:|---:|
| XGBoost | **57%** | **31%** |
| LightGBM | 30% | 31% |
| CatBoost | 27% | 25% |
| random forest | 16% | 24% |
| extra trees | 46% | 47% |

Two things that table adds beyond the earlier one. **The artefact is largely
XGBoost's**: CatBoost and LightGBM barely move, so this is specific to how
XGBoost's depth-wise growth at these settings accumulates gain on binary
columns. And **extra trees' preference is real** — 46% to 47%, unchanged under
SHAP — so a genuine data-type preference is distinguishable from an artefactual
one. Seven models that span only 0.08 Spearman are fitting nearly the same
function, so once the *function* is measured rather than the *fitting
procedure*, they largely agree: under SHAP all five tree models pick out the
same two leading blocks, the quantum columns of the nearest 10 flanking letters
(≈14%) and the guide's letter-pair quantum columns (≈13%), while the 3,383
middle-position indicators carry ≈10% between them.

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

So on these two screens the figure is a **range, 0.90–0.93** — and the more
important conclusion was that **two different enzymes cannot establish a ceiling
properly**: they differ in a way that is neither shared signal nor independent
noise, and the shared library pulls the estimate the other way again.

### The source paper's own agreement statistics, which are better evidence

An earlier version of this report said a trustworthy ceiling "needs replicates
nobody has published". **That was wrong twice over**, and the correction came
from Jacky reading the source paper more carefully than we had.

Guo et al. 2018 ran **two biological replicates per condition, by independent
transformations**, and published the *agreement between them* even though the
per-guide values were averaged away. Verified against the paper (*NAR* 46:7052,
Figure 2):

| comparison | n | statistic | what it measures |
|---|---:|---:|---|
| replicate vs replicate (Fig 2b) | 2 libraries | **R² > 0.78** | same library, same enzyme, independent transformations |
| genome-wide vs an independent **tiling** library (Fig 2c) | 901 shared guides | **R² = 0.771** | a *different* library, a separate experiment |
| screen vs individual colony counting (Fig 2d) | 15 sgRNAs | R² = 0.840 | an orthogonal assay |

**And the label we train on is a two-replicate average.** The paper states the
read counts "were averaged as the geometric mean". That matters, because the
attenuation argument above bounds prediction of *one noisy observation*; our
target is already the mean of two. Under Spearman–Brown the reliability of a
2-item mean is `2r / (1 + r)`, which is higher than `r`. The two-enzyme estimate
therefore erred in both directions at once: it treated enzyme biology as noise,
and it ignored that the target is an average.

Two revised estimates, and they bracket differently from the old pair:

| route | r | reliability | ceiling |
|---|---:|---:|---:|
| **tiling library (Fig 2c)** — correlation of two independent experiments' published scores, so this *is* the reliability of the quantity we predict, no correction needed | 0.878 | 0.878 | **0.937** |
| replicate agreement (Fig 2b) + Spearman–Brown for the 2-replicate mean | >0.883 | 0.938 | **0.969** |

**The usable range is therefore ≈0.94–0.97, not 0.90–0.93**, and the tiling
figure is the one to lead with: it is a test–retest of the exact published
quantity with a **different library**, which is precisely the objection the
two-enzyme estimate could not answer.

**Why the tiling figure is the better number, not just the smaller one.** The
two routes are not equally safe, and the direction of each error is known.

The replicate route needs the activity score to inherit the **read counts'**
reliability. It cannot. Guo's label is the `|Z|` of a **Cas9/dCas9 ratio** — two
count measurements, not one — and a difference of two measurements is never more
reliable than its components, usually much less, by exactly how much true signal
the two arms share. The dCas9 arm exists *to* share the library-abundance
structure, so that fraction is high. Taking each log-count's reliability at the
Fig 2b floor of 0.883:

| true signal shared between the Cas9 and dCas9 arms | reliability of the ratio score | ceiling |
|---:|---:|---:|
| 0% (what the replicate route assumes) | 0.883 | **0.968** |
| 30% | 0.841 | 0.956 |
| 50% | 0.791 | **0.940** |
| 70% | 0.694 | 0.905 |

So **0.968 is what you get only if the control arm shares nothing**, which is the
one thing it is designed not to do. It is an upper bound, and plausibly a loose
one.

The tiling route needs none of that. Both sides of Fig 2c are **already the
published score**, so r = 0.878 *is* its reliability — no propagation of count
noise, no Spearman–Brown step, no assumption about the ratio. Its own sampling
error is small: on 901 guides, r ∈ [0.863, 0.893], so the ceiling is
**0.937 [0.929, 0.945]**. Its bias is in the opposite, conservative direction —
a different library adds design differences that are not label noise, so the true
reliability is a little higher than 0.878.

**So: quote 0.94 as the working figure, with 0.94–0.97 as the range**, and say
which is which. Two residual caveats apply to both routes: the published values
are Pearson while the ceiling is quoted in Spearman, so the conversion
approximates; and Fig 2b is an inequality across ten libraries, so 0.883 is a
floor rather than an estimate.

**SLICER at 0.707–0.721 is therefore about 77% of the way to the working
figure of 0.94**, or 74% against the loose upper bound, against the ~80% the old
bracket implied. The headroom is larger than this report
previously claimed, which makes the saturation results in the rest of this part
more interesting rather than less: the model stops improving well short of a
ceiling that is further away than we thought.

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

### How we know the model is not simply failing to use what it has

Everything above says the *features* are exhausted. That only follows if the
model is actually using them, and "the model is underfitting" and "the features
are out of information" predict the same flat ablation table. The test that
separates them is a **residual test**, and it is the single most useful
diagnostic in this project.

Fit the champion, take its **out-of-fold** predictions — so every residual is for
a guide the model has never seen — and correlate every column we hold against
those residuals. A column that still tracks the errors is information the model
was given and failed to extract. A table of near-zeros means it is out of
material.

| what the model was given | strongest \|ρ\| with its residuals, among columns it **had** | among columns **held back** |
|---|---:|---:|
| all 6,232 published columns | **0.028** (none above 0.10) | **0.147** — the flank family |
| the same plus the flank family | **0.031** (none above 0.10) | — nothing left to hold back |

**The left column is the answer: not one of 6,480 available columns retains a
residual correlation above 0.031.** The model is not underfitting; there is
nothing left in those columns to extract. (Two independent runs of this
test are on file, `headroom.csv` and `representation_headroom.csv`; the maximum
is 0.028–0.035 depending on which arm is held back, and elsewhere this report
quotes the looser bound, 0.035.)

**And the right column is why the test can be trusted** — it has a demonstrated
positive control. When the flank family was held back, the test flagged it at
0.147, four times anything the model had, and adding it was then worth +0.080.
A diagnostic that only ever returns "nothing here" is not evidence; this one
found the one thing that was there.

### Is that still true of the reduced 427-column set?

It is a fair worry: the reduced representation throws away 5,805 columns, and
although they are all functions of the 20-mer, the 427 carry only its first and
second order. A trimer or tetramer indicator could in principle track something
mono- and dinucleotide indicators cannot express. Re-running the residual test on
the subset's own out-of-fold residuals (`representation_headroom.csv`):

| model | out-of-fold ρ | strongest \|ρ\| among columns it has | among the columns dropped |
|---|---:|---:|---:|
| all 6,232 | 0.531 | 0.028 | 0.147 *(the flank family)* |
| the 427 | 0.525 | 0.022 | 0.146 *(the flank family)* |
| 427 + the flank family | 0.597 | 0.028 | **0.070**, none above 0.10 |

Read the first two rows together: **dropping 5,805 columns did not change what
the model is missing.** In both cases the strongest residual correlate is the
same flank column at the same magnitude — the reduction did not create a blind
spot, it removed restatement.

The third row is the honest limit. With the flank family in, the dropped columns
top out at **0.070** — below the 0.10 threshold this project uses for "worth
chasing", but still **2.5× the strongest correlate among the columns the model
has.** So the reduced set does leave a little unused, consistent with its −0.01 ρ
cost, and the leftovers are identifiable: they are **higher-order k-mer
indicators**, exactly the information a mono-plus-dinucleotide encoding cannot
express.

That squares with the depth result above. There is a sliver of higher-order
sequence signal, and handing the model 1,146 trimer indicators to reach it costs
more in selection dilution than the sliver is worth. **Both facts are true: the
information exists, and taking it is not profitable.**

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

**The GC effect is not absent in human cells — it is a different shape, and we
measured it.** Fitting GC and GC² against the cut score on both matrices:

| | linear ρ | linear R² | quadratic R² | curved term | turning point | sits at |
|---|---:|---:|---:|---:|---:|---|
| human | +0.017 | 0.0002 | **0.0042** | p < 10⁻¹⁵ | GC **0.571** | the **60th** percentile |
| *E. coli* | −0.201 | 0.0407 | 0.0452 | p < 10⁻¹⁵ | GC 0.294 | the **1st** percentile |

Both have a significant curved term; what differs is **where the turning point
falls relative to the data.** In human guides it is interior — a genuine optimum
near GC 0.57, with decile means rising 0.235 → 0.262 and then falling to 0.192 at
GC 0.78. In *E. coli* it sits below the 1st percentile, so across the bulk of the
data the relationship is a monotonic decline. **The quadratic fit explains 22×
more human variance than the linear one**, which is why a linear ρ of +0.017 read
as "no relationship at all". So the honest statement is not that the mechanism
fails to transfer but that **it is replaced by one of a different shape.**

**And histones are the likely reason, though that part we have not tested.**
**Nucleosomes physically block Cas9** and bacteria have none, so in human cells a
major determinant is whether the target is *accessible* — a variable absent from
our data and not inferable from sequence. Targets in promoter regions, which are
kept open, cut better than intergenic ones.

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
4. **Interpretability measured, diagnosed, and then fixed** — importance
   rankings on this matrix disagree between models at ρ 0.05 because the matrix
   is redundant and because split gain additionally rewards a column for being
   binary. Measured with SHAP across seven model families the disagreement
   largely goes away (ρ 0.50 per column, 0.80 per block), which both identifies
   the artefact and supplies the method that avoids it. The consequence for the
   inherited paper stands: its biological claim was read off a gain ranking.
5. **A faithfulness probe shown to be uninterpretable on a redundant matrix**
   — the standard "drop the top features and refit" test returns ≈0 on the full
   matrix for every model, and −0.015 against −0.001 at random once the
   restatement is removed, at a matched share of attributed importance. It is a
   measurement artefact that this project itself reported as a finding for
   several weeks.
6. **The published feature set reduced to 7% of its columns** — 427 of 6,232 are
   statistically indistinguishable from all of them, and the matrix decomposes
   into sequence, non-sequence and restatement. This is the redundancy argument
   turned from an interpretation into a measurement.
7. **The flank effect characterised** — a gradient rather than a motif, peaking at
   250–500 letters, downstream-weighted, transferring across enzyme and organism
   at ~90% and across kingdoms at 0%.
8. **The ceiling argument tested rather than assumed**, twice. The two-enzyme
   derivation rests on a false independence assumption — the disagreement between
   WT-SpCas9 and eSpCas9 is predictable at ρ 0.546, so it is enzyme biology, not
   noise. And the source paper's own replicate and tiling-library agreement puts
   the limit at ≈0.94–0.97 rather than 0.90–0.93, once the published score is
   recognised as a two-replicate mean.
9. **Measured boundaries**: where the method stops, what the data cannot answer,
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
- **A ceiling precise enough to quote as a number.** The ≈0.94–0.97 range rests
  on Pearson R² values converted to a Spearman bound, one of them published only
  as an inequality. Narrowing it needs per-guide replicate values, which exist in
  SRA (PRJNA450978) but only as raw reads.

---

## Part 12 — Next steps and how to frame the write-up

### Next, in order of value for effort

1. **At data freeze, run crisprHAL's search to its own plateau on a free Colab
   GPU.** Both models are now tuned and the margin survives (+0.016, 5/5 folds),
   but theirs had 10 draws against our 627. Running *their* search until *their*
   patience fires is ~55 fits — **2–4 GPU-hours, one free Colab session** — not
   the 627-draw match that would cost real money. Do it once, after the rows,
   label and folds stop moving; before that it measures a moving target. It
   replaces an extrapolation with a measurement and changes no conclusion, so it
   is a reviewer-defence item, not a prerequisite.
2. **Put a test-set interval on the cross-organism transfers.** There are now
   three model seeds per cell, but all three share one test set, so the spread
   understates the uncertainty. A bootstrap over the test set is an hour and it
   is the weakest-supported table in the report.
3. **Find or generate a genome-wide screen in a distant bacterium.** The one open
   question analysis cannot close, and it needs wet-lab work.

### Closed since the last version

- **The upstream/downstream asymmetry is not a transcription effect.** Half of
  that test turned out to be impossible — the library is 99.97% template-strand,
  so there is no orientation to condition on — and the half that was possible
  found the asymmetry flat across expression bands (Part 5). An open question with
  one fewer candidate answer.
- **The human GC relationship is U-shaped**, with an interior optimum at GC ≈ 0.57
  and a quadratic fit explaining 22× more variance than the linear one (Part 10).
  "The mechanism does not transfer" is now "the mechanism is replaced by one of a
  different shape".
- **SHAP is adopted as the reported importance method**, and the figures that
  showed gain have been regenerated (Part 9, `results/attribution_summary_shap.csv`,
  `interpretability.csv`). It also removed the reason to switch the champion to
  LightGBM.
- **The published feature set is reduced to 427 columns** with no measurable loss,
  and the residual test confirms the reduction creates no blind spot (Part 4,
  Part 10).

### Could a GPU settle the comparison outright? Yes — and it should not be a budget item

The one thing limiting the head-to-head is that crisprHAL got 10 draws and
SLICER got 627, because their fits cost ~44 min of CPU each against our 2 s.

**The honest size of the job is much smaller than "match 627 draws".** The trace
analysis above is the reason: 40 draws captured **100%** of SLICER's eventual
gain, and the remaining 587 only confirmed the plateau. So an equal-footing
search for crisprHAL means running *its* search until *its* patience criterion
fires — on the order of **40–55 fits**, not 627.

| | fits | CPU-equivalent | on a rented T4/L4 |
|---|---:|---:|---:|
| matching SLICER draw-for-draw | 627 | ~460 h | 30–50 h |
| **running their search to its own plateau** | **~55** | **~40 h** | **2–4 h** |

Two to four GPU-hours fits inside a **single free Colab session**, and this
project already has a Colab notebook. So the correct recommendation is **not** to
fund GPU time; it is to run it free, once, at the right moment.

**A TPU would be the wrong instrument.** TPUs earn their keep on large dense
matrix multiplication under XLA. This model's recurrent branch is a bidirectional
GRU — a sequential loop, which is the shape that maps worst onto them. Fiddlier
to set up and probably slower per fit than a mid-range GPU.

#### When to run it, and when re-tuning is actually required

Re-tuning is a recurring cost if it is done while the inputs are still moving, so
the rule matters more than the budget:

| change | invalidates |
|---|---|
| rows, label, or fold scheme | **both** models' tuned configurations |
| the feature set | **SLICER's only** — crisprHAL never sees the matrix |
| documentation, new ablations, analyses outside the curated arm | **neither** |

So: **tune once, at data freeze.** Before that it measures a moving target. The
chosen configurations are recorded (`crisprhal_tuned_folds.csv`, `tuning.csv`),
so reproducing a number when the data has not changed is a **re-fit, not a
re-search**.

#### What it would and would not buy

It would replace the extrapolation — "their converged gain is probably ≈+0.012,
so the lead is probably ≈0.012 rather than 0.016" — with a measurement on equal
footing. It would **not** change the conclusion, and it would **not** change the
cost comparison, which is already labelled CPU-to-CPU.

**And the headline does not depend on it at all**, which is the strongest reason
not to spend money or urgency on it. The primary comparison should be
**untuned against untuned** — +0.0111, 5/5 folds, p = 0.0002 — because default
settings are what anyone installing either tool actually gets, and that
comparison carries no search-size caveat whatever. The tuned comparison is then
a *robustness check*: it shows the margin is not an artefact of whose defaults
suited this screen, and 10 draws against 627 is adequate for that purpose.
Sharpening +0.016 to +0.012 changes no claim in this report, and the framing
advice below says the margin should not be the headline in the first place.

### A note on tooling, deliberately parked

A guide-design tool or web UI is **not** a deliverable of this project. The
findings are the contribution, and the field already has several predictors, so a
tool without the findings would add little. Everything needed to build one later
exists — `featurise.matrix()` turns any 20-mer into the published representation,
and the model runs in under a second — but it is a separate piece of work and
should not compete with finishing the paper.

### How to frame it

**Not** "a better guide predictor". A 0.016 lead on one organism's screen invites
the one comparison this project loses — against a full-time lab with a GPU — and
stakes the paper on its least interesting number.

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
- and, as a by-product, a small measured lead over the state of the art
  (+0.016, 5/5 folds, both tuned) at about 1/74th of the CPU time — with the
  search asymmetry quantified rather than waved away.

### Honest caveats to carry into the paper

- **A small lead, not a rout.** +0.016 with both models tuned, 5/5 folds — but
  their search was 60× smaller than ours, which the trace analysis suggests
  accounts for about a quarter of it.
- **The speed advantage is CPU-to-CPU.** crisprHAL 2 as published is GPU-trained.
- **Cross-species is untested**, and the competition's position is no better,
  which does not make ours good.
- **Cross-enzyme evidence is narrow** — three point mutations on an identical
  library.
- **The ceiling is a bracket, not a number** — ≈0.94–0.97 from the source
  paper's own agreement statistics, revised upward from 0.90–0.93 once the
  two-enzyme assumption was dropped and the label's two-replicate averaging
  accounted for.
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
