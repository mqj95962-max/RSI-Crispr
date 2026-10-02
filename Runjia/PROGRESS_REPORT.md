> **Synced copy — do not edit here.**
> The source of truth is `docs/PROGRESS_REPORT.md` in Runjia's RSI09 working
> folder, which also holds the companion files this report cites (`RESULTS.md`
> for every number with its settings, `FINDINGS_LOG.md` for the chronological
> record). Those are not in this repository, so references to them below are
> plain names rather than links.
> Synced 2026-10-02.

# What determines whether a CRISPR guide works in *E. coli*

RSI09, *Enhancing sgRNA efficiency in CRISPR-Cas9 genome editing*.

Work covered: 8–15 September and 1–2 October 2026. Written to be read by
someone outside the project — a supervisor, a reviewer, a judge — with no
assumed background in molecular biology or machine learning. Restructured
2 October around *what is new* rather than the order things were tried, and
merged with the former `NOVELTY.md`.

Companion documents: `RESULTS.md` (every
number, with its settings), `FINDINGS_LOG.md` (the
chronological record, including superseded claims).

**One-paragraph version.** CRISPR-Cas9 cuts DNA at a site chosen by a short
"guide" RNA, and some guides work ten times better than others. This project
started from a 2023 paper that described each guide with 6,232 numbers,
including quantum-chemical properties of the DNA bases, and asked what *else*
determines whether a guide works. Nine further families of features were built
and measured. The useful answer turned out not to be a better predictor but a
sharper account of the problem: most well-motivated ideas **cannot** help, and
that is now provable in seconds rather than discoverable in a week; the one
thing that does help is the DNA *surrounding* the target, and it works as a
smooth compositional gradient rather than a pattern; two of our own positive
results were overturned by their own controls; and the interpretability claimed
across this literature does not survive being measured. As a by-product the
project reaches parity with the best published bacterial model on identical
data, with a model that trains in 0.3 seconds on a laptop.

---

## Contents

- [Part 0 — The words you'll need](#part-0--the-words-youll-need)
- [Part 1 — Three reference points: the source paper, the state of the art, and us](#part-1--three-reference-points-the-source-paper-the-state-of-the-art-and-us)
- [Part 2 — What made any of this possible: decoding the published matrix](#part-2--what-made-any-of-this-possible-decoding-the-published-matrix)
- [Part 3 — The one large effect, and the shape of it](#part-3--the-one-large-effect-and-the-shape-of-it)
- [Part 4 — The rule: when a new feature cannot possibly help](#part-4--the-rule-when-a-new-feature-cannot-possibly-help)
- [Part 5 — Three of our own positive results, overturned by their own controls](#part-5--three-of-our-own-positive-results-overturned-by-their-own-controls)
- [Part 6 — Interpretability, measured rather than asserted](#part-6--interpretability-measured-rather-than-asserted)
- [Part 7 — How much room is left: four independent ceilings](#part-7--how-much-room-is-left-four-independent-ceilings)
- [Part 8 — What is new here, and what is imported](#part-8--what-is-new-here-and-what-is-imported)
- [Part 9 — Where we corrected ourselves](#part-9--where-we-corrected-ourselves)
- [Part 10 — What to do next, and how to frame the write-up](#part-10--what-to-do-next-and-how-to-frame-the-write-up)

---

## Part 0 — The words you'll need

Skip this if you know them. Everything later uses these without re-explaining.

### The biology

**CRISPR-Cas9** — molecular scissors. **Cas9** is a protein that cuts DNA. It
does not know where to cut on its own.

**Guide RNA (sgRNA)** — the address label telling Cas9 where to cut. It carries
a 20-letter sequence, and Cas9 cuts where the DNA matches it. We mostly say
**guide**.

**Protospacer** — the 20 letters of DNA being targeted.

**PAM** — a mandatory three-letter suffix in the DNA immediately after the
target. It must read **NGG** (any letter, then G, then G) or Cas9 will not cut.

**DNA letters** — DNA is a string over A, C, G, T. **GC content** is the
fraction that are G or C. GC-rich DNA is physically harder to pull apart, which
turns out to matter.

**Flank / flanking DNA** — the DNA *surrounding* the target, on either side. The
guide does not match it, so the obvious assumption is that it is irrelevant. The
central finding of this project is that it is not.

***E. coli*** — a common gut bacterium, and the standard workhorse for this
kind of experiment.

**A screen** — an experiment testing tens of thousands of guides at once.
Cutting a bacterium's chromosome usually kills it, so you grow a mixed
population, sequence what survives, and a guide whose cells *disappeared* cut
well. The resulting efficiency number is the **cut score**.

Worth remembering: the screen measures *survival*, not cutting. Survival mixes
cutting with DNA repair, growth rate and sequencing noise, so the cut score is
a noisy, indirect measure of what we care about. Part 7 quantifies exactly how
noisy.

### The computation

**Feature** — one number describing a guide (its GC content, the temperature at
which its DNA comes apart). A model never sees the guide, only its features.
Choosing features is most of the work.

**Model** — a program that learns a pattern from examples: shown thousands of
guides and their cut scores, predict the score of a guide it has never seen.

**Spearman correlation (ρ, "rho")** — whether the model gets the **ordering**
right. 0 is random, 1 is perfect. **This is the number to watch**, because
nobody needs a guide's exact efficiency; they have ten candidates and want the
best one.

**R²** — how much of the variation a model explains. Reported here too, but it
answers a question nobody asks. Watch for places below where R² improves and ρ
does not — that gap is itself a finding.

**Pick percentile** — the most honest reading. Take ten candidate guides, let
the model choose, and ask where that guide really falls. Picking blindly gives
the 50th percentile; our current model gives the **76th**; a perfect model would
give the 91st.

**Cross-validation** — hide part of the data, train on the rest, test on the
hidden part, repeat. The discipline that stops you fooling yourself.

**Ablation** — add one group of features, re-measure, see whether the model
actually got better. How every idea here was tested rather than assumed.

**Quantum chemical features** — the source paper's idea and this project's
starting point. Physics calculations describing the electrons in each DNA
letter: how tightly neighbouring letters stack, how much energy their bonds
hold.

**Gradient boosting (LightGBM, XGBoost)** — models built from many small
decision trees. Fast, run on an ordinary laptop, features have names you can
inspect.

**Neural network (CNN, BiGRU)** — a more flexible family that learns from raw
data. Usually needs a GPU and is much harder to interpret. The best published
competitor is one of these.

---

## Part 1 — Three reference points: the source paper, the state of the art, and us

Everything in this report is positioned against three fixed points, so it is
worth putting them side by side first.

| model | year | organism | training guides | ρ |
|---|---|---|---:|---:|
| Guo et al. | 2018 | *E. coli* | ~61,000 | 0.542 |
| **Noshay et al.** — the paper this project extends | 2023 | *E. coli* | 40,468 | Pearson 0.502, R² 0.249 |
| crisprHAL 1 | 2023 | *E. coli* | 40,308 | 0.627 |
| **crisprHAL 2** — the state of the art | 2026 | *E. coli* | 33,495 curated | **0.697** |
| DeepCC9 | 2026 | **human** | 55,604 | 0.861 |

Two warnings about that table. **The human number is not comparable** — human
screens count edits at a single site with a clean readout, while bacterial
screens measure survival. crisprHAL showed human models applied to bacteria
score ρ −0.2 to 0.1, worse than useless. And **the paper this project builds on
is not the leader**: Noshay et al.'s Pearson 0.502 sits below the 2018
baseline.

So the three reference points are:

- **The source paper** — 6,232 features per guide, mostly quantum chemistry,
  describing the 20-letter target. Its framing is explainable AI.
- **crisprHAL 2 (the bar)** — a convolutional + recurrent neural network,
  GPU-trained, reading 378 raw nucleotides around each guide: 189 upstream, the
  20-letter target, 3 for the PAM, 166 downstream.
- **This project** — the source paper's features, plus nine families of our
  own, in a gradient-boosted tree.

### Where this project stands

| | ρ | tested on |
|---|---:|---|
| the source paper, as published | 0.502 (Pearson) | 40,468 guides |
| our reproduction of it | 0.527 | 13,880 |
| + flanking-DNA features | 0.609 | 13,880 |
| all nine feature families | 0.641 | 13,880 |
| **flanking DNA, on the curated dataset** | **0.707 ± 0.008** | **33,567** |
| *crisprHAL 2, as published* | *0.697* | *33,495* |
| *crisprHAL 2, re-run by us on our own folds* | *0.697 ± 0.007* | *33,567* |
| the ceiling any model could reach (Part 7) | ~0.90 | — |

**The comparison is now against a model we ran ourselves, not a published
number.** This was the largest outstanding caveat in the September write-up. We
imported crisprHAL 2's architecture, hyper-parameters and 48 training epochs
unchanged from their repository and trained it on our rows and our folds.

The harness checks out twice: on *their* shipped split it reproduces their
published hold-out of 0.695 at **0.694**, and its mean over our folds is
**0.6971** against their published 0.697. Then, paired on identical rows,
identical folds and the same label:

| | ρ |
|---|---:|
| ours (LightGBM, base + flanking DNA) | **0.7078 ± 0.0078** |
| crisprHAL 2, re-run here | 0.6971 ± 0.0067 |
| paired difference | **+0.0107, ahead in 5 of 5 folds** |

**State this as a small but consistent advantage, not a victory.** +0.011 ρ is
a fraction of one pick percentile and would not change anyone's guide choice.
Their arm is one seed against our two, their hyper-parameters were tuned by
them on their own split rather than re-tuned here, and crisprHAL 2 generalises
across organisms and nucleases while this project is *E. coli* SpCas9 only.

What makes it worth reporting is the cost: LightGBM fits in **0.3 seconds on a
CPU**, against **125 minutes** for crisprHAL 2's five folds on the same
machine — roughly a 25,000-fold difference — and every feature it uses has a
name.

### What the three actually differ in

The useful contrast is not accuracy, it is what each one is doing.

**Information: we and crisprHAL 2 largely see the same DNA.** They feed their
network 189 nt upstream and 166 downstream; almost everything our flank
features measure lies inside that window. We are not showing the model DNA they
never see. (Two small exceptions: our windows reach 1 kb, and we apply the
source paper's quantum-chemical tables to the flanking bases, which nobody has
done.)

**Representation: completely different, and this is the scientific point.**
They give raw bases to a pattern-detector. We give a tree the *composition* of
windows — GC fraction, purine fraction, longest homopolymer run, averaged
quantum tensors at 50, 250, 500 and 1,000 nt. Part 3 shows that difference is
worth a factor of four, and explains why.

**The question asked: entirely different.** crisprHAL 2 is a prediction paper —
build the best bacterial predictor and ship it as a tool. This is a
**mechanism-and-method** paper: what determines cutting efficiency, and how do
you tell in advance whether a proposed determinant can possibly help? Most of
Part 8 has no counterpart in crisprHAL 2 because it answers questions that
paper never asks.

---

## Part 2 — What made any of this possible: decoding the published matrix

This is an enabling contribution rather than a biological one, but nothing else
in the report could have happened without it.

### The problem

The source paper described 13,880 guides with 6,232 numbers each, and two
things made that dataset nearly unusable by anyone else:

- **It never recorded the DNA sequence of any guide** — only derived numbers.
- **5,887 of the 6,232 columns were named `V1`, `V2`, `V3`…**, with no
  description of what they measured.

So the dataset applied only to the guides its authors had already processed.
You could not use their method on a new guide.

### What we found

**The DNA sequence is hidden in the numbers.** One column records the number of
electrons per DNA letter. Across all 13,880 guides it takes exactly four values
— 42, 48, 50, 56 — the valence-electron counts of C, T, A and G. That column
*is* the sequence in a different alphabet; read position by position, it
reconstructs the guide.

**It checks out against the real genome.** 13,879 of 13,880 guides decoded (one
row has a missing value), every one was located in the *E. coli* genome with
the required NGG PAM in the right place, and 13,877 appear exactly once. As an
independent check, locating where bacterial DNA copying begins from GC skew
alone put it at position 3,923,620 against a textbook value near 3,923,800.

**The quantum tables are recoverable exactly.** Each quantum column is a fixed
lookup: a given short run of letters always yields the same number. All 4 single
letters, 16 pairs, 64 triples and 256 quadruples appear, **with no
contradictions anywhere**. We have not approximated their method; we have
recovered it.

**The anonymous columns were identified.** 5,853 of the 5,887 `V####` columns
are "is there letter X at position Y" indicators — none ambiguous, **verified
against 1.17 million rebuilt values at 100% agreement**.

### Why it matters

**99.0% of the published matrix (6,169 of 6,232 columns) now regenerates from a
guide's 20 letters alone.** A resource that worked only for its authors works
for anyone, on any bacterial Cas9 screen. The head-to-head in Part 1 was
impossible without it — it is what let us compute the source paper's
representation for crisprHAL 2's guides. And, as Part 4 shows, the *way* it
decodes is also what makes most feature ideas provably useless.

---

## Part 3 — The one large effect, and the shape of it

### The effect

The DNA *surrounding* the target should not matter: the guide does not match it
and Cas9 does not read it. It matters more than anything else we measured —
**+0.080 ρ** on the original data and **+0.164** on the curated version of the
same data.

Three things are settled about it:

- **It reaches a few hundred letters.** The effect peaks at 250–500 letters out
  and has faded by 1,000. No point looking further.
- **The flank wants the opposite of the target.** GC-rich *targets* cut worse
  (−0.20); GC-rich *flanks* cut slightly **better** (+0.05). Opposite
  directions at different distance scales, so these are two physical effects,
  not one.
- **Downstream matters 3.4× more than upstream**, and on clean data the distant
  signal matters *more* than the nearby one.

### The shape: a gradient, not a pattern

This is the project's main biological claim, and it is a claim about *form*
rather than magnitude.

We built a convolutional + recurrent network of the same kind as crisprHAL 2
and gave it raw DNA:

| what the model reads | ρ |
|---|---:|
| neural net, the 20-letter target only | 0.498 |
| gradient boosting, the published features | 0.527 |
| neural net, target + 100 letters either side | 0.516 |
| gradient boosting, published features + our flank features | **0.609** |

Raw flanking DNA gains the network **+0.019**. Windowed composition of the same
region gains the simpler model **+0.164** on the same rows and label — four
times as much from the same DNA.

**Why: the flank effect is a smooth average, not a motif.** There is no short
recurring pattern to find. A convolution with a 5-letter kernel is built to
detect local patterns and is the wrong instrument for "how GC-rich are the next
500 letters"; a windowed mean computes exactly that. crisprHAL 2's architecture
cannot express such an average directly — it must approximate it — which is a
plausible mechanistic reason a GPU network on 33,495 rows ends up level with a
0.3-second tree.

**That claim is ours.** crisprHAL 2 establishes *that* long context helps; we
establish *what shape* the signal has, and therefore how to encode it. Their
architecture cannot easily test it.

### It is not a property of this enzyme — and the species question cannot be answered with the data that exists

Everything above was measured on one screen: *E. coli*, WT-SpCas9. So the
obvious objection is that the effect belongs to that enzyme or that genome
rather than to DNA. crisprHAL ships three further screens. Working out *which
of them can answer* turned out to be most of the work.

**The clean result: a different enzyme, same organism.** eSpCas9 is an
engineered high-fidelity version of Cas9, screened genome-wide in *E. coli*
(all 59,489 guides located across all five chromosome arcs). Tested under the
strictest scheme this project has — hold out a whole fifth of the chromosome
and make the model extrapolate:

| what the model sees | ρ | gain |
|---|---:|---:|
| published features only | 0.685 | — |
| + the immediate ±10 letters | 0.728 | +0.043 |
| + composition of 50–1000 letter windows | 0.761 | **+0.076** |
| + all flank features | **0.789** | **+0.104** |

Every fold improves (p = 7.4e-06), and the long-range windows again beat the
immediate context (+0.033, 5/5 folds) — the same ordering found in *E. coli*
WT. Loose and strict cross-validation agree to three decimals.

**So the smooth-gradient finding is not an artifact of SpCas9.** Within
*E. coli*, it survives changing the enzyme — in size *and* in shape.

**The species question: the available data cannot answer it.** The one
other-organism screen is *C. rodentium*, and before trusting it we measured
where its guides sit. They are not spread across the chromosome:

| screen | genomic spread | density |
|---|---|---:|
| eSpCas9, *E. coli* | all five arcs, genome-wide | — |
| TevSpCas9, *C. rodentium* | **229 kb — 4.3% of the chromosome** | **110 guides per kb** |

It densely tiles a single locus, and that is the authors' own design — their
paper says a 236 kb fragment was screened. Our measurement matches.

The consequence is specific. Across 236 kb, flank composition barely varies, so
there is almost nothing for a *correlation* to work with: flank GC correlates
**−0.022** with cutting there, against +0.159 and +0.138 in the two *E. coli*
screens. Family A still helps (+0.060, surviving 10 kb blocking), but its two
halves become equal rather than long-range-dominant, and no individual
long-range feature clears the project's 0.035 noise level.

Read on its own, that says the gradient is absent in *C. rodentium*. The
transfer experiment below shows it is not absent — it is **invisible to a
correlation measured inside a 236 kb window**, which is a different thing.

**What did replicate across both organisms** is the mechanism from the next
section: GC-rich *targets* cut worse everywhere, and most strongly of all in
*C. rodentium* (−0.379 against −0.147 here). That is a claim about the physics
of prying a duplex open, so it should be indifferent to the organism, and it is.

So the standing limitation is **narrowed rather than retired**: the flank
effect is a property of DNA rather than of SpCas9, demonstrated genome-wide
under the strictest CV; whether its long-range half crosses species is still
open, and answering it needs a genome-wide screen in another bacterium.

A fourth screen, TevSaCas9, was deliberately **not** used: its enzyme reads a
different PAM, so this pipeline's position labels would have been misaligned
against it. It would have produced numbers.

### A model trained on one organism does rank another's guides

The tests above all train and test within a single screen, which answers "are
these features informative here?" — not "does a model built here work there?".
That second question is the one a reader assumes, so it was measured directly:
fit on every guide of one screen, predict every guide of another.

| trained on | tested on | published features only | + flank features | that screen's own model |
|---|---|---:|---:|---:|
| *E. coli* (WT) | ***C. rodentium*** | 0.626 | **0.700** | 0.764 |
| ***C. rodentium*** | *E. coli* (WT) | 0.478 | **0.628** | 0.707 |
| *E. coli* (eSp) | ***C. rodentium*** | 0.587 | 0.656 | 0.764 |

**Cross-organism transfer works, at about 90% of a locally-trained model.** The
guide sequences are completely different and the organisms are different, so
this is real generalisation.

**And the flank features are what carry it.** They improve every direction, and
most in the hardest one: going from *C. rodentium* to *E. coli*, the published
features alone manage 0.478, and adding flank features recovers it to 0.628.
The flank encoding is the part of the model that survives changing organism.

A fourth pairing — *E. coli* WT to *E. coli* eSpCas9 — reaches 0.689, but those
two screens use **the same guide library**, so the model has seen every test
sequence before. It isolates the enzyme change cleanly and says nothing about
unfamiliar DNA; it is not a generalisation result.

**Which half of the flank features makes the crossing?** The long-range
windows, in both directions — and overwhelmingly in the harder one:

| trained on | tested on | base | + local ±10 nt | + long-range windows | + both |
|---|---|---:|---:|---:|---:|
| *E. coli* | *C. rodentium* | 0.626 | 0.661 | **0.665** | 0.700 |
| *C. rodentium* | *E. coli* | 0.478 | 0.537 | **0.586** | 0.628 |

In the second row, 88 long-range columns are worth +0.107 against 260 local
columns' +0.058. That model was trained on guides from a single 236 kb window —
the very window where a correlation finds nothing — and its long-range features
still rank guides across an entire other genome.

**So a relationship can be learnable from data too narrow to reveal it.**
Measuring a correlation needs spread in the data; fitting a function needs
rather less. Discovery and verification have different requirements, and
treating them as the same is what made the earlier reading look conclusive.

**Does the *E. coli* model then beat *C. rodentium*'s own model on
*C. rodentium*?** No — but the part that matters does. Both are out-of-sample
on the same guides:

| predicting *C. rodentium* | its own model | the *E. coli* model |
|---|---:|---:|
| published features only | **0.704** | 0.626 |
| + flank features | **0.764** | 0.700 |
| *what the flank features added* | *+0.060* | ***+0.074*** |

The local model wins overall by 0.064, and **all of that advantage is in the
base features** — it knows its own screen's label scale, nuclease and
sequence-to-score mapping. But the *E. coli*-trained model gets **more out of
*C. rodentium*'s flanks than *C. rodentium*'s own model does.**

That suggested the flank relationship might be **better** learned from a
genome-wide screen in the wrong organism than from a narrow one in the right
organism — and so that importing it should improve the local model. We tested
that, and it is wrong on both counts.

The comparison itself was flawed: the two gains sit on different starting
points (0.626 against 0.704), and the model with more room to improve shows a
bigger gain for free. So we ran the clean version — train a flank-only model on
*E. coli*, use its prediction as a single extra column for a *C. rodentium*
model that keeps its own labels:

| what the *C. rodentium* model gets | ρ |
|---|---:|
| published features only | 0.704 |
| **+ one imported column from the *E. coli* flank model** | **0.751** |
| + its own 348 flank features | 0.764 |
| + both | **0.764** |

**Adding the import on top of local flank features is worth +0.0002 —
nothing.** Whatever the *E. coli* model knows about flanks, *C. rodentium*'s
own 236 kb of flanking DNA already supplies. The hybrid idea fails.

**But a single imported number does 78% of the work of 348 local ones.** That
is a compression result rather than new information — and it makes the shared-
relationship claim stronger, not weaker. If the flank-to-efficiency function
were organism-specific, a model fitted on *E. coli* could not stand in for 348
locally-computed columns at 78% strength. It can. The relationship really is
common to both organisms; the narrow screen simply is not short of it.

### What the competition's cross-species claim actually rests on

Since "they generalise across species and we do not" has been this project's
standing concession, it is worth checking. Two things turn out to be wrong with
how it was stated here.

**It is the earlier paper's claim.** The cross-species results belong to
crisprHAL 1 (*Nat Commun*, 2023). crisprHAL 2 — the model Part 1 benchmarks
against — is a data-curation paper, described by its own repository as
rebuilding two prior *E. coli* datasets.

**And "three species" is doing a lot of work:**

| organism | what the data actually is | guides |
|---|---|---:|
| *E. coli* | training and held-out test | 45,010 / 7,821 |
| *C. rodentium* | a genuine screen in that organism | 31,796 |
| *S. enterica* | **a 2 kb piece of one of its genes, cloned into *E. coli*** | **~300** |

The *Salmonella* test is that organism's **DNA placed inside an *E. coli*
cell** — about 300 guides across 2 kb. That tests whether the model copes with
unfamiliar sequence, which is worth knowing, but it is not a screen in another
organism: there is no *Salmonella* chromosome and no *Salmonella* cell. And the
*C. rodentium* data is the same set measured above, covering 4.3% of a
chromosome.

So the fair statement is not that they have solved cross-species prediction and
we have not. It is that **nobody has the dataset the question needs** — a
genome-wide Cas9 screen in a bacterium other than *E. coli* — and this project
is the one that measured why the existing substitutes cannot stand in for it.

*(Verified against the* Nat Commun *full text: a 236 kb* C. rodentium
*fragment, and 296 sgRNAs against a 2 kb* katG *fragment carried on a plasmid.
The one point still unverified from a PDF is the crisprHAL 2 attribution — the
PeerJ full text would not load here — so confirm that before a manuscript
relies on it.)*

### The mechanism: which step is the bottleneck

Two numbers point the same way: GC-rich targets cut worse (−0.20), and a guide
that binds its target *more* strongly also cuts worse (−0.21).

The second is counter-intuitive — surely binding better is good? It identifies
the rate-limiting step. Before the guide can pair with its target, the DNA
double helix must be **pried open**. If grabbing on were limiting, stronger
binding would help; because stronger binding *hurts*, prying open must be
limiting — and GC-rich DNA is harder to pry open.

In the field's language: cutting here is **strand-invasion limited, not
hybridisation limited**.

### Two published findings confirmed, one shown not to transfer

- **A known penalty replicates in bacteria.** A G immediately after the PAM is
  the single most predictive flanking feature (−0.108), matching a ~12.6%
  penalty reported in human cells.
- **Published RNA-folding thresholds replicate** in a dataset their authors
  never used: 4.0% and 8.2% of guides fall below the two published cut-offs,
  and both groups are less active, as predicted.
- **A human-derived "optimal energy window" does not transfer.** Guides inside
  the range reported as optimal in human cells are, if anything, *less* active
  in bacteria (−0.133).

---

## Part 4 — The rule: when a new feature cannot possibly help

Most of the nine families moved the model by less than 0.01. For a week that
looked like failure. It is the most transferable result in the project.

### The rule

Every one of the 20 positions in the published matrix carries indicators for at
least three of its four possible letters — and if it is not three of them, it
must be the fourth. So **the matrix already contains the guide's sequence
completely**, losslessly.

That has a hard consequence. **Any new feature calculated purely from the
guide's sequence adds exactly zero information.** It can only restate, in a
different shape, something the model already had — however good the biology
sounds.

### The measurement

For each new feature: how much of it is recoverable from the guide's letters
alone? (If almost all, it was never new.)

| feature family | recoverable from sequence alone | what it actually gained |
|---|---:|---:|
| CRISPRoff binding energy | **98%** | +0.007 |
| DNA duplex stability | **85%** | +0.020 |
| RNA folding shape | 30% | +0.007 |
| DNA methylation sites | 14% | +0.004 |
| **flanking DNA** | **~0%** | **+0.080** |
| chromosome DNA abundance | ~0% | +0.045 |
| chromosome 3D packing | ~0% | +0.016 |
| gene transcription | ~0% | +0.015 |

Why the top two are so high is obvious in hindsight: duplex stability is
*defined* as a sum over neighbouring letter pairs, and the binding energy is
dominated by the same term. Importing them was importing the matrix back into
itself.

**What predicts whether a family helps is not the quality of the biology — it
is whether the feature describes anything outside the 20 letters.** Every
family that does helped; every family that is a sequence calculation failed.
`sgrna.diagnose --what redundancy` runs this check in seconds, and should be run
**before** building any future feature. This generalises well beyond CRISPR:
any paper adding "thermodynamic features" to a one-hot sequence model is
subject to it.

### The rule has an exception, and we found its limit

The rule says sequence-derived features add no *information*. It does not say
they are useless, because a model can be helped by the same information in a
better-shaped variable. Duplex stability is the proof: 85% recoverable, yet
worth +0.020 over raw GC content.

So "a physical reparameterisation can beat raw composition" is a real
possibility — and in October we tested it properly and it failed. **Family I**
describes the flanks by how readily they physically bend and deform
(dinucleotide bending and stiffness scales) rather than by letter composition,
including a feature for *coherent curvature* — bends spaced one helical turn
apart reinforce each other, which letter-counting cannot express.

| arm | ΔR² | Δρ |
|---|---:|---:|
| family I alone | +0.038 | +0.036 |
| family I, labels scrambled | −0.002 | −0.002 |
| family A alone | +0.087 | +0.078 |
| **family A + family I together** | **+0.088** | **+0.079** |

**On top of family A, family I is worth +0.0007 R²** against a seed-to-seed
spread of ±0.0032, and its best feature correlates with the model's remaining
errors at 0.028 — below the 0.035 noise level of features the model already
uses. The coherent-curvature feature, the one part that is not expressible as
composition, sits at 0.025.

So family A's windowed composition already spans what a physical shape scale
can say. **The "better-shaped variable" argument has now been tested twice and
gone both ways** — a success for duplex stability, a failure for flank shape.
Reporting only the success would overstate it, and the pair together is more
informative than either: what matters is not whether a variable is *physical*
but whether its shape adds anything to the encoding already present.

Family I is the seventh family shown to be unable to help, and the first where
the redundancy and residual diagnostics predicted the outcome before the
ablation confirmed it.

---

## Part 5 — Three of our own positive results, overturned by their own controls

This is the part to show a mentor. Each of these was a result we wanted.

### 1. The "supercoiling" effect is not supercoiling

The largest non-sequence effect came from a dataset measuring how twisted the
DNA is in each chromosomal region: **+0.047 R²**, with a clean scrambled
control. That would be a genuinely novel biological finding.

Then we ran the experiment's own **negative control** — the version of the
measurement with the biological part deliberately removed, which should contain
nothing. **It predicted cutting just as well.**

| what the model was given | ΔR² |
|---|---:|
| the whole family | +0.047 |
| raw read density, twisting stripped out | **+0.045** |
| *only* the twisting-specific ratios | **+0.004** |

**97% of the effect is how much DNA was present and sequenced in that region**,
not torsion. Anyone using a ChIP-style genomic track as a model feature needs
this control, and running it is not standard practice.

### 2. And the read-density effect is mostly a sequencing artifact

That left "a measured profile of local chromosomal DNA abundance predicts
cutting" — a weaker but still unreported claim, with the mechanism declared
open. In October we closed it, in three steps.

| block | columns | ΔR² | Δρ |
|---|---:|---:|---:|
| the 100 kb positional component | 16 | +0.013 | **−0.002** |
| the local within-bin component | 16 | **+0.036** | **+0.033** |
| **GC + 25-mer uniqueness, reference genome only** | **8** | **+0.038** | +0.030 |

- **It is not replication timing.** In a growing population, regions near the
  origin of DNA copying are physically present in more copies. That component
  is real, but gives **no ranking improvement at all** (−0.002 ρ) — a clean
  example of R² rising while the only number that matters does not.
- **It is GC bias and mappability.** Local read density tracks window-matched GC
  content at up to **−0.70**, and the untagged background track tracks 25-mer
  uniqueness at **−0.69**: sequencing is GC-biased, and reads map ambiguously in
  repeats. **Eight columns computed from the reference genome reproduce the
  entire 16-column measured gain** — so the GEO download was never needed, and
  this effect transfers free to any organism with a reference genome.
- **Against family A it nearly vanishes.** Paired across the same 15 folds,
  genome-computed GC and mappability add **−0.0002 ρ** on top of family A
  (8/15 folds, p = 0.84). The measured tracks retain **+0.0061 ρ** (15/15
  folds, p = 2.5e-05) — small, unambiguous, and not reproducible from sequence.
  That residual is the only genuinely unexplained part, about a tenth of what
  the standalone +0.045 implied.

**And the way our earlier reasoning was wrong is itself worth recording.** In
September we eliminated GC bias because *GC correlates only +0.03 to +0.06 with
the cut score* — true, but the wrong test. GC's correlation with the **read
density track** is an order of magnitude larger. A feature can route a signal
into a model without correlating with the outcome itself. **Test a candidate
cause against the mediator, not only against the outcome.**

### 3. Four feature families are one variable

Three chromosome-position families correlate with each other at 0.77–0.96,
because all three are built from counting sequencing reads and inherit the same
gradient. Individually they look like +0.047, +0.017 and +0.016 — three
findings. **Together they give +0.051**, barely more than the best alone. Two
thirds of the apparent evidence was one observation counted three times.

Part 5.2 makes it four: the abundance effect is largely flank composition, which
is family A. One explanation now does the work of four.

---

## Part 6 — Interpretability, measured rather than asserted

The source paper's framing is explainable AI, and every competing model claims
interpretability on the strength of one ranking of important features. We tested
whether that claim survives measurement. It does not.

**On this matrix, XGBoost's two standard importance methods rank its own top
features at ρ −0.02** — statistically independent. A mechanism read off one plot
is one plotting choice away from a different mechanism.

That is load-bearing rather than pedantic: it undercuts the explainable-AI
framing of the very paper this project extends, and it is, as far as we can
tell, **the first measured rather than asserted interpretability comparison in
this literature**. Five measurements replace the assertion — selection
stability, concentration, method agreement, faithfulness, direction consistency
(`RESULTS.md` §11).

Two consequences we acted on:

- **We switched the main model from XGBoost to LightGBM.** Same accuracy, far
  better behaved on exactly this measure, so the explanations can be trusted.
- **Feature importance is not evidence.** 20 of the *scrambled* flank features
  were still selected by the model as "important" (19.8 of the top 300). A
  feature being chosen as important proves nothing on its own — which is
  precisely how a paper ends up with a mechanism story built on noise.

---

## Part 7 — How much room is left: four independent ceilings

"Have we wrung the cloth dry?" has four separate answers, and they say
different things.

**1. The measurement itself caps any model at ~0.90.** Two independent screens
of the same guide library, 33,567 shared guides, agree with each other at only
**ρ 0.810**. Standard attenuation reasoning makes the highest correlation any
model could reach against a single observed score **√0.810 ≈ 0.90**. The two
screens used different enzymes, so part of that disagreement is real biology
rather than noise, making 0.90 a *conservative* floor on the ceiling. Two
supporting numbers: 85% of cut-score variance is *within* gene rather than
between genes, so the label really measures guide quality rather than which gene
was hit; and the published score agrees with crisprHAL's independent
re-derivation at ρ 0.945.

**2. The features we have are exhausted.** Correlate every feature against the
model's own errors: if a feature correlates with what the model got wrong, the
model is leaving information unused.

| features | count | max correlation with the errors |
|---|---:|---:|
| the ones the model used | 6,480 | **0.035** |
| the ones held back | 289 | 0.155 |

Nothing the model was given retains a correlation above 0.10. **It is not
underfitting; it is out of material.** A better model on these features will not
help.

**3. The model class is not the constraint.** Sixteen kinds of model on
identical splits and features span **0.116 ρ in total — 3.3 percentile points
in the guide you would actually pick.** A straight line comes within 0.03 of the
best, which is a statement about the biology: the effects mostly *add up* rather
than interacting. Nearest-neighbours fails, which is what you expect if the
answer is a sum of many small independent effects rather than a landscape with
neighbourhoods. Stacking sixteen models wins +0.002 at 48× the cost. And Part 3
showed a neural network on raw DNA merely ties.

**4. More data will not help.** The learning curve flattens hard: R² 0.241 at
4,000 guides, 0.296 at 13,880, 0.306 at 20,000, **0.307 at 26,000**. The last
6,000 guides bought +0.001. This assay saturates near 20,000 guides, so anyone
planning a bigger screen to get a better predictor should read that line first.

### So where is the remaining ~0.2?

Not in the label, not in the current features, not in the model, not in the
guide's representation, and — as of Part 5.2 — **not in chromosomal abundance
either**, which was September's nomination for "largest unexplained effect" and
is now mostly accounted for.

What is left is specific: the held-back column above, at 0.155, is read density,
nucleoid-protein dependence and origin distance. Part 5.2 shows most of that is
GC and mappability already inside family A, leaving **≈+0.007 R² of genuinely
unexplained local signal**. Identifying it needs data that is not on disk —
local protein occupancy, nucleotide pools, something assay-specific.

**That is a much better position than a noise ceiling.** The honest summary: the
*sequence* cloth is dry — rows, models and representations of the DNA letters
are all exhausted — and what remains is a small, real, local effect plus
whatever lives in features nobody has built yet.

---

## Part 8 — What is new here, and what is imported

### New

**1. A published feature matrix decoded and made reusable.** 99.0% of 6,232
columns now regenerate from any 20-letter guide, verified at 100% agreement on
1.17 million rebuilt values (Part 2). Without it the head-to-head could not have
been run.

**2. A rule for when a feature cannot possibly help, and a test costing
seconds.** The lossless-encoding argument plus its measurement, including the
exception and the limit of that exception (Part 4). **The most transferable
thing the project has produced.**

**3. Controls that overturn two of our own positive results**, and the
methodological lesson behind the second — test a candidate cause against the
mediator, not the outcome (Part 5).

**4. An effect shown to be computable from a reference genome alone.** The
chromosomal abundance gain needs no ChIP-seq download: eight genome-derived
columns reproduce it (Part 5.2). Practical, and it transfers to any organism
with a sequenced genome.

**5. Interpretability measured rather than asserted** — and found inconsistent
in the method the source paper relies on (Part 6).

**6. A ceiling for the whole field**, measured four independent ways, with the
remaining headroom localised rather than merely named (Part 7).

**7. The shape of the flank effect** — compositional rather than motif-borne,
long-range windows worth more than immediate context on a clean label,
downstream worth 3.4× upstream, and strand-invasion rather than hybridisation
limiting (Part 3). This is *why* a 0.3-second tree matches a GPU network.

**8. The flank effect shown to be a property of DNA, and to be what makes a
model portable between organisms.** +0.104 on a genome-wide screen with a
different nuclease under contiguous-arc CV, with the long-range-dominant shape
reproduced; and in cross-organism model transfer, ~90% retention in both
directions with the flank features contributing +0.074 to +0.150 of it
(Part 3). Paired with a measurement showing the one available other-organism
screen covers 4.3% of a chromosome — a dataset-suitability result that applies
equally to the published cross-species evidence.

**9. Parity with the state of the art at a thousandth of the compute**, from a
feature set whose source paper reported R² 0.249 — now confirmed against a
re-run of the competing model rather than its published number (Part 1).

### Imported, and openly so

Family **B** is the CRISPRoff authors' energy model, **C** is ViennaRNA, **F** a
regex for methylation motifs, **H** a genome scan; the chromosome-position
families read published GapR-seq, Hi-C and RegulonDB data with standard tools;
family **I**'s parameter scales are published dinucleotide tables. **Using them
is not a contribution.**

What is a contribution is that most of them were measured against a control and
found empty, with Part 4 explaining why in a way that generalises. **A negative
result with a mechanism behind it is worth more than a positive result without
one.**

---

## Part 9 — Where we corrected ourselves

Kept visible rather than quietly edited.

| what we said earlier | the correction |
|---|---|
| "cleaning the data is worth +0.003, i.e. nothing" | true for the basic model, **false for the flank features**, which gain +0.164 on clean data against +0.080 on the original. The noise was *hiding* the signal, not damping it |
| "the nearby letters are two thirds of the flank effect" | true on noisy data; on clean data the **distant averages are the larger half** (0.660 vs 0.618) |
| "we are ~0.054 behind the best published model" | that compared different guide sets; on identical guides it is **parity** |
| "the abundance effect is not GC bias; by elimination it is genuine chromosome state" | GC was tested against the *label* instead of against the read-density track. It **is** GC bias plus mappability, reproducible from the genome alone, and already inside family A |
| "DNA abundance is the single largest unexplained effect left" | it decomposes into replication timing (no ranking gain), GC/mappability (already in family A) and **≈+0.007 R² that is actually new** |

Four genuine bugs were also found and fixed:

1. An energy calculation used the wrong one of two similar quantities, so
   **zero** guides appeared to fall in the relevant range and the feature was a
   constant. After the fix: 63.7%. This exposed the failed-transfer finding in
   Part 3.
2. Seventeen features per window were near-identical copies of plain GC content
   (correlation 0.998). Replaced by four deliberately different ones — **348
   features instead of 408, four times the reach, and a better result.**
3. A single unlocatable guide formed a one-row test set and crashed the
   strictest cross-validation.
4. A restart check trusted a summary file rather than counting completed work,
   so a long run could permanently block itself.

Two data defects in family I's imported parameter table are recorded in
`features/i_shape.py`: one column with corrupt values, excluded, and two scales
that were not reverse-complement symmetric as published, symmetrised — without
which every feature would have depended on which DNA strand a guide happened to
target.

### How we guard against fooling ourselves

**Against the model cheating.** The screen places about 20 guides in every
gene, so a model given "where on the chromosome is this guide" can succeed by
memorising neighbourhoods. Three progressively stricter splits were added: by
whole gene, by 100,000-letter block, and by hiding one of five contiguous
chromosome arcs. Under strict splitting two families lost 42% and 36% of their
apparent value. **The flank features lost nothing** — +0.0867 under both the
loosest and the strictest scheme.

**Against reading meaning into noise.** Every family was re-run with its values
scrambled. Those land between −0.002 and +0.001, so anything below about ±0.004
is noise.

---

## Part 10 — What to do next, and how to frame the write-up

### Next

Two of September's four priorities are now closed — the crisprHAL re-run
(Part 1) and flank shape (Part 4) — and the third was substantially answered
(Part 5.2). What remains:

**Finish the transfer question before opening anything else.** It is now half
done — the effect is known to survive a new enzyme and a new organism, but only
as a feature family, and only as a headline number. Taking it to the same depth
as the *E. coli* work means three more steps, in order:

1. **Decompose family A in the new screens** exactly as Part 3 decomposed it
   here: the immediate ±10 nt against the 50–1000 nt windows. If the long-range
   half is again the larger one, then the *shape* of the effect transfers, not
   just its sign — and the shape is the mechanistic claim.
2. **Check the mechanism replicates**: do GC-rich targets still cut worse, does
   downstream still outweigh upstream, does the distance profile still peak at
   250–500 nt? A mechanism that holds in two organisms is a much stronger claim
   than a gain that holds in two organisms.
3. ~~The true model transfer~~ — **done** (Part 3): ~90% retention in both
   directions between *E. coli* and *C. rodentium*, with the flank features
   carrying the portability.
4. **For the remaining species question, find a genome-wide screen outside
   *E. coli*.** What is still unshown is whether the long-range gradient is
   *detectable* in another organism, as opposed to transferable into one. That
   needs a dataset nobody in this literature has, so it is a search rather than
   an analysis.

Only then:

4. **Identify the ≈+0.007 R² residual** in local abundance — the only genuinely
   unexplained signal left. Candidates: local protein occupancy, nucleotide
   pools, something assay-specific. Needs data not currently on disk.
5. **Write it up.** Parts 2–7 are complete, measured, and self-consistent.

### How to frame it

**Not** "a better sgRNA predictor". Even at parity that framing invites the one
comparison this project loses — against a full-time lab with a GPU.

Frame it as: **what determines whether a CRISPR guide works in *E. coli*, and
how to tell in advance whether a proposed explanation can possibly help.** Then
every result is load-bearing: a published dataset decoded and made usable by
anyone; a proof that seven well-motivated feature families could not have
helped, with a test that predicts it in seconds *and* a measured exception that
shows where the rule stops; two of our own positive results overturned by their
own controls; an effect shown to need no sequencing data at all; the first
measured interpretability comparison in this literature, which finds the source
paper's own method internally inconsistent; a ceiling showing the field has
~0.2 ρ of headroom and where it is not; the flank effect shown to be a smooth
gradient rather than a pattern — which is *why* a 0.3-second model keeps up;
and, as a by-product, parity with the state of the art at a thousandth of the
compute.

### Honest caveats to carry into the paper

- **One organism for the long-range claim, two nucleases.** The flank effect
  is shown to be enzyme-independent within *E. coli* (Part 3), but the species
  question is untested here: the only other-organism screen available covers
  4.3% of one chromosome. The competing line of work has cross-species evidence
  this project does not — though on inspection it is one confined screen plus a
  2 kb cloned fragment, and it belongs to crisprHAL 1 rather than to the
  crisprHAL 2 model benchmarked in Part 1 (see below).
- **The +0.011 margin over crisprHAL 2 rests on one seed of their model**, with
  their hyper-parameters rather than re-tuned ones. Consistent across all five
  folds, but small; "comparable performance at a thousandth of the cost" is the
  defensible phrasing, not "better".
- **Loose cross-validation flatters the chromosome-position families**, which is
  why they are reported under the strict scheme. The flank features survive it
  unchanged.
- **The ~0.90 ceiling rests on two screens that used different enzymes**, so
  part of their disagreement is biology rather than noise. It is a conservative
  floor on the ceiling, not the ceiling.
- **Part of the label's reproducible structure is systematic** — same library,
  same growth conditions, same locus effects — which sequence alone can never
  capture. So 0.90 bounds what is predictable in principle, not what is
  predictable from sequence.
- **The +0.007 residual is small enough to deserve replication** before any
  mechanistic story is attached to it.

---

## What exists now

- **31 Python modules, about 7,300 lines**, in `src/sgrna/`.
- **53 results files** in `results/` — every table here traces to one, and the early `superseded/` runs were deleted on 2026-10-02 once re-measured.
- **A 41-cell notebook** running the whole pipeline locally through Colab,
  saving progress after every step.
- **Three documents**: this report (merged with the former `NOVELTY.md` on
  2 October), `RESULTS.md` (all numbers with their settings), `FINDINGS_LOG.md`
  (chronological, including superseded claims).

---

## Sources

- [crisprHAL 2 — Better data for better predictions, PeerJ 2026](https://peerj.com/articles/20706/)
- [crisprHAL — A generalizable Cas9/sgRNA prediction model, Nat Commun 2023](https://www.nature.com/articles/s41467-023-41143-7)
- [Noshay et al. — Quantum biological insights into CRISPR-Cas9 sgRNA efficiency, NAR 2023](https://academic.oup.com/nar/article/51/19/10147/7279034)
- [DeepCC9 — An interpretable deep learning framework, Bioinformatics 2026](https://academic.oup.com/bioinformatics/article/42/7/btag483/8723703)
- [Improved prediction of bacterial CRISPRi guide efficiency, Genome Biology 2023](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03153-y)
- [BoostMEC — CRISPR-Cas9 cleavage efficiency through boosting models, BMC Bioinformatics 2022](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04998-z)
