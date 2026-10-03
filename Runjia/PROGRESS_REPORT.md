> **Synced copy — do not edit here.**
> The source of truth is `docs/PROGRESS_REPORT.md` in Runjia's RSI09 working
> folder, which also holds the companion files this report cites (`RESULTS.md`
> for every number with its settings, `FINDINGS_LOG.md` for the chronological
> record). Those are not in this repository, so references to them below are
> plain names rather than links.
> Synced 2026-10-03.

# RSI09 progress report

*Enhancing sgRNA efficiency prediction in CRISPR-Cas9 genome editing.*

Covers 8 September – 3 October 2026. Written to be read without a background in
molecular biology or machine learning: every technical term is defined in
Part 1, and each result states what it means in plain words before giving the
number.

The model built here is called **GuideGauge** throughout — provisional. The name
says what it does rather than what was found with it: it *gauges* how well a
candidate **guide** will cut. A tool should outlive any one result, so naming it
after the finding would be a mistake. Alternatives if it does not stick:
**CasGauge**, **LociScore**, **TERRAIN**.

Companion documents: `RESULTS.md` (every
number, with the protocol for each), `FINDINGS_LOG.md`
(chronological record, including every claim we later had to withdraw).

---

## Contents

- [Part 0 — Why this problem is worth working on](#part-0--why-this-problem-is-worth-working-on)
- [Part 1 — The words you'll need](#part-1--the-words-youll-need)
- [Part 2 — Where the project stands](#part-2--where-the-project-stands)
- [Part 3 — What made it possible: decoding the published dataset](#part-3--what-made-it-possible-decoding-the-published-dataset)
- [Part 4 — The one large effect, and what kind of thing it is](#part-4--the-one-large-effect-and-what-kind-of-thing-it-is)
- [Part 5 — The rule, and a blind test of it](#part-5--the-rule-and-a-blind-test-of-it)
- [Part 6 — Three of our own results, overturned by their own controls](#part-6--three-of-our-own-results-overturned-by-their-own-controls)
- [Part 7 — Controls: what we ran, and what we could not](#part-7--controls-what-we-ran-and-what-we-could-not)
- [Part 8 — The model: why LightGBM, and why the choice barely matters](#part-8--the-model-why-lightgbm-and-why-the-choice-barely-matters)
- [Part 9 — How much room is left, and can the limit be raised](#part-9--how-much-room-is-left-and-can-the-limit-be-raised)
- [Part 10 — What is new here, and what is imported](#part-10--what-is-new-here-and-what-is-imported)
- [Part 11 — Where we corrected ourselves](#part-11--where-we-corrected-ourselves)
- [Part 12 — Next steps and how to frame the write-up](#part-12--next-steps-and-how-to-frame-the-write-up)
- [Appendix A — Every feature set, grouped, with its data source](#appendix-a--every-feature-set-grouped-with-its-data-source)
- [Appendix B — Papers referred to](#appendix-b--papers-referred-to)

---

## Part 0 — Why this problem is worth working on

### What CRISPR is used for

CRISPR-Cas9 is a way of cutting DNA at a chosen place. Once cut, a cell's own
repair machinery either breaks the gene (useful for finding out what a gene
does) or pastes in a replacement sequence (useful for fixing or changing it).
Four kinds of use, roughly in order of how established they are:

1. **Research.** By far the largest use. If you want to know what a gene does,
   you break it and see what changes. Essentially all of modern genetics runs
   on this.
2. **Medicine.** Casgevy (exagamglogene autotemcel) was approved in late 2023
   for sickle-cell disease and beta-thalassaemia — the first approved CRISPR
   therapy — and has since been extended to younger patients. Many more are in
   trials.
3. **Agriculture and industry.** Crop traits; and **engineering bacteria** to
   manufacture things — insulin, fragrances, biofuels, drug precursors. This is
   a large existing industry, and it runs on *E. coli* and a handful of other
   microbes.
4. **Antimicrobials and diagnostics.** Using Cas9 to kill specific bacteria,
   and Cas enzymes as the detection step in diagnostic tests.

In every one of these, somebody has to choose a guide. Choosing badly wastes an
experiment; at scale it wastes a screen of tens of thousands of guides.

### Is a bacteria-specific result worth having?

This is a fair challenge, and we can now answer it with our own measurement
rather than an opinion. **We tested whether a model trained on bacteria
transfers to human cells, and it does not — at all.** In both directions the
correlation is slightly *negative*: −0.048 going one way, −0.017 the other
(Part 9). So nothing here should be sold as relevant to human gene therapy.

What it *is* relevant to:

- **Bacterial engineering**, which is the use-case above that actually runs in
  *E. coli*, and where guide choice is a daily practical problem.
- **Mechanism.** The reason the human transfer fails is itself a finding: the
  determinants of cutting are *different* in the two settings, and we can say
  which ones and by how much. A negative result with a measured boundary is
  more useful than a hedge.
- **Method.** The most transferable thing this project produced — the rule in
  Part 5 for telling in advance whether a proposed feature can help — is not
  about bacteria at all. It applies to any model built on DNA sequence.

So: the *predictions* are bacterial and should be advertised as such. The
*methods and the reasoning* are general. That split is worth being explicit
about in the write-up rather than leaving a reader to guess.

---

## Part 1 — The words you'll need

Skip if you know them. Nothing later re-explains them.

### The biology

**DNA** — a string in a four-letter alphabet: A, C, G, T. It normally exists as
a double helix, two strands zipped together, where A pairs with T and C with G.

**GC content** — the fraction of letters that are G or C rather than A or T.
G–C pairs are held together by three hydrogen bonds and A–T pairs by two, so
**GC-rich DNA is harder to pull apart**. This single fact explains a surprising
amount of what follows.

**CRISPR-Cas9** — molecular scissors. **Cas9** is a protein that cuts DNA but
does not know where.

**Guide RNA (sgRNA)**, or just **guide** — the address label. It carries a
20-letter sequence, and Cas9 cuts wherever the DNA matches it.

**Protospacer** — the 20 letters of DNA being targeted; the part the guide
matches.

**PAM** — a three-letter signal in the DNA immediately after the target, which
must read **NGG** (any letter, G, G) or Cas9 will not cut at all. A mandatory
suffix. Different Cas enzymes require different PAMs, which matters in Part 4.

**Flank / flanking DNA** — the DNA *surrounding* the 20-letter target. The guide
does not match it and Cas9 does not read it, so the obvious expectation is that
it is irrelevant. The central finding of this project is that it is not.

***E. coli*** — a gut bacterium, and the standard workhorse organism.
***C. rodentium*** — a mouse pathogen, a close relative; both are in the family
*Enterobacteriaceae*. How close they are turns out to matter (Part 4).

**Nucleoid-associated proteins** — bacteria pack their DNA with proteins such as
HU and H-NS. This is *not* the same as the **histones** and **nucleosomes** that
package DNA in humans, animals and plants. A nucleosome is a spool of DNA
wrapped around eight histone proteins, and it physically blocks Cas9. Bacteria
have no histones. This difference becomes important in Part 9.

**Chromatin** — the general term for DNA plus its packaging proteins.
**Chromatin accessibility** means how physically exposed a stretch of DNA is.

**Supercoiling** — DNA can be over- or under-twisted, like a coiled phone cable.
Over-twisted (positively supercoiled) DNA is harder to open. Part 6 is about a
dataset meant to measure this.

### The experiment the data comes from

**A screen** — an experiment testing tens of thousands of guides at once. In
bacteria, cutting the chromosome usually kills the cell, so you grow a mixed
population, sequence what's left, and a guide whose cells **disappeared** was a
guide that cut well.

**Cut score** — the efficiency number that comes out of a screen.

**The catch, which runs through everything:** the screen measures *survival*,
not cutting. Survival mixes cutting with DNA repair, with growth rate, and with
sequencing noise. The cut score is an indirect and noisy measure of the thing we
care about.

**Sequencing read depth** — how many times a given stretch of DNA was read. It
varies for uninteresting technical reasons (GC-rich fragments amplify less
efficiently; repeated sequence is ambiguous to map), and Part 6 is about
mistaking that variation for biology.

**Data cleaning / curation** — here, specifically: throwing away guides whose
cut score came from too few sequencing reads to be reliable. A guide read 5
times has a much noisier score than one read 500 times. The dataset we compare
against discards **45%** of the original guides this way. It does not change any
remaining guide's score; it removes the untrustworthy rows.

**CRISPRi** — a related technique using a *disabled* Cas9 (dCas9) that binds
without cutting, switching a gene off instead of breaking it. It measures a
different quantity, so CRISPRi datasets cannot be used here. Worth knowing
because most bacterial screens outside *E. coli* are CRISPRi (Part 4).

### The modelling

**Feature** — one number describing a guide: its GC content, the temperature at
which its DNA comes apart, and so on. A model never sees the guide itself, only
its features. Choosing features is most of the work.

**Feature set / family** — a group of features built from one idea and one data
source. This project has nine, listed in Appendix A.

**Model** — a program shown thousands of guides with their measured scores, which
then predicts the score of a guide it has not seen.

**Spearman correlation (ρ, "rho")** — how well the model gets the **ordering**
right. 0 is random, 1 is perfect. **This is the number to watch**, because
nobody needs a guide's exact efficiency: they have ten candidates and want the
best one.

**R²** — how much of the variation a model explains. Reported too, but it answers
a question nobody asks.

**Pick percentile** — the most honest reading. Take ten candidate guides, let the
model choose, and ask where that guide really falls. Blind choice gives the 50th
percentile; GuideGauge gives the **76th**; a perfect model would give the 91st.

**Cross-validation** — hide part of the data, train on the rest, test on the
hidden part, repeat. The discipline that stops you fooling yourself.

**Grouped cross-validation** — the stricter version used here for anything
positional: instead of hiding random guides, hide a whole contiguous chunk of
chromosome. Needed because the screen puts ~20 guides in every gene, so a model
told *where* a guide sits can score well by memorising "guides around here do
about this well" — which looks like skill and is not.

**Ablation** — the experiment used throughout: add one feature set, re-measure,
see whether the model actually got better.

**Permutation control** — re-run a feature set with its rows **shuffled**, so the
features are intact but no longer describe the right guide. Whatever that scores
is the set's noise floor. A gain inside the floor is not a gain. Part 7.

**Hyperparameters** — a model's dials, set before it sees data: how many trees,
how deep, how fast it learns. Not learned from the data. **Tuning** means
searching for good values.

**Quantum chemical features** — the previous paper's idea and this project's
inheritance. Physics calculations describing the electrons in each DNA letter:
how tightly neighbouring letters stack, how much energy their bonds hold.

**Gradient boosting (LightGBM, XGBoost)** — a family of model built from many
small decision trees, each correcting the previous one's mistakes. Fast, runs on
a normal laptop, and its features have names you can inspect.

**Neural network (CNN, BiGRU)** — a more flexible family that learns from raw
data. A **CNN** (convolutional neural network) scans for short recurring
patterns. Usually needs a GPU (a specialised chip) and is hard to interpret. The
best published competitor is one of these.

---

## Part 2 — Where the project stands

### The three reference points

| | what it is | Spearman ρ | guides |
|---|---|---:|---|
| Noshay et al. 2023 | the paper GuideGauge inherits its features from | 0.502 (Pearson) | 40,468 |
| **crisprHAL 2** (2026) | the best published bacterial model | **0.697** | 33,495 curated |
| **GuideGauge** | this project | **0.707** | 33,567 curated |
| **GuideGauge, tuned** | same, after searching its dials (below) | **0.718** | 33,567 curated |
| the apparent limit | see Part 9 | ~0.90 | — |

### The head-to-head, and an honest correction to it

Earlier numbers in this project were measured on different guides from the model
they were compared against, which makes the comparison meaningless. That was
fixed in two steps.

**Step one: same guides.** Because of Part 3, the published feature set can be
computed for *anyone's* guides. Running GuideGauge on crisprHAL 2's own 33,567
curated guides, with their label and their protocol, gives **ρ 0.7078**.

**Step two: their model, our folds.** Rather than trusting their published
number, we re-ran their actual model — their architecture, their
hyperparameters, their 48 training epochs, their code imported unchanged — on
our exact data splits. It scored **0.6971 ± 0.0067**. As a check that our
harness is faithful rather than flattering, we also ran it on their own split
and got **0.6940** against their published 0.695.

So GuideGauge is ahead by **+0.0107**, in 5 of 5 folds, paired t-test p = 0.0006.

**And here is the correction.** Neither model was tuned: ours ran on
hand-chosen settings, theirs on the settings its authors shipped. We asked how
much tuning could move our number, searching 12 configurations **inside each
training fold** (so the choice never sees the test data). Result:

| | ρ |
|---|---:|
| GuideGauge, default settings | 0.7078 |
| GuideGauge, tuned | **0.7181** |
| gain from tuning | **+0.0103**, 5/5 folds, p = 0.0015 |

**Tuning alone moves our score by as much as our entire lead.** So the +0.0107
margin cannot support "GuideGauge is better" — it is the same size as an effect we
know is available to whichever model gets tuned, and tuning theirs is not
affordable here (one configuration costs ~105 minutes of CPU for five folds, so
a twelve-point search is a day of compute).

**The defensible claim is parity.** Two models, built on completely different
principles, land within a hundredth of each other on identical data, and the
difference between them is no larger than the difference between tuned and
untuned versions of either. Writing it as a win would not survive a referee who
asks the question you just asked.

### What the two actually differ in: cost

Measured on one machine, CPU to CPU, with thread count pinned and five
repetitions (`results/timing_controlled.csv`):

| | GuideGauge | crisprHAL 2 |
|---|---|---|
| one fold | **17.1 s** (spread 0.7 s) | ~21 min |
| five folds, end to end | **1.4 min** | **105 min** |
| peak memory | 4.03 GB | 4.55 GB |
| features | named quantities you can look up | learned, inside a network |

So **about 74× less time for equal accuracy, at essentially the same memory.**
The saving is time, not footprint — do not describe the tree pipeline as light
on RAM, because it is not; it holds a 33,567 × 6,517 table in memory.

Two caveats that belong with that number. **crisprHAL 2 as published is
GPU-trained**, so this shows their architecture needs ~74× more *CPU* time, not
that a laptop beats a GPU. And **most of our 17 s is not the model** — 14 s of
it is a feature-selection step whose only output is a ranking. Swapping that
step's XGBoost for LightGBM makes the fold 6.6 s (≈190× against crisprHAL) with
accuracy indistinguishable (+0.0008).

**That swap was tried and then reverted**, for a reason worth recording. The
headline is fine under it — the head-to-head goes 0.7078 → 0.7092 — but the
*baseline* moves from 0.2937/0.5278 to 0.2897/0.5251, and that baseline is this
project's reproduction check: the line that tells a reader our harness
reproduces the published figure exactly. Trading a verifiable anchor for ten
seconds is a bad deal, so the default stays XGBoost and LightGBM is available as
a setting (`config.SELECTOR`) for exploratory sweeps where the anchor does not
matter.

> **This figure has been wrong twice, so it is worth saying how.** An earlier
> version claimed "0.3 seconds, ~25,000× faster". The 0.3 s was real but was the
> *final model fit alone*, excluding the selection step that produces its
> inputs — and it was then compared against crisprHAL's entire five-fold
> training rather than against one of ours. Uncontrolled measurements on this
> machine also scattered between 13 and 116 seconds per fold, because it runs
> out of RAM and swaps; that is why the numbers above come from a pinned,
> repeated run instead. **State what a measured number measures: a fit time is
> not a training time.**

### Reading 0.71 practically

Given ten candidate guides and GuideGauge's favourite, that guide lands at the
**76th percentile** of true efficiency — against the 50th for picking blindly
and the 70th for the model we started from. Real, useful, and not a
transformation.

---

## Part 3 — What made it possible: decoding the published dataset

### The problem

The project began from a 2023 paper describing 13,880 guides with 6,232 numbers
each, mostly quantum chemistry. Two things made it nearly unusable by anyone
else:

- **it never recorded the DNA sequence of any guide**, only derived numbers; and
- **5,887 of the 6,232 columns were named `V1`, `V2`, `V3`…** with no
  description.

So the method could only ever be applied to guides its authors had already
processed.

### What we found

**The sequence is hidden in the numbers.** One column records the electron count
of the DNA letter at each position. Across all 13,880 guides it takes exactly
four values — 42, 48, 50, 56 — and those are the valence-electron counts of C,
T, A and G. That column *is* the sequence in a different alphabet. Reading it
position by position reconstructs every guide.

**It checks out against the real genome.** 13,879 of 13,880 decoded (one row has
a missing value); every one found in the *E. coli* genome with the required NGG
PAM in the right place; 13,877 appearing exactly once. As an independent check,
a standard method for locating where bacterial DNA copying starts put it at
position 3,923,620 — the textbook value is near 3,923,800.

**The quantum tables are recoverable exactly.** Each quantum column is a fixed
lookup: a given short run of letters always gives the same number. All 4 single
letters, 16 pairs, 64 triples and 256 quadruples appear in the data **with no
contradictions anywhere**. So the parameterisation is recovered, not
approximated, and can be applied to DNA the paper never touched — which is
exactly what Part 4's feature set does.

**The anonymous columns were identified.** 5,853 of the 5,887 `V####` columns are
"is there letter X at position Y" indicators. None ambiguous, and verified
against **1.17 million rebuilt values at 100% agreement**.

### What the decoding also tells us: 6,232 columns, about 20 variables

The decoding is usually described as a reuse result. It is also a statement
about what the matrix *contains*, and that turns out to matter more.

A 20-letter guide holds at most 40 bits of information. The matrix spends 6,232
numbers describing it, in three overlapping ways: 5,853 "is letter X at
position Y" indicators, 316 quantum-chemistry columns that are fixed per-letter
lookups, and a few summaries (GC, melting temperature) that are weighted sums
of the same positions.

We measured this directly. Two different ways of picking the "best" 300 columns
agree on only **114 of them** — yet both cover **all 20 positions of the guide,
100% overlap**, at about nine columns per position. They are not choosing
different information. They are choosing different spellings of the same twenty
letters.

Three results reported later are consequences of this, and should be read as
one finding rather than three: sixteen model classes span only 0.116 ρ (Part 8),
quadrupling the feature cap changes nothing (Part 8), and nothing the model
holds correlates with its own errors above 0.035 (Part 9). All three say there
are few independent dimensions here to exploit.

**So "6,232 features" is not a rich description.** It is an elaborate spelling
of a 20-letter word. That is exactly why the flanking-DNA features mattered:
they were the first thing added that described DNA *outside* the guide, and so
the first thing that could add information rather than restate it.

### Running the model backwards: what does it think a perfect guide is?

Because the features are computed from the guide, the model can be inverted —
searched for the sequence it scores highest.

| | sequence | predicted | GC |
|---|---|---:|---:|
| the model's ideal guide | `ACTGCACAAAAGATGTCTTT` | 38.76 | **35%** |
| the model's worst guide | `AAAAAAATTGCCCCCCGGGG` | −14.14 | 55% |
| *real guides, 1st–99th percentile* | — | *0.98 – 41.73* | *mean 51%* |

Two sanity checks pass: the ideal guide scores *inside* the range really
observed rather than off the end of it, and it is **GC-poor** — which is the
strand-invasion mechanism from the previous section, recovered by a completely
independent route.

**And that guide exists nowhere in the *E. coli* genome.** This is the honest
limit of inversion here: nobody picks a guide from all possible sequences, they
pick one from the few valid sites inside the gene they want to cut. The global
optimum is unreachable by construction, so inversion is useful for *diagnosis*
— reading back what the model learned, in a form a biologist can check — and
not for design. The design question is the constrained one: of the guides
available in my target, which is best? That is ranking.

### Why it matters

**99.0% of the published dataset — 6,169 of 6,232 columns — now regenerates from
a guide's 20 letters alone.** A resource that worked only for its authors works
for anyone, on any bacterial CRISPR screen.

Everything downstream depends on this: the head-to-head in Part 2, the
cross-organism tests in Part 4, and the human boundary in Part 9 are all only
possible because the representation became portable.

---

## Part 4 — The one large effect, and what kind of thing it is

This part makes a single argument in six steps: the DNA *around* a cut site
predicts cutting; the signal is a smooth gradient rather than a pattern; it
belongs to the DNA rather than to the enzyme; it crosses between organisms;
it stops at the boundary of the kingdom; and it has a mechanism we can name.

### 1. The surrounding DNA predicts cutting, and it is the largest effect found

The guide does not match the flanking DNA and Cas9 does not read it, so the
expectation is nothing. Instead, features describing the flanks are worth
**+0.080 Spearman** on the original data and **+0.164** on the cleaned version —
larger than every other idea tested in this project combined.

The features are deliberately simple: for windows of 50, 250, 500 and 1,000
letters on each side, the GC fraction, the purine fraction, the longest run of a
repeated letter, and averages of the recovered quantum tables. Plus the identity
of each of the 10 letters immediately either side.

### 2. The signal is a gradient, not a pattern — and that is testable

If the flank effect were a short recurring motif, a convolutional neural network
should find it better than hand-computed averages do, because that is exactly
what CNNs are built for. We built one with the same architecture family as the
best published competitor and gave it raw DNA letters.

On the **same rows and the same label**, adding ±100 letters of raw flanking
sequence to the network is worth **+0.019 ρ**, while the hand-computed windowed
composition is worth **+0.080** — about **four times as much**, from the same
region of DNA.

> **A correction to an earlier version of this report**, which paired the
> network's +0.019 against the +0.164 measured on the *cleaned* dataset and
> called it "four times". Those two numbers come from different data, so the
> ratio was meaningless. **Compared like with like — both on the original
> 13,880 guides — it is +0.019 against +0.080, a factor of 4.2.** The
> conclusion is unchanged; the arithmetic supporting it was wrong.

So the effect has no motif to find. It is something closer to "how GC-rich are
the next 500 letters" — a smooth average that a windowed mean computes exactly
and a 5-letter pattern detector has to approximate badly. **We used the right
instrument for the shape of the signal, and that is the whole of GuideGauge's
advantage.**

Two more facts about the shape: it **peaks at 250–500 letters** and has faded by
1,000 (so there is no point looking further out — this question is closed), and
**downstream matters about 3.4× more than upstream**.

### 3. It is a property of the DNA, not of this particular enzyme

A result from one enzyme in one organism could be an idiosyncrasy. The crisprHAL
project ships four screens, and three are usable:

| screen | organism | enzyme | guides | base | + flanks | gain |
|---|---|---|---:|---:|---:|---:|
| WT-SpCas9 | *E. coli* | SpCas9 | 33,567 | 0.544 | **0.707** | +0.163 |
| eSpCas9 | *E. coli* | eSpCas9 | 59,489 | 0.685 | **0.789** | +0.104 |
| TevSpCas9 | *C. rodentium* | TevSpCas9 | 25,210 | 0.704 | **0.764** | +0.060 |

All three under the strict grouped cross-validation, all gains winning 5 of 5
folds (p = 7×10⁻⁶ and 3×10⁻⁵ for the two new ones).

The fourth screen, **TevSaCas9**, was deliberately excluded. Its enzyme requires
a different PAM (NNGRRT instead of NGG), so only 45% of its target sites have the
NGG our pipeline assumes, and the position labels for the other 55% would be
silently shifted. Our earlier phrasing — "it would have produced numbers" — was
too cryptic. **What it means: the pipeline would have run to completion without
any error and returned plausible-looking results that were measuring the wrong
positions.** That is more dangerous than a crash, and it is why the screen was
left out rather than included with a caveat.

> **A limitation we should state rather than let a reader find.** eSpCas9 is
> WT-SpCas9 with three point mutations, and the two screens use the **same guide
> library** — 100% of the same target sequences. So "two different enzymes agree"
> is a much weaker statement than it sounds: they are near-identical enzymes on
> identical DNA. The genuinely different enzyme is TevSpCas9, which is a fusion
> protein needing its own extra recognition motif — and it is also in a
> different organism, so that arm confounds the two variables. **The cross-enzyme
> evidence is real but narrow.**

### 4. A model trained in one organism ranks another organism's guides

The tests above all train and test within one screen. A stronger question is
whether a *trained model* crosses. Setting it up as a 2×2 needs one enzyme per
organism; we use *E. coli* WT-SpCas9 and *C. rodentium* TevSpCas9, which is the
only pairing the available data allows.

**Spearman ρ, GuideGauge with flank features, every cell on the full target screen:**

| | tested on *E. coli* | tested on *C. rodentium* |
|---|---:|---:|
| **trained on *E. coli*** | **0.708** | 0.700 |
| **trained on *C. rodentium*** | 0.628 | **0.764** |

Read the off-diagonal against the diagonal below it: going *E. coli* →
*C. rodentium* keeps **92%** of what a locally-trained model achieves; going
back the other way keeps **89%**. The guides share no sequence and the organisms
share no chromosome, so this is genuine transfer.

**The flank features are the portable part.** Without them the same table reads
0.626 and 0.478 — the inherited published representation transfers poorly, and
nearly all the retention comes from the flank block. Crossing out of
*C. rodentium*, 88 long-range columns beat 260 short-range ones (+0.107 against
+0.058), so it is specifically the *gradient* that travels.

**Why the enzymes differing between cells is not fatal, and the note it needs.**
Ideally both cells of a row would use the same enzyme. They cannot: no screen
uses the same nuclease in both organisms. The reason this is tolerable is
measurable rather than assumed — within *E. coli*, a model trained on
WT-SpCas9 and tested on eSpCas9 retains **87%**, which is the same band as the
cross-organism retentions (89–92%). So changing enzyme and changing organism
cost about the same amount, and neither dominates the table. That said, the
WT→eSp arm is **not** a generalisation test (those two screens share all their
guides, so every test sequence was in training) and must not be quoted as one.
It bounds enzyme sensitivity, nothing more.

### 5. The species question is still open, and the reason is the data

Both new screens are *Enterobacteriaceae* — *E. coli* and *C. rodentium* are in
the same family, so **two species this close is a weak test of generality.**
Worse:

**The *C. rodentium* screen covers 4.3% of one chromosome.** Its 25,210 guides
sit in a 229 kb span at 110 guides per kb — which is the authors' stated design
(a 236 kb fragment), not a defect we found. Inside a window that narrow the
long-range gradient cannot even be *seen*: flank GC correlates −0.022 with the
label there, against +0.159 and +0.138 in the two *E. coli* screens. So that
screen can confirm a transferred model but cannot discover the effect.
Discovery and verification have different data requirements.

**And no better dataset exists.** We searched for a genome-wide Cas9 *cutting*
efficiency screen in a more distant bacterium — a different phylum, say
*Bacillus* (Firmicutes) or *Mycobacterium* (Actinobacteria). There isn't one.
What exists in those organisms is **CRISPRi** — disabled Cas9 that silences
rather than cuts — which measures a different quantity and cannot substitute.

The competition is in the same position, and this is worth getting right because
we previously misattributed it. The cross-species claim belongs to **crisprHAL 1**
(*Nat Commun* 2023), not to crisprHAL 2 (a data-curation paper). Their two
non-*E. coli* validations are a **236 kb fragment** in *C. rodentium* and
**296 guides in 2 kb of *S. enterica* DNA cloned onto a plasmid inside *E.
coli***. Both are confined fragments; one is foreign sequence in an *E. coli*
cell.

**So the honest position is not "they solved cross-species and we did not". It is
that nobody has the dataset the question needs, and this project is the one that
measured why the existing substitutes cannot stand in for it.** Testing a
different phylum requires generating a genome-wide screen there — a wet-lab
project, not an analysis one.

### 6. The mechanism: which step is the bottleneck

Two correlations point the same way. GC-rich *targets* cut worse (−0.20 to
−0.15), and a guide that binds its target *more* strongly also cuts worse
(−0.21).

The second is the informative one. Before a guide can pair with its target, the
DNA double helix has to be **pried apart**. If grabbing on were the slow step,
stronger binding would help. Because stronger binding *hurts*, the slow step
must be the prying apart — and GC-rich DNA, with its extra hydrogen bond per
pair, is harder to pry apart.

In the field's language: cutting here is **strand-invasion limited, not
hybridisation limited.** This replicates in all three bacterial screens, and is
*strongest* in *C. rodentium* (−0.379) — so unlike the flank gradient, this one
does not depend on a screen's genomic span.

---

## Part 5 — The rule, and a blind test of it

Six of the nine feature sets improved the model by less than 0.01. For most of
the first week that looked like failure. It is the most transferable result the
project produced.

### The rule

Every one of the 20 target positions in the published dataset carries indicators
for at least three of its four possible letters — and if a position is not three
of them, it must be the fourth. So **the dataset already encodes the guide's
sequence completely and without loss.**

That has a hard consequence. **A feature computed purely from the guide's
20 letters adds exactly zero new information.** It can only restate, in a shape
the model may find easier, something already present. However good the biology
sounds.

### Testing it the right way round

Explaining six failures after the fact with a rule invented after seeing them is
weak. So here is the rule applied **as a prediction**, with the reasoning written
out before the results column is read.

The prediction is mechanical, in two questions:

1. **Can this feature be computed from the 20 target letters alone?** If yes, it
   carries no new information, and should give at most a small gain from being a
   more convenient shape.
2. **If it needs outside information, is that information actually about the
   locus?** If it comes from a measurement whose variation is technical rather
   than biological, it should give nothing useful.

| feature set | computable from the 20 letters? | **predicted** | **measured Δρ** | prediction correct? |
|---|---|---|---:|---|
| `b_energy` — binding energy | yes, it is a sum over letter pairs | nothing new; small shape gain | +0.007 | ✓ |
| `d_mechanics` — duplex stability | yes, also a sum over letter pairs | nothing new; small shape gain | +0.020 | ✓ |
| `c_folding` — RNA folding | yes — the guide folds on its own | nothing new; small shape gain | +0.007 | ✓ |
| `f_methylation` — methylation motifs | yes, it is a text search | nothing new; small shape gain | +0.003 | ✓ |
| `i_shape` — DNA bendability | **no** — it describes the flanks | genuinely new ⇒ should help | +0.036 alone, **+0.001 on top of flanks** | ✗ (see below) |
| `h_offtarget` — near-matches elsewhere | needs the whole genome | new information ⇒ should help | +0.000 | ✗ (see below) |
| `a_flank` — flanking composition | **no** | genuinely new ⇒ should help | **+0.080** | ✓ |
| `d_supercoiling` — chromosome position | **no**, needs a measurement | new ⇒ should help | +0.045 | ✓ but see Part 6 |
| `e_nucleoid` — 3D packing | **no**, needs a measurement | new ⇒ should help | +0.016 | ✓ but see Part 6 |
| `g_transcription` — gene activity | **no**, needs a measurement | new ⇒ should help | +0.015 | ✓ but see Part 6 |

**Seven of ten predictions correct, and the three misses are each informative
rather than noise.**

- **`i_shape` is the rule's limit, and the most interesting failure.** DNA
  bendability is genuinely outside the 20 letters, so the rule predicted it
  should help, and alone it does (+0.036). But stacked on the flank features it
  adds **+0.001** against a seed-to-seed spread of ±0.003. The reason is that
  the rule as stated is too coarse: being outside the 20-mer is necessary but
  not sufficient — the information also has to be outside *everything already
  in the model*, and position-resolved composition over four window sizes
  already spans what an average over 2-letter steps can express. **The refined
  rule: ask whether a feature is computable from the features you already
  have, not just from the sequence.**
- **`h_offtarget` fails for a different reason** — the information is real but
  the effect isn't there. Its gain (+0.0003) is *below its own noise floor*
  (Part 7), so it is indistinguishable from shuffled data. Near-matches
  elsewhere in the genome simply do not affect this measurement.
- **The three positional sets passed the rule and then failed a different
  test** — not redundancy but causation, which is Part 6.

### Why `d_mechanics` is the exception worth reporting

Among the "computable from the 20 letters" group, `d_mechanics` gains +0.020 —
small, but three times the others and well clear of its floor. A physical
reparameterisation of information already present *can* help a tree, because it
says *where* along the guide the duplex is weak rather than how GC-rich it is
overall, and a tree would otherwise have to learn that from scratch.

So "a better-shaped version of the same information" is a claim that has now
been tested twice and gone both ways: it won for `d_mechanics` (+0.020 over raw
GC) and lost for `i_shape` (−0.000 against the composition it was derived from).
Reporting both is more honest, and more useful, than quoting only the success.

### The practical upshot

Run the redundancy check **before** building a feature. It takes seconds
(`python -m sgrna.diagnose --what redundancy`) and would have saved four of the
nine feature sets from being built at all. The measured recoverabilities:
binding energy **98%** recoverable from the guide alone by a straight line;
duplex stability **85%**; RNA folding 30%; methylation motifs 14%; everything
that helped, **≈0%**.

---

## Part 6 — Three of our own results, overturned by their own controls

### 1. The "supercoiling" effect is not supercoiling

Our largest non-sequence gain (+0.047 R²) came from a dataset measuring DNA
twisting across the chromosome. That would have been a novel biological finding.

Then we ran the source experiment's own **negative control** — a version of the
measurement with the biological part deliberately removed, which should contain
nothing. **It predicted cutting just as well.**

| what the model was given | ΔR² |
|---|---:|
| the whole feature set | +0.047 |
| read depth only, with the twisting part removed | **+0.045** |
| **only** the twisting-specific part | **+0.004** |

97% of the effect is **how much DNA was sequenced in that region**, not how
twisted it is. **The write-up must not say "supercoiling predicts sgRNA
efficiency."** It says "a measured profile of local DNA read depth does".

### 2. And that read-depth effect is mostly a sequencing artefact

Having renamed the effect honestly, we asked what it actually is. Two
candidates, both computable from the reference genome with no experiment at all:
how **unique** the local sequence is (repeated DNA is ambiguous to map), and its
**GC content** (GC-rich fragments amplify less efficiently during library
preparation).

They correlate with the read-depth tracks at up to **−0.70**. And a feature set
of **8 columns computed from the genome alone reproduces the entire gain**
(+0.038 against +0.036).

> **What those 8 columns are**, since they carry a lot of weight. Two
> quantities, each at four window sizes (100, 500, 2,000 and 10,000 letters
> centred on the guide): (a) the fraction of positions in the window whose
> surrounding 25-letter sequence occurs **exactly once** in the *E. coli*
> genome — 97.35% of the genome passes this, and the failures are repeats; and
> (b) the GC fraction of the window. 2 × 4 = 8. No sequencing data is involved;
> both are read straight off the reference genome.

So the GEO download, the ChIP track, the untagged control and the rifampicin arm
were all unnecessary. Convenient, too: a genome-only feature transfers to any
organism with a reference sequence, at no data cost.

**And then it collapses as well.** Stacked on the flank features, the
genome-derived version adds **−0.0002** (8 of 15 folds, p = 0.84) — nothing. The
sequence-intrinsic explanation of read depth *is* windowed flank composition
under another name.

What survives is small and real: the *measured* read depth still adds **+0.0061
ρ, winning 15 of 15 folds** (p = 2.5×10⁻⁵), and that part is **not** reproducible
from the genome. About a tenth of what the standalone +0.045 implied. Our
earlier claim that this was "the single largest unexplained effect left" was
overstated and has been withdrawn.

### 3. Four feature sets are one variable

This is why grouping matters. Four of the nine sets trace back to the same
underlying quantity:

**The position group** — `d_supercoiling` (read depth), `e_nucleoid` (3D
packing from Hi-C), `g_transcription` (gene activity), and the genome-derived
mappability-and-GC block. Their leading features correlate with each other at
0.77–0.96, because Hi-C contacts and ChIP coverage are *both* read counts and
both inherit the same GC and mappability bias. Individually they look like
+0.047, +0.017 and +0.016 — three findings. **Together: +0.051.** Two thirds of
the apparent evidence was one observation counted three times. And all of it
then collapses into the flank features.

**The thermodynamic group** — `b_energy`, `c_folding`, `d_mechanics` all
describe how hard the duplex is to open, and all correlate 0.90–0.96 with a
plain GC-content column the dataset always had.

Appendix A groups all nine sets this way, which is how they should appear in the
paper: **two or three distinct ideas, not nine.**

---

## Part 7 — Controls: what we ran, and what we could not

Two different kinds of control get confused, and the distinction matters.

### Permutation controls — now run on every family that has a gain

A permutation control re-runs a feature set with its rows shuffled: the columns
are intact and just as numerous, but no longer describe the right guide.
Whatever that scores is the set's noise floor.

**Until 3 October only four sets had one**, and this project was quoting a single
floor measured on `a_flank` as though it applied to all of them. It does not —
the floor depends on how many columns a set contributes and how easily they
survive selection. The missing five were run, and the delta recomputed **paired
fold by fold** (previously a set run on three seeds was being differenced
against a baseline run on five, mixing the effect with seed-to-seed spread):

| feature set | real Δρ | its own floor | ratio | verdict |
|---|---:|---:|---:|---|
| `a_flank` | +0.0795 | −0.0009 | 88× | clear |
| `d_mechanics` | +0.0201 | +0.0001 | 153× | clear |
| `b_energy` | +0.0074 | −0.0001 | 92× | clear |
| `c_folding` | +0.0074 | +0.0010 | 7.7× | clear |
| `f_methylation` | +0.0034 | +0.0011 | 3.0× | **marginal** |
| `h_offtarget` | +0.0003 | +0.0004 | 0.8× | **indistinguishable from noise** |

Two new conclusions: **`h_offtarget` should be reported as a null result**, not a
tiny positive one; and **`f_methylation` is only 3× its own floor** and should
be described as marginal.

One sobering detail from the `a_flank` control: **20 of its shuffled columns
were still chosen by the model as "important"**. A feature being selected
proves nothing on its own.

### Negative controls from the source experiment — available for exactly one set

A permutation control asks "is this better than noise?". It cannot ask "is my
biological interpretation right?" — for that you need a control *built into the
original experiment*, like the untagged no-antibody track that demolished the
supercoiling story in Part 6.

**This kind of control exists for only one of our data sources.** Being explicit
about it:

| data source | mock/negative control available? |
|---|---|
| GapR-seq (supercoiling) | **yes** — untagged, no-antibody, and rifampicin arms. Used, and it overturned the result. |
| Hi-C (3D packing) | no mock arm in the published data |
| RegulonDB / PRECISE-1K (transcription) | not that kind of data — annotation, not a measurement with a control |
| the reference genome (flanks, methylation, off-target, shape) | not applicable — no experiment to control |

So the answer to "has that control been run on everything?" is: **the
permutation control now has, and the source-experiment control could only ever
apply to one set, where it was run and did change the conclusion.** For
`e_nucleoid` the equivalent check was indirect — we showed its signal correlates
0.77–0.96 with the read-depth tracks that the untagged control had already
discredited, which is weaker evidence and should be described as such.

---

## Part 8 — The model: why LightGBM, and why the choice barely matters

### Sixteen models, tested on identical data

| | ρ | | ρ |
|---|---:|---|---:|
| stacked combination | **0.611** | PLS | 0.575 |
| LightGBM | 0.609 | small neural net | 0.572 |
| XGBoost | 0.608 | extra trees | 0.556 |
| CatBoost | 0.601 | random forest | 0.550 |
| hist gradient boosting | 0.600 | LightGBM, robust loss | 0.500 |
| support vector machine | 0.586 | nearest neighbours | 0.495 |
| ridge regression (a straight line) | 0.577 | | |

The whole range is 0.116, and in terms of the guide you would actually pick,
**3.3 percentile points**. Three readings:

- **A straight line gets within 0.03 of the best.** That is a statement about
  the biology: the effects mostly add up rather than interacting in complicated
  ways.
- **Nearest neighbours fails**, so guides with similar features do not have
  similar scores — consistent with many small independent effects rather than a
  landscape with neighbourhoods.
- **Combining models wins by +0.002 at 48× the cost**, and its internal weights
  give the random forest a weight of −0.003. Not worth adopting.

### So why LightGBM?

**Not for accuracy — for interpretability.** It ties XGBoost on score (0.609 vs
0.608) and behaves far better when asked to explain itself:

| | LightGBM | XGBoost |
|---|---:|---:|
| do its two importance methods agree? | **+0.44** | **−0.02** |
| same features chosen on a different split? | 0.72 | 0.49 |
| features needed for half the importance | 151 | 330 |

Since the accuracy is identical, there is no cost to preferring the model whose
explanation is reproducible. But the row that reads worst — XGBoost's two
standard ways of ranking its own features agreeing at **−0.02**, i.e. not at
all — raises a question sharp enough to deserve its own test.

### What it means that the two models disagree about their own features

Put the facts together and they look contradictory. XGBoost's two importance
methods rank its features independently. LightGBM's agree. The two models pick
substantially **different** features. And they score **the same**. So either one
of them is wrong about the biology, or something else is going on.

We tested it (`src/sgrna/attribution.py`, `results/attribution*.csv`). Let each
model choose its own 300 features on the same data, then compare the top 50:

| | result |
|---|---|
| columns both models chose | **14 of 50** |
| overall correlation of the two importance rankings | +0.69 |
| held-out score | 0.627 vs 0.622 — **tied** |

So they really do choose differently, and it really does cost nothing. Two
explanations were on the table, and the test separates them.

**What is not the explanation: "they found different biology."** For each column
only one model picked, we asked how well its closest counterpart in the *other*
model's set stands in for it. Median absolute correlation **0.456**, against
**0.069** for randomly chosen columns. So the two sets sit in the same
correlated neighbourhood of the matrix — nowhere near random — which is not
what two genuinely different findings would look like.

**But the simple version of redundancy is not the whole explanation either.**
Only 22% of the unshared picks have a counterpart above 0.7. If every column
had a near-twin, substitution would be one-to-one and that number would be
high. It is not, and the honest reading is that a tree does not need a single
substitute for a dropped column — it can rebuild the same function from several
weakly-correlated ones. Pairwise correlation therefore *understates* how
replaceable a column is, and this test can only put a floor under it.

**What the test does show, and we did not anticipate: the two models prefer
different *kinds* of column, and that preference is algorithmic, not
biological.** Dividing each model's importance by what the column describes:

| kind of column | XGBoost | LightGBM |
|---|---:|---:|
| target: position/letter indicators (binary, 0 or 1) | **49%** | 18% |
| target: quantum descriptors (continuous) | 19% | **33%** |
| flanking-DNA windows (continuous) | 9% | **25%** |
| flanking DNA, nearest 10 letters | 11% | 17% |
| other published columns | 12% | 6% |

The disagreement is almost entirely along one axis: **XGBoost leans on the
binary indicators, LightGBM on the continuous columns.** That is a known
consequence of how the two build trees — LightGBM sorts continuous values into
histogram bins and grows leaf-by-leaf, which makes a continuous column cheap to
split on repeatedly, while XGBoost's depth-wise growth at these settings finds
the thousands of binary indicators competitive.

So a large part of what an importance ranking reflects is **the splitting
algorithm's affinity for a column's data type** — not how much that column
matters to the outcome. That is the clearest single reason not to read mechanism
off an importance plot.

**And it explains the −0.02 as well.** Gain and permutation importance answer
different questions. Gain asks how much splits on a column improved the fit,
and credits whichever of several interchangeable columns a tree happened to use
first. Permutation asks how much the prediction degrades when that column is
destroyed — and if the others can rebuild it, the answer is "barely", no matter
how much gain it was credited with. The two therefore diverge most when
importance is spread thinly across many interchangeable columns, which is
exactly XGBoost's situation here: 330 columns for half the importance, most of
them binary indicators each carrying a sliver. LightGBM concentrates on 151
mostly continuous columns that are individually harder to replace, so
destroying one does measurably hurt, and the two methods line up.

**The decisive evidence that neither has found "the true features" is
faithfulness.** Drop each model's top 20 features and refit: it costs barely
more than dropping 20 at random (+0.008 for LightGBM, +0.003 for XGBoost). If
the top of the ranking were load-bearing, removing it would hurt. It does not.

Put plainly: **method agreement measures how concentrated and individually
irreplaceable an attribution is. It does not measure whether the attribution is
right.** LightGBM's +0.44 makes its ranking *reproducible*, which is a real
advantage and the reason to use it — but reproducible is not the same as
correct, and the earlier wording in this report ("so the story is stable",
implying the story is therefore trustworthy) claimed too much.

### What to trust instead

Three kinds of evidence survived every model we tried, and they are what the
biological claims in this report actually rest on:

1. **Group-level ablations defined by a hypothesis before fitting** — add the
   flanking-DNA block, measure, compare against its own shuffled control. The
   answer does not depend on which model or which columns.
2. **The sign and size of a single named quantity** — target GC against score
   is −0.201, flank GC is +0.159. Anyone can recompute these in one line; no
   model is involved.
3. **Transfer** — does the relationship still hold in a different screen,
   enzyme, or organism? Part 4. A spurious attribution does not survive this.

None of the biology in Part 4 came from an importance ranking, which in
hindsight was lucky rather than principled. It should be stated as policy in
the paper.

### A consequence for the paper this project extends

Noshay et al. read their biological conclusion — that quantum-chemical
properties at the 3′ end of the guide are what matter — off the importance
ranking of an iterative random forest, on this exact matrix. Our result says
that the allocation of importance *between kinds of column* on this matrix
swings from 49% to 18% depending on which algorithm you ask.

This does **not** show their conclusion is wrong. In fact LightGBM puts more
weight on the quantum descriptors than on anything else (33%), which is
consistent with it. What it shows is that **an importance ranking on this matrix
cannot establish the claim on its own** — the same evidence, read through a
different algorithm, gives a different answer. The claim needs the kind of
support listed just above, and that is a precise, constructive criticism rather
than a dismissal.

### How LightGBM uses the features, and what happens if we give it more

A gradient-boosted tree model builds many small decision trees in sequence.
Each tree asks a handful of yes/no questions ("is downstream GC above 54%?"),
each question splitting the guides into two groups whose average scores differ
as much as possible. Each new tree is fitted to the *errors* of the ones before
it, so the model improves by correcting itself. Crucially, **each tree only uses
a few features**, and features that never produce a good split are never used.

That has two consequences the project tested:

**Giving it more features does almost nothing.** The model is capped at the top
300 features by default. Raising the cap:

| features allowed | ρ |
|---:|---:|
| 300 | 0.6068 |
| 600 | 0.6079 |
| 1,200 | 0.6082 |

**+0.001 for four times as many.** Beyond a few hundred, the extra columns are
correlated with ones already in and produce no new splits worth making.

**And it is not leaving information unused.** If a model were failing to exploit
a feature it holds, that feature would still correlate with the model's
remaining errors. Of the 6,480 features it is given, **none** has a correlation
above 0.10 with its own errors (the largest is 0.035). It has squeezed its
features dry. So the limit is the features, not the model — which is exactly why
model choice costs so little.

### The fixed-selection problem, which is a fair objection

In the sixteen-model comparison, the top-300 features were chosen **once, by
gradient boosting**, and then handed to every model. So each alternative was
being judged on features picked to suit a tree. A nearest-neighbour model or a
linear model might prefer different features entirely, and this design cannot
tell.

Two things to say about that:

1. **It is deliberate, because it isolates one variable.** With selection fixed,
   the comparison measures *the predictor*. Let each model select its own and
   you are comparing predictor-plus-selector pairs, which is a different and
   less interpretable experiment.
2. **We also ran the other experiment**, with each model selecting using its own
   importance scores — that is where the interpretability table above comes
   from. The ordering of the models does not change.

What neither version removes is that *the feature set itself* was engineered
while looking at tree-model results, so it is plausibly tree-friendly. The clean
answer to that is in Part 4: the neural network was given **raw DNA**, not our
features, and still did not win. The comparison that would settle it completely
— letting a neural network design its own features from scratch on these rows —
is what crisprHAL 2 is, and Part 2 is that comparison.

A note on the inherited work: the original 2023 paper selected its features with
an **iterative random forest**, not with gradient boosting, and read its
biological conclusions off that model's importance ranking. Our result that two
importance methods on this matrix can disagree at −0.02 applies directly to that
reasoning, which is one of the sharper things this project has to say.

---

## Part 9 — How much room is left, and can the limit be raised

### The apparent ceiling, and a correction to how we justified it

Two screens of the same guide library — WT-SpCas9 and eSpCas9 — agree with each
other at **ρ 0.810** on 33,567 shared guides. Standard reasoning (each
measurement is signal plus independent noise, so their correlation *is* the
reliability) puts the best achievable score against a single observed cut score
at about **√0.810 ≈ 0.90**.

> **We previously called 0.90 a *conservative* estimate on the grounds that the
> enzymes differ, so part of the disagreement is real biology rather than noise.
> That reasoning is probably backwards.** eSpCas9 is WT-SpCas9 with three point
> mutations, and the two screens use the **same library** — identical target
> sequences. So the 0.810 is largely two measurements of nearly the same
> quantity on exactly the same DNA, which is the *best case* for agreement. Two
> genuinely independent repeats, with a different library prep and a different
> enzyme, would likely agree **less** — making the true ceiling **lower** than
> 0.90, not higher. The number should be read as an **optimistic upper bound**,
> and it would be worth recomputing from true biological replicates if any are
> published.

Either way, GuideGauge at 0.707–0.718 and crisprHAL 2 at 0.697 are both some way
below it. The limit is not what either model is currently hitting.

### Why the learning curve saturates

Four times the training data buys **+0.008 ρ** (0.536 → 0.544), and the curve is
flat between 20,000 and 26,000 guides. That surprises people, so here is why it
should be expected.

More data helps a model in two ways: it lets the model estimate its parameters
more precisely, and it lets the model support a more complex hypothesis. Both
have already run out here.

- **Precision has run out.** The signal is largely additive (a straight line
  gets within 0.03 of the best model). Estimating a few hundred additive
  coefficients from 13,000 examples is already comfortable; another 13,000
  sharpens them by a negligible amount.
- **Complexity cannot be bought.** More rows would let a model fit a richer
  function only if the extra structure existed in the features — and the
  residual test says it does not: nothing the model holds still correlates with
  its own errors above 0.035.
- **And the label noise does not shrink.** Every new row carries the same
  measurement noise as the old ones. More rows average out noise *in the fitted
  parameters*, not noise in the test labels you are scored against. That part
  of the gap is fixed by the assay, not by the sample size.

So the saturation is not a surprise; it is what an additive model with
exhausted features and a noisy label is supposed to do.

### Can the ceiling be raised?

Yes, but not by anything in this repository — the ceiling is a property of the
*measurement*, so raising it means measuring differently.

1. **Average over replicates.** The ceiling applies to predicting a *single*
   noisy observation. Predicting the mean of three independent screens is a
   much easier target, and the apparent ceiling rises accordingly. This is the
   cheapest route and needs no new technique, only repetition.
2. **Measure cutting rather than survival.** The current label is depletion from
   a growing population, which mixes cutting with repair and growth. A direct
   readout — sequencing the broken ends, or measuring editing at defined sites
   — would remove whole categories of noise. More work, much better label.
3. **Report the residual properly.** Our own analysis says ~0.2 of Spearman is
   unexplained and that it is *not* in the features we have, *not* in the model
   class, *not* in the label noise, and *not* in the guide's own representation.
   What is left is features nobody has built — and after `i_shape` failed, the
   honest statement is that we do not currently know what they are.

### Where the remaining room is not

Worth stating, because each was a live hypothesis that got closed:

- not in **more rows** (+0.008),
- not in **model class** (0.116 across sixteen),
- not in **a neural network on raw sequence** (ties on the guide, loses on the
  flanks),
- not in **more features** (+0.001 for 4× the cap),
- not in **DNA shape** (+0.001 on top of flanks),
- not in **chromosome position** (collapses into the flanks),
- and not in **the label being dirty** (cleaning it helps by revealing the flank
  effect, not by raising the ceiling).

### Where everything here stops: human cells

The sharpest boundary we measured. The source paper also published a human
dataset sharing 6,216 columns with the bacterial one, so a model crosses with
**no change of representation at all** — a failure cannot be blamed on
mismatched features.

| | ρ |
|---|---:|
| *E. coli* model on *E. coli* | 0.531 |
| human model on human | 0.404 |
| ***E. coli* model on human** | **−0.048** |
| **human model on *E. coli*** | **−0.017** |

**Cross-kingdom transfer is useless in both directions** — very slightly worse
than guessing. And the mechanism inverts, which is why it is negative rather
than merely weak:

| | *E. coli* | human |
|---|---:|---:|
| GC content → cut score | **−0.201** | **+0.017** |
| melting temperature → cut score | −0.201 | +0.017 |

**This answers the question directly: the GC penalty does not generalise to
humans.** In *E. coli* it is the project's strongest single mechanism; in human
data it is inert. An *E. coli*-trained model applies "avoid GC-rich targets",
which in human cells is simply not a rule.

**Do histones explain this?** Partly, and the published literature supports the
mechanism, though we have not tested it ourselves:

- **Nucleosomes physically block Cas9.** Human DNA is wrapped around histone
  proteins in nucleosomes, and Cas9 cannot easily reach DNA inside one. Bacteria
  have no histones at all. So in human cells a large determinant of whether a
  guide works is *whether the target is accessible* — a variable that does not
  exist in our data and cannot be inferred from the sequence.
- **The published human effect of GC is non-monotonic, not absent.** A
  well-cited analysis found that very high *and* very low GC targets both work
  less well, with roughly 40–60% being the useful range, and that a *linear*
  association between GC and efficiency was not statistically significant. Our
  +0.017 is a linear correlation, so it is consistent with a real U-shape being
  invisible to the measure we used.
- **Accessibility dominates in a way it cannot in bacteria.** The same analysis
  found targets in promoter regions — which are kept open — cut better than
  targets in intergenic regions.

So the fair summary is: **in human cells, accessibility is a first-order
determinant and target GC is at best a weak non-linear one; in bacteria, there
are no nucleosomes and GC is a strong monotonic one.** Bacteria are not
"chromatin-free" — they pack DNA with proteins like HU and H-NS, which is what
our `e_nucleoid` set was about — but that packaging is not nucleosomal and does
not occlude targets the same way. Confirming the U-shape in the human data we
hold would be a cheap and worthwhile addition.

---

## Part 10 — What is new here, and what is imported

### Imported, and openly so

`b_energy` is the CRISPRoff authors' energy model; `c_folding` is ViennaRNA;
`f_methylation` is a text search; `h_offtarget` is a genome scan; `i_shape` uses
published dinucleotide scales; and the position sets read published GapR-seq,
Hi-C and RegulonDB/PRECISE-1K data with standard tools. Sources in Appendix A.
**Using them is not a contribution.**

What *is* a contribution is that **seven of them were measured against controls
and found unable to help**, with Part 5's rule explaining why in a way that
generalises.

### New

1. **The published dataset was decoded and made portable** — 99% of 6,232
   columns regenerate from any 20-letter guide, verified at 100% on 1.17 million
   values. Everything else here depends on it, and it is a reproducibility
   contribution in its own right.
2. **A rule for when a feature cannot possibly help, with a cheap test**, applied
   as a prediction rather than an excuse (Part 5) — and with its limit found.
3. **A control that overturned our own positive result** (Part 6), of a kind not
   standard in this field.
4. **Interpretability measured rather than asserted** — and then *diagnosed*.
   The disagreement between importance methods turns out to be a consequence of
   redundancy plus each algorithm's preference for binary or continuous columns,
   not of one model being wrong (Part 8). The consequence is that an importance
   ranking on this matrix cannot carry a biological claim on its own —
   including the inherited paper's.
5. **The flank effect characterised, not just reported** — a gradient rather than
   a motif, peaking at 250–500 letters, downstream-weighted, transferring across
   enzyme and organism at ~90%, and stopping dead at the kingdom boundary.
6. **Measured boundaries**, which is rarer than it should be: where the method
   stops (human cells, 0% transfer), what the data cannot answer (a different
   phylum, and why), and how much room is left.

---

## Part 11 — Where we corrected ourselves

Kept visible rather than quietly edited. Eleven corrections, of which three are
from the most recent round.

| claim | corrected to |
|---|---|
| "cleaning the label is worth +0.003, i.e. nothing" | true for the baseline, **false for the flank features**, which gain +0.164 on clean data against +0.080 |
| "the nearby letters are two thirds of the flank effect" | true on noisy data; on clean data the **distant averages are the larger half** |
| "we are 0.054 behind crisprHAL 2" | compared different guides; on the same guides it is parity |
| "parity, but against their published figure" | **their model re-run on our folds** scores 0.6971; we are ahead by +0.0107 in 5/5 folds |
| "the read-depth effect is genuine local chromosome state" | GC was tested against the **label** instead of against the **read-depth track**, where it correlates up to −0.70. It is library-prep bias plus mappability |
| "read depth is the largest unexplained effect left" | decomposes into replication timing (no rank gain), GC/mappability (already in the flanks), and ≈+0.007 R² that is actually new |
| "the flank effect transfers to another organism" | the cross-*enzyme* half holds; the only other-organism screen covers 4.3% of one chromosome |
| "LightGBM fits in 0.3 s, ~25,000× faster" | the 0.3 s was the **final fit only**, and was compared against their whole five-fold run. Controlled, like for like: **1.4 min against 105 min, about 74×**, at comparable memory |
| **"we are ahead by +0.0107"** | **tuning our own model alone is worth +0.0103, so the margin is within the tuning effect. The claim is parity** |
| **"0.90 is a conservative ceiling because the enzymes differ"** | **the two screens share their library and near-identical enzymes, so 0.90 is an optimistic upper bound, not a conservative one** |
| **"windowed composition beats a CNN four times over (+0.164 vs +0.019)"** | **those were different datasets. Like for like it is +0.080 vs +0.019 — still 4.2×, but the stated pair was wrong** |

### Bugs found, and the pattern in them

Seven in total; the five instructive ones:

1. An energy calculation used the wrong one of two similar quantities, so **zero**
   guides appeared to fall in the relevant range and the feature was a silent
   constant. After the fix, 63.7%.
2. Seventeen features per window turned out to be near-copies of plain GC
   content (correlation 0.998). Replaced by four deliberately different ones:
   **348 features instead of 408, four times the reach, and a better result.**
3. **An off-by-one in PAM detection scored a perfect 1.000 while being wrong** —
   testing positions +20 and +21 instead of +21 and +22 also passes on every
   guide. Caught only by demanding the detector return a value another module
   independently hard-codes. *A validation that a wrong answer can pass is not a
   validation* — and this was the second bug of exactly that shape.
4. The flank builder **silently returned zero features for every row** when a
   flag it depends on was unset. No error.
5. A resume check trusted a summary file rather than counting completed work, so
   a finished experiment's summary row was never written — which is how one
   permutation control appeared not to exist for two weeks when its data was on
   disk all along.

The pattern worth a line in the methods: **four of the seven returned confident,
plausible, wrong numbers rather than failing.** Those are the ones that need
cross-checks against an independently computed value, not better error handling.

---

## Part 12 — Next steps and how to frame the write-up

### Next

The order is settled: **harden the *E. coli* methodology, finish the paper,
then port the method to human.** Human editing is the medically important
target, and doing *E. coli* properly is what makes starting there cheap.

**Before the paper is written**

1. **Tune both models, or state parity.** The single highest-value item. Our
   tuning gain (+0.0103) is the size of the whole margin, so either crisprHAL 2
   gets the same search — about a day of CPU — or the paper claims parity and
   says why. *Parity is a perfectly good result; an unsupported lead is not.*
2. **Recompute the ceiling from true replicates** if any exist. The current
   0.90 rests on two screens that share a library, which makes it an optimistic
   bound rather than a conservative one.
3. **Add seeds to the single-seed arms.** The transfer, hybrid and human
   results are seed 41 only. The conclusions are large enough to survive, but a
   reviewer will ask.
4. **Confirm the human GC relationship is U-shaped** in the matrix we already
   hold. Cheap, and it upgrades "the mechanism does not transfer" to "the
   mechanism is replaced by a different one".

**Known-open, needing data nobody has**

5. **The ≈+0.007 R² abundance residual** — the only genuinely unexplained
   signal left. Needs local protein occupancy or similar.
6. **A genome-wide screen in a distant bacterium.** The one question analysis
   cannot close; worth saying plainly that it needs wet-lab work, because no
   one else has the dataset either.

**Then: the human phase**

7. **Port the method, not the model.** Cross-kingdom prediction measures ≈0 in
   both directions (−0.048 and −0.017), because the strand-invasion mechanism
   that dominates in bacteria is absent in human cells. What transfers is the
   *approach*: the redundancy test, testing a candidate cause against the
   mediator rather than the outcome, measuring interpretability instead of
   asserting it, and the controls that overturned two of our own results.
8. **Open with the obvious hypothesis.** This project's central finding is that
   DNA *outside* the 20-mer matters. There is no reason its compositional form
   should hold in human cells — but the human analogue of "context outside the
   target" is **chromatin accessibility and nucleosome positioning**, which is
   known to matter and which the published human matrix does not contain. That
   is a well-posed first question, and the whole apparatus for answering it
   already exists. The human feature matrix is already on disk at
   `data/raw/human_feature_matrix.csv`.

### A note on tooling, deliberately parked

Everything needed for a guide-design tool already exists — locate the target,
enumerate its NGG sites, compute base and flank features, rank. The inversion
result above says why that is the right shape for a tool and a global "optimal
guide" is not: guides are chosen from the few valid sites inside a target gene,
so the useful output is a ranking of real candidates with a reason attached.

**It is not a priority.** The findings are the contribution, the field already
has predictors, and a tool without the findings adds little. Noted here so the
option is not forgotten.

### How to frame it

**Not** "a better guide predictor". At parity, that framing invites exactly the
comparison this project loses — against a full-time lab with a GPU.

Frame it as: **what determines whether a CRISPR guide works in *E. coli*, and how
to tell in advance whether a proposed explanation can possibly help.** Then
every result is load-bearing:

- a published dataset decoded, verified, and made usable by anyone;
- a rule predicting which feature ideas cannot work, applied as a prediction and
  correct 7 times in 10, with its limit identified by the three misses;
- our own positive result overturned by the source experiment's own control;
- the first *measured* rather than asserted interpretability comparison in this
  literature, finding the inherited paper's own method internally inconsistent;
- the flank effect characterised as a gradient and shown to transfer across
  enzyme and organism at ~90% and across kingdoms at 0%;
- and, as a by-product, parity with the state of the art at roughly 1/74th of
  the CPU time.

### Honest caveats to carry into the paper

- **Parity, not a win.** The +0.0107 margin is inside the tuning effect.
- **Our speed advantage is CPU-to-CPU.** crisprHAL 2 as published is
  GPU-trained; we did not beat a GPU.
- **Cross-species is untested.** Both new screens are *Enterobacteriaceae*, and
  the *C. rodentium* one covers 4.3% of one chromosome. The competition's
  position is no better, but that does not make ours good.
- **Cross-enzyme evidence is narrow.** eSpCas9 differs from WT-SpCas9 by three
  point mutations on an identical guide library.
- **The ceiling is an optimistic bound**, for the same library-sharing reason.
- **Everything is bacterial.** Measured, not hedged: 0% transfer to human cells.
- **`h_offtarget` is a null result** and `f_methylation` is marginal, by their own
  permutation floors.
- **Loose cross-validation flatters positional features**, which is why they are
  reported under the strict scheme.

---

## Appendix A — Every feature set, grouped, with its data source

Nine sets, but only **four distinct ideas**. Grouping them this way is how they
should appear in the paper.

### Group 1 — the guide's own sequence (all redundant by Part 5's rule)

| set | what it measures | how it is computed | data source | Δρ |
|---|---|---|---|---:|
| `b_energy` | energy of guide–DNA binding, split into components | the CRISPRoff energy model, run on the guide | [CRISPRoff](https://github.com/RTH-tools/crisproff) repository, imported unchanged | +0.007 |
| `c_folding` | whether the guide RNA folds up on itself | ViennaRNA folding of spacer, and of spacer+scaffold together | ViennaRNA library | +0.007 |
| `d_mechanics` | how easily the duplex comes apart, position by position | nearest-neighbour thermodynamics (SantaLucia & Hicks 2004 parameters) | published parameter tables | +0.020 |
| `f_methylation` | Dam/Dcm methylation motifs at the PAM and seed | text search for `GATC` and `CCWGG` | reference genome | +0.003 |

### Group 2 — the surrounding DNA (the one that worked)

| set | what it measures | how it is computed | data source | Δρ |
|---|---|---|---|---:|
| `a_flank` | composition of the flanks at four scales | GC, purine fraction, longest letter run, and averages of the recovered quantum tables, over windows of 50/250/500/1000 letters each side; plus letter identity at the 10 nearest positions | reference genome + the quantum tables recovered in Part 3 | **+0.080** |
| `i_shape` | physical bendability and stiffness of the flanks | nine dinucleotide-step scales averaged over windows, plus a phased bend sum at the 10.5-letter helical repeat | published dinucleotide shape scales (two of which we had to correct — one scale was corrupt, two were not strand-symmetric as published) | +0.036 alone, **+0.001** on top of `a_flank` |

### Group 3 — where the guide sits on the chromosome (all one variable, all collapse into Group 2)

| set | what it measures | how it is computed | data source | Δρ |
|---|---|---|---|---:|
| `d_supercoiling` | DNA twisting, nominally | read density of GapR ChIP tracks in windows, plus ratios against the controls | GEO **GSE152880** (GapR-seq), including its untagged and rifampicin control arms | +0.045, of which +0.004 is twisting-specific |
| `e_nucleoid` | 3D crowding and dependence on packaging proteins | contact counts from a Hi-C matrix; nucleoid-protein dependence | Lioy et al. 2018 Hi-C matrices | +0.016 |
| `g_transcription` | distance to a promoter, strand, expression | promoters re-located by matching the 80-letter sequence each entry ships with (avoiding a coordinate-system mismatch), joined to expression | RegulonDB + PRECISE-1K | +0.015 |
| `z_mapgc` | the sequencing artefact behind all of the above | 25-letter uniqueness fraction and GC fraction, at four window sizes — **the 8 columns of Part 6** | reference genome only, no experiment | +0.038 alone, **−0.000** on top of `a_flank` |

### Group 4 — the rest of the genome

| set | what it measures | how it is computed | data source | Δρ |
|---|---|---|---|---:|
| `h_offtarget` | copy number and burden of near-matching sites | genome-wide scan for sequences within a few letters of the target | reference genome | +0.000 — **a null result by its own floor** |

Label sources: the published `cut.score` from Guo et al.'s 2018 *E. coli*
depletion screen (13,880 guides, via Noshay et al.'s matrix), and the
read-count-filtered re-derivation shipped with crisprHAL (33,567 guides). Genomes:
*E. coli* NC_000913.2, and *C. rodentium* ICC168 NC_013716.1 (5,346,659 letters,
54.7% GC) — the only new download this round, made because the transfer question
could not be asked without it.

---

## Appendix B — Papers referred to

The two that would not load for us earlier, in case they open for you:

- crisprHAL 2 — *Better data for better predictions: data curation improves deep
  learning for sgRNA/Cas9 prediction*, PeerJ 2026.
  [PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC12903899/) ·
  [PeerJ, which did load](https://peerj.com/articles/20706/)
- DeepCC9 — *An interpretable deep learning framework uncovers features governing
  CRISPR-Cas9 genome-editing efficiency*, Bioinformatics 2026.
  [PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC13384063/) ·
  [Oxford Academic, which did load](https://academic.oup.com/bioinformatics/article/42/7/btag483/8723703)

The rest:

- [Noshay et al. — *Quantum biological insights into CRISPR-Cas9 sgRNA efficiency*, NAR 2023](https://academic.oup.com/nar/article/51/19/10147/7279034) — the matrix this project decodes
- [crisprHAL 1 — *A generalizable Cas9/sgRNA prediction model*, Nat Commun 2023](https://www.nature.com/articles/s41467-023-41143-7) — the cross-species claims belong here, not to crisprHAL 2
- [Liu et al. — *Sequence features associated with the cleavage efficiency of CRISPR/Cas9*, Sci Rep 2016](https://www.nature.com/articles/srep19675) — the human GC U-shape and the chromatin-accessibility effect
- [*Nucleosomes impede Cas9 access to DNA*, eLife 2016](https://elifesciences.org/articles/12677) and [*Nucleosomes inhibit Cas9 cleavage in vivo*, PNAS 2018](https://www.pnas.org/content/115/38/9351) — why human and bacterial determinants differ
- [*Improved prediction of bacterial CRISPRi guide efficiency*, Genome Biology 2023](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03153-y) — the closest adjacent problem, and a CRISPRi (not cutting) screen
- [FDA — approval of the first CRISPR therapy](https://www.fda.gov/news-events/press-announcements/fda-approves-first-gene-therapies-treat-patients-sickle-cell-disease)
