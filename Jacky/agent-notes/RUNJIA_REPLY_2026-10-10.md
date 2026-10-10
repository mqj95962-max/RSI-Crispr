# Reply to Runjia — 10 Oct 2026

His commit: `e4a71fc` (*accept Jacky's intergenic orientation result; pick between
the two ceiling routes*). Orientation is in Part 5; ceiling framing updated.

---

## 1. Process — `Runjia/sgrna/` edits

**Agree.** From now on, any change under `Runjia/sgrna/` will be called out in:

- the PR body (bullet list of files + one-line intent), **and**
- a line in `Jacky/agent-notes/STATUS.md` under “Pipeline edits touching Runjia/sgrna”.

Additive / backwards-compatible still fine; the issue was silent presence on
`main` while the follow-up said “not pushed.” That wording was stale after
`7b8266e` / the follow-up merge — sorry.

Convention locked in `STATUS.md` (this commit).

---

## 2. Ceiling — accept 0.937 as the working number

**Do not push back.** His ratio argument is right, and it uses *our* provenance.

Short version of why:

- Fig 2b + Spearman–Brown → 0.968 only if the *activity score* is about as
  reliable as a single arm’s read counts.
- Our label is `|Z|` of a **Cas9 / dCas9 ratio**. A contrast is never more
  reliable than its arms; shared signal between arms (exactly what the dCas9
  control is for — library abundance, and much of the gene-level structure)
  makes the ratio *less* reliable.
- So 0.968 is an optimistic bound that assumes ~0% shared signal. His shared-
  fraction table (→ 0.956 / 0.940 / 0.905) shows how fast that falls.
- Fig 2c needs none of that: both sides are already the published score →
  reliability = 0.878 → ceiling **0.937 [0.929, 0.945]**, and library-design
  differences bias it *down* (conservative).

We already agreed on the substance: 0.90–0.93 was too low; “replicates nobody
published” was wrong. The only dispute was which statistic to lead with.
**Lead with ≈0.94 (tiling); keep 0.94–0.97 as the range.** SLICER ~77% of the
way. That matches his Part 10 update.

Optional one-liner for him if useful: “Agreed — Fig 2c is the right lead;
Fig 2b+SB stays as an optimistic upper edge of the range, not the working
figure.”

---

## Suggested short reply (paste-ready)

> Thanks — glad the intergenic orientation landed in Part 5.
>
> On process: agreed. Anything that touches `Runjia/sgrna/` will be listed in
> the PR body and in `STATUS.md` so your sync can’t silently overwrite it.
> Sorry the follow-up still said “not pushed” after those three files were
> already on main — and thanks for catching / adopting the `guide_overlap` fix.
>
> On the ceiling: I won’t push back. The Cas9/dCas9 ratio argument is right
> given the provenance we traced; Fig 2b+SB to 0.968 only holds if the control
> shares nothing, which it doesn’t. Happy with ≈0.94 working and 0.94–0.97 as
> the range, tiling as the lead. Substance agreement stands.

---

## Still open after his acceptance (unchanged priorities)

| item | status |
|---|---|
| Bootstrap CI on orientation Δ (opposite − same) | nice next if we deepen 2B.3 |
| Better orientation label (TSS/operon, not nearest gene) | medium |
| 2A.1 dCas9 abundance / SRA | optional, blocked on counts |
| Distant bacterium screen | wet-lab / his Part 12 |
