"""Family I -- the physical shape of the DNA around the target, not its letters.

Family A established that the DNA *outside* the 20-mer carries the one large
effect in this dataset, and that it behaves as a smooth compositional gradient
rather than a motif. `a_flank` describes that DNA by its letter composition.
This module asks the next question: is the effect better described by what the
flanking DNA *physically does* -- how readily it bends, how floppy it is, how
much force a given deformation costs?

The mechanism this is reaching for: Cas9 does not act on an isolated 20-mer.
It binds, bends and locally unwinds a stretch of duplex, and the DNA on either
side has to accommodate that distortion. Two flanks with identical GC content
can differ substantially in how cheaply they deform, because deformability is
a property of neighbouring base *pairs* (dinucleotide steps), not of base
counts.

Only the flanks are described here, and that is a deliberate restriction
rather than an omission. `diagnose.redundancy` showed the published matrix to
be a lossless encoding of the 20 nt protospacer, so any quantity computed from
the protospacer alone -- however good its mechanism -- is a deterministic
function of columns already present and adds exactly zero information. Shape
descriptors of the protospacer would be in that position. Shape descriptors of
the flanks are not.

Scales come from `data/reference/dinucleotide_shape_scales.tsv`; see that file
for provenance and for the two anomalies it carries. Two decisions are made
here rather than in the data, so they stay visible:

  * `stiff_roll` is excluded. Four of its sixteen steps are 0.02 where the
    rest are 16-26, which is a data artifact, not physics.
  * `bend_major` and `bend_minor` are symmetrised over reverse complements.
    As published they are not RC-symmetric, yet a duplex step property must
    be: the same physical step read from the other strand cannot bend
    differently. Without this, every feature would depend on which strand the
    guide happens to target.

Features, per flank side and window:

  mean      the average shape value over the window -- the gradient-style
            summary that worked for family A
  sd        how variable the shape is across the window -- a uniform stretch
            and an alternating stiff/floppy one can share a mean
  curv      a phased bend sum: each step's net bend toward the major groove is
            treated as a vector rotating with the ~10.5 bp helical repeat and
            summed, so bends that reinforce each other every turn add up and
            bends that cancel do not. This is the standard way intrinsic
            curvature is estimated from a bend scale, but with a bending
            *propensity* scale in place of measured wedge angles, so read it
            as a proxy for coherent curvature, not as degrees of bend.

Also computed: the 10 bp immediately outside the protospacer on each side --
the DNA that has to distort first -- as its own window.

Caveat for the write-up: these scales are themselves sequence-derived, so they
are functions of the flank sequence. They are not new information about the
flanks; they are a *reparameterisation* of it, and the claim being tested is
that a 2 nt-step physical scale is a better-shaped variable than letter
composition. `d_mechanics` beating raw GC (+0.020 against the same
information) is the precedent for that being possible, and the redundancy
diagnostic is what decides it.
"""

from __future__ import annotations

import csv
import pathlib

import numpy as np
import pandas as pd

from .. import config
from ..io_utils import save_block

PREFIX = "eng.shape"

SCALE_PATH = (
    pathlib.Path(__file__).resolve().parents[3]
    / "data" / "reference" / "dinucleotide_shape_scales.tsv"
)

# See the module docstring: corrupt in the source table.
EXCLUDED_SCALES = ("stiff_roll",)

# Not RC-symmetric as published, so averaged with their reverse complements.
SYMMETRISE = ("bend_major", "bend_minor")

# Flank windows, in nt outward from the protospacer/PAM edge. 10 nt is the DNA
# that distorts first; 250 nt is where family A's gradient was still measurable.
WINDOWS = (10, 25, 50, 250)

HELICAL_REPEAT = 10.5  # bp per turn of B-DNA

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def load_scales(path: pathlib.Path = SCALE_PATH) -> dict[str, dict[str, float]]:
    """Dinucleotide scales, with the documented exclusion and symmetrisation."""
    with open(path) as fh:
        rows = [line for line in fh if not line.startswith("#")]
    scales: dict[str, dict[str, float]] = {}
    for rec in csv.DictReader(rows, delimiter="\t"):
        name = rec.pop("scale")
        if name in EXCLUDED_SCALES:
            continue
        table = {step: float(val) for step, val in rec.items()}
        if name in SYMMETRISE:
            table = {
                step: (val + table[_revcomp(step)]) / 2.0
                for step, val in table.items()
            }
        scales[name] = table
    return scales


SCALES = load_scales()
SCALE_NAMES = tuple(sorted(SCALES))


def profile(seq: str, scale: dict[str, float]) -> np.ndarray:
    """The scale's value at every dinucleotide step along `seq`, 5' -> 3'."""
    seq = seq.upper()
    return np.array(
        [scale.get(seq[i : i + 2], np.nan) for i in range(len(seq) - 1)],
        dtype=float,
    )


def phased_bend(seq: str, repeat: float = HELICAL_REPEAT) -> float:
    """Coherent curvature proxy: |sum of net bend as a phased vector| per step.

    Each step contributes its major-minor groove bend preference, rotated by
    the helical phase at that position. Bends spaced one helical turn apart
    point the same way and reinforce; bends half a turn apart cancel. Divided
    by the number of steps so windows of different length stay comparable.
    """
    major = profile(seq, SCALES["bend_major"])
    minor = profile(seq, SCALES["bend_minor"])
    net = major - minor
    if net.size == 0 or np.all(np.isnan(net)):
        return float("nan")
    phase = 2.0 * np.pi * np.arange(net.size) / repeat
    x = np.nansum(net * np.cos(phase))
    y = np.nansum(net * np.sin(phase))
    return float(np.hypot(x, y) / net.size)


def features_for_guide(upstream: str, downstream: str) -> dict:
    """Shape summaries of the flanking DNA only -- never the protospacer."""
    rec: dict[str, float] = {}
    for side, flank in (("up", upstream), ("dn", downstream)):
        for width in WINDOWS:
            # Upstream runs toward the protospacer, so its last `width` nt are
            # the ones adjacent to it; downstream's first `width` nt are.
            seg = flank[-width:] if side == "up" else flank[:width]
            tag = f"{side}{width}"
            usable = len(seg) >= 2 and all(c in "ACGT" for c in seg.upper())
            for name in SCALE_NAMES:
                key = f"{PREFIX}.{name}.{tag}"
                if not usable:
                    rec[f"{key}.mean"] = np.nan
                    rec[f"{key}.sd"] = np.nan
                    continue
                prof = profile(seg, SCALES[name])
                if np.all(np.isnan(prof)):
                    rec[f"{key}.mean"] = np.nan
                    rec[f"{key}.sd"] = np.nan
                else:
                    rec[f"{key}.mean"] = float(np.nanmean(prof))
                    rec[f"{key}.sd"] = float(np.nanstd(prof))
            rec[f"{PREFIX}.curv.{tag}"] = (
                phased_bend(seg) if usable else np.nan
            )
    return rec


def build(index: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in index.iterrows():
        rec: dict[str, float] = {config.ID_COL: row[config.ID_COL]}
        up = row.get("upstream")
        dn = row.get("downstream")
        if isinstance(up, str) or isinstance(dn, str):
            rec.update(
                features_for_guide(
                    up if isinstance(up, str) else "",
                    dn if isinstance(dn, str) else "",
                )
            )
        rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    from ..io_utils import load_guide_index

    idx = load_guide_index()
    block = build(idx)
    path = save_block(block, "i_shape")
    print(f"i_shape: {block.shape[0]:,} rows x {block.shape[1] - 1:,} features -> {path}")


if __name__ == "__main__":
    main()
