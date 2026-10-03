"""Family D1 -- how hard the target duplex is to open, position by position.

Cas9 invades double-stranded DNA. Before an R-loop can form the duplex has to
breathe open next to the PAM, so the *local* stability of the target -- not
its average GC content -- is the physically relevant quantity. Guo et al.
(2018) already reported that the melting temperature of the guide:target
duplex was their single strongest predictor in E. coli, but the matrix carries
only one global Tm column.

This module computes the unified nearest-neighbour thermodynamics of
SantaLucia & Hicks (2004) for DNA:DNA, then reads them out the way the
mechanism suggests: as a sliding profile, as the PAM-proximal versus
PAM-distal split, and over the flanks that have to be pushed aside.

The three-window Tm split (positions 1-5, 6-13, 14-20) is the one used by the
Doench/Azimuth rule set, kept here so E. coli and human models can be compared
on the same footing.

Parameters: SantaLucia & Hicks, Annu Rev Biophys Biomol Struct 33:415 (2004).
Units: dH kcal/mol, dS cal/(mol K), dG37 kcal/mol.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .. import config
from ..io_utils import save_block

PREFIX = "eng.mech"

# Unified nearest-neighbour parameters, 1 M NaCl.
NN_DH = {
    "AA": -7.9, "TT": -7.9, "AT": -7.2, "TA": -7.2,
    "CA": -8.5, "TG": -8.5, "GT": -8.4, "AC": -8.4,
    "CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2,
    "CG": -10.6, "GC": -9.8, "GG": -8.0, "CC": -8.0,
}
NN_DS = {
    "AA": -22.2, "TT": -22.2, "AT": -20.4, "TA": -21.3,
    "CA": -22.7, "TG": -22.7, "GT": -22.4, "AC": -22.4,
    "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2,
    "CG": -27.2, "GC": -24.4, "GG": -19.9, "CC": -19.9,
}
NN_DG37 = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
}
INIT_DH = {"GC": 0.1, "AT": 2.3}
INIT_DS = {"GC": -2.8, "AT": 4.1}
INIT_DG37 = {"GC": 0.98, "AT": 1.03}

R = 1.987  # cal / (mol K)
DEFAULT_STRAND_CONC = 0.5e-6  # mol/L; only shifts Tm by a constant offset

# The Doench/Azimuth window split, in 1-based protospacer coordinates.
TM_WINDOWS = {"w1_5": (1, 5), "w6_13": (6, 13), "w14_20": (14, 20)}

SLIDING_WINDOW = 5


def _init_key(base: str) -> str:
    return "GC" if base in "GC" else "AT"


def duplex_thermo(seq: str) -> tuple[float, float, float]:
    """(dH, dS, dG37) for a perfectly matched DNA duplex."""
    seq = seq.upper()
    if len(seq) < 2 or any(c not in "ACGT" for c in seq):
        return (np.nan, np.nan, np.nan)
    dh = INIT_DH[_init_key(seq[0])] + INIT_DH[_init_key(seq[-1])]
    ds = INIT_DS[_init_key(seq[0])] + INIT_DS[_init_key(seq[-1])]
    dg = INIT_DG37[_init_key(seq[0])] + INIT_DG37[_init_key(seq[-1])]
    for i in range(len(seq) - 1):
        step = seq[i : i + 2]
        dh += NN_DH[step]
        ds += NN_DS[step]
        dg += NN_DG37[step]
    return dh, ds, dg


def melting_temperature(seq: str, conc: float = DEFAULT_STRAND_CONC) -> float:
    """Tm in degrees C for a non-self-complementary duplex."""
    dh, ds, _ = duplex_thermo(seq)
    if np.isnan(dh):
        return np.nan
    tm = (dh * 1000.0) / (ds + R * np.log(conc / 4.0)) - 273.15
    return float(tm)


def stacking_profile(seq: str) -> np.ndarray:
    """dG37 of each dinucleotide step, 5' -> 3'."""
    seq = seq.upper()
    return np.array(
        [NN_DG37.get(seq[i : i + 2], np.nan) for i in range(len(seq) - 1)]
    )


def sliding_dg(seq: str, window: int = SLIDING_WINDOW) -> np.ndarray:
    """dG37 of every `window`-nt sub-duplex -- a local breathing profile."""
    out = []
    for i in range(len(seq) - window + 1):
        out.append(duplex_thermo(seq[i : i + window])[2])
    return np.asarray(out, dtype=float)


def features_for_guide(proto: str, pam: str, upstream: str, downstream: str) -> dict:
    rec: dict[str, float] = {}

    dh, ds, dg = duplex_thermo(proto)
    rec[f"{PREFIX}.proto.dH"] = dh
    rec[f"{PREFIX}.proto.dS"] = ds
    rec[f"{PREFIX}.proto.dG37"] = dg
    rec[f"{PREFIX}.proto.tm"] = melting_temperature(proto)

    # Doench-style window split.
    for name, (lo, hi) in TM_WINDOWS.items():
        sub = proto[lo - 1 : hi]
        rec[f"{PREFIX}.tm.{name}"] = melting_temperature(sub)
        rec[f"{PREFIX}.dG.{name}"] = duplex_thermo(sub)[2]

    # PAM-proximal (seed) versus PAM-distal halves.
    seed, distal = proto[10:], proto[:10]
    rec[f"{PREFIX}.dG.seed"] = duplex_thermo(seed)[2]
    rec[f"{PREFIX}.dG.distal"] = duplex_thermo(distal)[2]
    rec[f"{PREFIX}.dG.seed_minus_distal"] = (
        rec[f"{PREFIX}.dG.seed"] - rec[f"{PREFIX}.dG.distal"]
    )

    # Local breathing profile across protospacer + PAM: where is it weakest,
    # and how close is that weak point to the PAM, where invasion starts?
    span = proto + (pam or "")
    prof = sliding_dg(span)
    if prof.size and not np.all(np.isnan(prof)):
        rec[f"{PREFIX}.slide.min"] = float(np.nanmin(prof))
        rec[f"{PREFIX}.slide.max"] = float(np.nanmax(prof))
        rec[f"{PREFIX}.slide.mean"] = float(np.nanmean(prof))
        rec[f"{PREFIX}.slide.sd"] = float(np.nanstd(prof))
        # dG37 is negative for a stable duplex, so the *largest* value marks
        # the easiest place to melt. Distance is counted back from the PAM,
        # which is where strand invasion starts.
        rec[f"{PREFIX}.slide.weakest_from_pam"] = float(
            len(prof) - 1 - int(np.nanargmax(prof))
        )
        rec[f"{PREFIX}.slide.strongest_from_pam"] = float(
            len(prof) - 1 - int(np.nanargmin(prof))
        )

    stacks = stacking_profile(span)
    if stacks.size and not np.all(np.isnan(stacks)):
        rec[f"{PREFIX}.stack.mean"] = float(np.nanmean(stacks))
        rec[f"{PREFIX}.stack.min"] = float(np.nanmin(stacks))
        rec[f"{PREFIX}.stack.sd"] = float(np.nanstd(stacks))

    # Flanks the enzyme has to push aside.
    for label, seg in (("up25", upstream[-25:]), ("dn25", downstream[:25]),
                       ("up100", upstream[-100:]), ("dn100", downstream[:100])):
        if seg and all(c in "ACGT" for c in seg):
            rec[f"{PREFIX}.flank.{label}.dG37"] = duplex_thermo(seg)[2]
            rec[f"{PREFIX}.flank.{label}.dG_per_bp"] = (
                duplex_thermo(seg)[2] / len(seg)
            )
        else:
            rec[f"{PREFIX}.flank.{label}.dG37"] = np.nan
            rec[f"{PREFIX}.flank.{label}.dG_per_bp"] = np.nan
    return rec


def build(index: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in index.iterrows():
        rec: dict[str, float] = {config.ID_COL: row[config.ID_COL]}
        proto = row.get("protospacer")
        if isinstance(proto, str) and "N" not in proto:
            rec.update(
                features_for_guide(
                    proto,
                    row.get("pam") if isinstance(row.get("pam"), str) else "",
                    row.get("upstream") if isinstance(row.get("upstream"), str) else "",
                    row.get("downstream") if isinstance(row.get("downstream"), str) else "",
                )
            )
        rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    from ..io_utils import load_guide_index

    idx = load_guide_index()
    block = build(idx)
    path = save_block(block, "d_mechanics")
    print(f"d_mechanics: {block.shape[0]:,} rows x {block.shape[1] - 1:,} features -> {path}")


if __name__ == "__main__":
    main()
