"""Family C -- sgRNA folding, done on the transcript that is actually made.

The published matrix has one structure column, `sgRNA.structuresgRNA.raw`, and
earlier experiments in this project tested structural *proxies*. Neither is
the quantity the literature identifies.

Moreb & Lynch (2022, The CRISPR Journal) meta-analysed 39 guide libraries,
including a 1.2-million-guide self-targeting E. coli set, and found two
separate, sigmoidal failure modes:

  * the spacer folding back on itself -- activity collapses once the spacer's
    own MFE drops below about -5 kcal/mol;
  * the spacer base-pairing with the constant scaffold and disrupting the
    repeat:anti-repeat duplex -- activity falls once that duplex stability
    drops below about -15 kcal/mol.

Both are threshold effects on a minority of guides, which is exactly the shape
a boosted tree captures and a linear term does not. The single most common way
this feature family gets written off as useless is folding the bare 20-mer;
the sgRNA that E. coli transcribes is spacer + scaffold, so that is what gets
folded here.

Columns
-------
eng.fold.spacer_mfe          MFE of the 20 nt spacer alone
eng.fold.full_mfe            MFE of spacer + SpCas9 scaffold
eng.fold.scaffold_delta      full MFE minus the scaffold's own MFE
eng.fold.duplex_stability    best spacer:scaffold hybrid energy (RNAduplex)
eng.fold.seed_accessibility  mean unpaired probability over the seed
eng.fold.spacer_paired_frac  fraction of spacer bases paired in the full fold
eng.fold.below_*             the two published thresholds, as indicators
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

from .. import config
from ..io_utils import save_block

PREFIX = "eng.fold"

SPACER_MFE_THRESHOLD = -5.0      # Moreb & Lynch 2022
DUPLEX_THRESHOLD = -15.0         # Moreb & Lynch 2022
SEED_LENGTH = 10                 # PAM-proximal nucleotides, i.e. spacer 11..20


def _rna() -> "module":
    try:
        import RNA
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "ViennaRNA's Python bindings are needed: pip install viennarna"
        ) from exc
    return RNA


@lru_cache(maxsize=1)
def _scaffold_mfe() -> float:
    RNA = _rna()
    return float(RNA.fold(config.SPCAS9_SCAFFOLD)[1])


def _to_rna(dna: str) -> str:
    return dna.upper().replace("T", "U")


@lru_cache(maxsize=200_000)
def _fold(seq: str) -> tuple[str, float]:
    RNA = _rna()
    structure, energy = RNA.fold(seq)
    return structure, float(energy)


@lru_cache(maxsize=200_000)
def _duplex(a: str, b: str) -> float:
    RNA = _rna()
    return float(RNA.duplexfold(a, b).energy)


@lru_cache(maxsize=200_000)
def _unpaired_probabilities(seq: str) -> tuple[float, ...]:
    """Per-base probability of being unpaired, from the partition function.

    `bpp()` hands back a 1-indexed upper-triangular matrix as nested lists.
    Summing it with a Python double loop costs ~9,000 iterations per guide and
    dominates the whole family, so it goes through numpy in one step instead.
    """
    RNA = _rna()
    fc = RNA.fold_compound(seq)
    fc.pf()
    bpp = np.asarray(fc.bpp(), dtype=float)  # (n+1, n+1), row/col 0 unused
    paired = bpp.sum(axis=0) + bpp.sum(axis=1)
    return tuple(np.clip(1.0 - paired[1 : len(seq) + 1], 0.0, 1.0))


def features_for_spacer(spacer_dna: str) -> dict[str, float]:
    spacer = _to_rna(spacer_dna)
    scaffold = config.SPCAS9_SCAFFOLD
    full = spacer + scaffold

    _, spacer_mfe = _fold(spacer)
    full_struct, full_mfe = _fold(full)
    duplex = _duplex(spacer, scaffold)

    unpaired = np.asarray(_unpaired_probabilities(full))
    seed_slice = unpaired[len(spacer) - SEED_LENGTH : len(spacer)]
    spacer_struct = full_struct[: len(spacer)]
    paired_frac = sum(ch != "." for ch in spacer_struct) / len(spacer)

    return {
        f"{PREFIX}.spacer_mfe": spacer_mfe,
        f"{PREFIX}.full_mfe": full_mfe,
        f"{PREFIX}.scaffold_delta": full_mfe - _scaffold_mfe(),
        f"{PREFIX}.duplex_stability": duplex,
        f"{PREFIX}.seed_accessibility": float(np.mean(seed_slice)),
        f"{PREFIX}.spacer_accessibility": float(np.mean(unpaired[: len(spacer)])),
        f"{PREFIX}.spacer_paired_frac": float(paired_frac),
        f"{PREFIX}.below_spacer_threshold": int(spacer_mfe < SPACER_MFE_THRESHOLD),
        f"{PREFIX}.below_duplex_threshold": int(duplex < DUPLEX_THRESHOLD),
    }


def build(index: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in index.iterrows():
        rec: dict[str, float] = {config.ID_COL: row[config.ID_COL]}
        proto = row.get("protospacer")
        if isinstance(proto, str) and "N" not in proto:
            rec.update(features_for_spacer(proto))
        rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    from ..io_utils import load_guide_index

    idx = load_guide_index()
    block = build(idx)
    path = save_block(block, "c_folding")
    print(f"c_folding: {block.shape[0]:,} rows x {block.shape[1] - 1:,} features -> {path}")
    for col, thresh in (
        (f"{PREFIX}.below_spacer_threshold", SPACER_MFE_THRESHOLD),
        (f"{PREFIX}.below_duplex_threshold", DUPLEX_THRESHOLD),
    ):
        n = int(block[col].sum())
        print(f"  {n:,} guides fall past the published threshold {thresh} kcal/mol ({col})")


if __name__ == "__main__":
    main()
