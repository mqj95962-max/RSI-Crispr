"""Family F -- Dam and Dcm methylation at the target site.

E. coli methylates its own chromosome at two motifs, and unlike CpG
methylation in mammals it does so essentially everywhere:

    GATC    Dam methyltransferase, N6-methyladenine
    CCWGG   Dcm methyltransferase, 5-methylcytosine  (W = A or T)

Until recently SpCas9 was treated as methylation-blind. Two 2026 results
sharpen that. A machine-learning study of SaCas9 in bacteria found roughly
ten-fold lower activity when a Dam GATC overlaps the PAM-proximal positions,
confirmed by comparing wild-type with Dam-deficient E. coli; and structural
work on ThermoCas9 showed that methylating a PAM cytosine blocks *binding*
outright rather than catalysis. The sensitive zone in both cases is the PAM
and its immediate flank, not the body of the protospacer.

Whether SpCas9 shows the same effect in this screen is a one-regex question
with a clean mechanistic prediction attached, which makes it unusually good
value: the whole family costs almost nothing to compute and the expected
result is a small effect confined to an identifiable subset of guides -- a
readable SHAP story even if the global R-squared barely moves.

Columns
-------
eng.meth.dam.*   GATC overlap with protospacer / seed / PAM / flank, counts,
                 and distance from the cut site to the nearest site
eng.meth.dcm.*   the same for CCWGG
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .. import config
from ..io_utils import save_block

PREFIX = "eng.meth"

DAM = re.compile(r"(?=(GATC))")
DCM = re.compile(r"(?=(CC[AT]GG))")

SEED_START = 11   # 1-based protospacer position where the PAM-proximal seed begins
FLANK_NT = 20     # how far either side of the PAM to look for a nearby site


def _site_starts(seq: str, pattern: re.Pattern) -> list[int]:
    """0-based start offsets of every (overlapping) motif occurrence."""
    return [m.start() for m in pattern.finditer(seq)]


def _motif_features(name: str, pattern: re.Pattern, motif_len: int,
                    proto: str, pam: str, upstream: str, downstream: str) -> dict:
    rec: dict[str, float] = {}

    # A window that covers everything the enzyme touches, with coordinates
    # measured relative to protospacer position 1.
    up = upstream[-FLANK_NT:] if upstream else ""
    dn = downstream[:FLANK_NT] if downstream else ""
    window = up + proto + pam + dn
    proto_start = len(up)
    pam_start = proto_start + len(proto)

    starts = _site_starts(window, pattern)

    def overlaps(lo: int, hi: int) -> int:
        """Does any motif occurrence intersect [lo, hi) in window coords?"""
        return int(any(s < hi and s + motif_len > lo for s in starts))

    rec[f"{PREFIX}.{name}.in_protospacer"] = overlaps(proto_start, pam_start)
    rec[f"{PREFIX}.{name}.in_seed"] = overlaps(
        proto_start + SEED_START - 1, pam_start
    )
    rec[f"{PREFIX}.{name}.in_distal"] = overlaps(
        proto_start, proto_start + SEED_START - 1
    )
    rec[f"{PREFIX}.{name}.in_pam"] = overlaps(pam_start, pam_start + len(pam))
    rec[f"{PREFIX}.{name}.in_pam_flank"] = overlaps(
        pam_start, pam_start + len(pam) + 5
    )
    rec[f"{PREFIX}.{name}.n_in_window"] = float(len(starts))
    rec[f"{PREFIX}.{name}.n_in_protospacer"] = float(
        sum(proto_start <= s < pam_start for s in starts)
    )

    # Cas9 cuts between protospacer positions 17 and 18.
    cut = proto_start + 17
    if starts:
        rec[f"{PREFIX}.{name}.nearest_to_cut"] = float(
            min(abs((s + motif_len / 2) - cut) for s in starts)
        )
    else:
        rec[f"{PREFIX}.{name}.nearest_to_cut"] = float(len(window))
    return rec


def build(index: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in index.iterrows():
        rec: dict[str, float] = {config.ID_COL: row[config.ID_COL]}
        proto = row.get("protospacer")
        if isinstance(proto, str) and "N" not in proto:
            pam = row.get("pam") if isinstance(row.get("pam"), str) else ""
            up = row.get("upstream") if isinstance(row.get("upstream"), str) else ""
            dn = row.get("downstream") if isinstance(row.get("downstream"), str) else ""
            rec.update(_motif_features("dam", DAM, 4, proto, pam, up, dn))
            rec.update(_motif_features("dcm", DCM, 5, proto, pam, up, dn))
            rec[f"{PREFIX}.any_in_seed"] = float(
                max(rec[f"{PREFIX}.dam.in_seed"], rec[f"{PREFIX}.dcm.in_seed"])
            )
        rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    from ..io_utils import load_guide_index

    idx = load_guide_index()
    block = build(idx)
    path = save_block(block, "f_methylation")
    print(f"f_methylation: {block.shape[0]:,} rows x {block.shape[1] - 1:,} features -> {path}")
    for col in (f"{PREFIX}.dam.in_seed", f"{PREFIX}.dam.in_pam_flank",
                f"{PREFIX}.dcm.in_seed"):
        if col in block:
            n = int(np.nansum(block[col]))
            print(f"  {n:,} guides ({n / len(block):.1%}) have {col}")


if __name__ == "__main__":
    main()
