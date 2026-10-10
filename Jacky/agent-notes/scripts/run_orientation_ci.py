#!/usr/bin/env python3
"""Grouped bootstrap CI on intergenic orientation Δ asymmetry.

Reuses cached intergenic X/y from run_intergenic_orientation.py (builds if
missing). For each orientation stratum, gets OOF predictions for
upstream_only and downstream_matched arms, then block-bootstraps 100 kb
groups to CI:

  asymmetry = ρ(y, pred_dn_matched) - ρ(y, pred_up)
  delta     = asymmetry(opposite_gene) - asymmetry(same_as_gene)

Primary set: s6_good. Sensitivity: s4_hq.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(os.environ.get("RSI09_ROOT", "/workspace/Runjia")).resolve()
os.environ["RSI09_ROOT"] = str(ROOT)
sys.path.insert(0, str(ROOT))

# Ensure project markers exist before importing config-heavy modules
(ROOT / "external_data").mkdir(parents=True, exist_ok=True)
(ROOT / "data" / "raw").mkdir(parents=True, exist_ok=True)

from sgrna import config  # noqa: E402
from sgrna.asymmetry import _side_of  # noqa: E402
from sgrna.run_ablation import cross_validate  # noqa: E402

# Import builders from the sibling script
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_intergenic_orientation import (  # noqa: E402
    OUT, bin100k_groups, build_matrix, load_intergenic_catalog, load_sets,
    locate_and_orient,
)

N_BOOT = 2000
SEED = 41
ALPHA = 0.05
MIN_N = 50


def _rho(y, p):
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    m = np.isfinite(y) & np.isfinite(p)
    if m.sum() < 3:
        return np.nan
    r = spearmanr(y[m], p[m]).statistic
    return float(r) if np.isfinite(r) else np.nan


def arm_columns(names, family_of, seed: int = 41):
    names = np.asarray(names)
    family_of = np.asarray(family_of)
    base = np.where(family_of == "base")[0]
    flank = np.where(family_of != "base")[0]
    side = np.array([_side_of(n) for n in names[flank]])
    up_idx = flank[side == "up"]
    dn_idx = flank[side == "dn"]
    rng = np.random.default_rng(seed)
    dn_matched = rng.choice(dn_idx, size=min(len(up_idx), len(dn_idx)),
                            replace=False)
    return {
        "upstream_only": np.sort(np.concatenate([base, up_idx])),
        "downstream_matched": np.sort(np.concatenate([base, dn_matched])),
    }, names, family_of


def oof_for_arm(X, y, names, family_of, cols, groups, tag: str, arm: str):
    res = cross_validate(
        X[:, cols], y, list(np.asarray(names)[cols]),
        np.asarray(family_of)[cols],
        seeds=(SEED,), groups=groups,
        tag=f"orientci_{tag}_{arm}", experiment=arm,
        verbose=True, return_oof=True,
    )
    return np.asarray(res["oof_pred"], float)


def bootstrap_delta(y_opp, up_opp, dn_opp, g_opp,
                    y_same, up_same, dn_same, g_same,
                    n_boot: int = N_BOOT, seed: int = 0):
    """Independent within-stratum group bootstrap of Δ = asym_opp - asym_same."""
    rng = np.random.default_rng(seed)

    def one_stratum(y, up, dn, groups):
        y = np.asarray(y, float)
        up = np.asarray(up, float)
        dn = np.asarray(dn, float)
        groups = np.asarray(groups, dtype=object)
        ok = np.isfinite(y) & np.isfinite(up) & np.isfinite(dn)
        y, up, dn, groups = y[ok], up[ok], dn[ok], groups[ok]
        uniq = np.unique(groups)
        point = _rho(y, dn) - _rho(y, up)
        stats = []
        for _ in range(n_boot):
            picked = rng.choice(uniq, size=len(uniq), replace=True)
            # expand with multiplicity
            parts_y, parts_up, parts_dn = [], [], []
            for g in picked:
                m = groups == g
                parts_y.append(y[m])
                parts_up.append(up[m])
                parts_dn.append(dn[m])
            yb = np.concatenate(parts_y)
            if len(yb) < MIN_N:
                continue
            stats.append(_rho(yb, np.concatenate(parts_dn))
                         - _rho(yb, np.concatenate(parts_up)))
        return point, np.asarray(stats, float), int(len(uniq))

    # Use separate RNGs so strata draws are independent but reproducible
    rng_opp = np.random.default_rng(seed)
    rng_same = np.random.default_rng(seed + 1)

    def boot_asym(y, up, dn, groups, rng):
        y = np.asarray(y, float); up = np.asarray(up, float); dn = np.asarray(dn, float)
        groups = np.asarray(groups, dtype=object)
        ok = np.isfinite(y) & np.isfinite(up) & np.isfinite(dn)
        y, up, dn, groups = y[ok], up[ok], dn[ok], groups[ok]
        uniq = np.unique(groups)
        point = _rho(y, dn) - _rho(y, up)
        stats = np.empty(n_boot, float)
        for i in range(n_boot):
            picked = rng.choice(uniq, size=len(uniq), replace=True)
            parts_y, parts_up, parts_dn = [], [], []
            for g in picked:
                m = groups == g
                parts_y.append(y[m]); parts_up.append(up[m]); parts_dn.append(dn[m])
            yb = np.concatenate(parts_y)
            if len(yb) < MIN_N:
                stats[i] = np.nan
            else:
                stats[i] = (_rho(yb, np.concatenate(parts_dn))
                            - _rho(yb, np.concatenate(parts_up)))
        return point, stats, int(len(uniq))

    p_opp, s_opp, n_g_opp = boot_asym(y_opp, up_opp, dn_opp, g_opp, rng_opp)
    p_same, s_same, n_g_same = boot_asym(y_same, up_same, dn_same, g_same, rng_same)
    delta_point = p_opp - p_same
    delta_stats = s_opp - s_same
    delta_stats = delta_stats[np.isfinite(delta_stats)]
    lo = float(np.quantile(delta_stats, ALPHA / 2))
    hi = float(np.quantile(delta_stats, 1 - ALPHA / 2))
    return dict(
        asymmetry_opposite=p_opp,
        asymmetry_same=p_same,
        delta=delta_point,
        delta_lo=lo,
        delta_hi=hi,
        delta_se=float(np.std(delta_stats, ddof=1)),
        n_boot=int(len(delta_stats)),
        n_groups_opposite=n_g_opp,
        n_groups_same=n_g_same,
        p_delta_le_0=float(np.mean(delta_stats <= 0)),
        alpha=ALPHA,
        method=("independent within-stratum 100kb group bootstrap of "
                "OOF Spearman asymmetry; delta = opp - same"),
    )


def run_set(tag: str):
    print(f"\n=== {tag} ===", flush=True)
    lib = load_intergenic_catalog()
    sets = load_sets(lib)
    idx = locate_and_orient(sets[tag])
    X, y, names, meta = build_matrix(idx, tag)
    family_of = meta["family_of"]
    groups = bin100k_groups(np.asarray(meta["left"]))
    relative = np.asarray(meta["relative"], dtype=object)
    arms, names, family_of = arm_columns(names, family_of, seed=SEED)

    oof = {}
    for stratum in ("opposite_gene", "same_as_gene"):
        mask = relative == stratum
        print(f"  stratum {stratum} n={mask.sum()}", flush=True)
        for arm, cols in arms.items():
            key = f"{stratum}:{arm}"
            print(f"  OOF {key}", flush=True)
            pred = oof_for_arm(
                X[mask], y[mask], names, family_of, cols, groups[mask],
                tag=f"{tag}_{stratum}", arm=arm)
            oof[key] = dict(
                y=y[mask], pred=pred, groups=groups[mask],
                n=int(mask.sum()),
            )

    # point asymmetries from OOF
    for stratum in ("opposite_gene", "same_as_gene"):
        a = (_rho(oof[f"{stratum}:downstream_matched"]["y"],
                  oof[f"{stratum}:downstream_matched"]["pred"])
             - _rho(oof[f"{stratum}:upstream_only"]["y"],
                    oof[f"{stratum}:upstream_only"]["pred"]))
        print(f"  OOF asymmetry_matched {stratum}: {a:.4f}", flush=True)

    ci = bootstrap_delta(
        oof["opposite_gene:downstream_matched"]["y"],
        oof["opposite_gene:upstream_only"]["pred"],
        oof["opposite_gene:downstream_matched"]["pred"],
        oof["opposite_gene:downstream_matched"]["groups"],
        oof["same_as_gene:downstream_matched"]["y"],
        oof["same_as_gene:upstream_only"]["pred"],
        oof["same_as_gene:downstream_matched"]["pred"],
        oof["same_as_gene:downstream_matched"]["groups"],
        n_boot=N_BOOT, seed=0,
    )
    ci.update(
        set=tag,
        n_opposite=oof["opposite_gene:downstream_matched"]["n"],
        n_same=oof["same_as_gene:downstream_matched"]["n"],
    )
    print(json.dumps(ci, indent=2), flush=True)
    return ci


def main():
    # Preconditions
    if not config.QCT_TABLES.exists():
        raise SystemExit(
            f"Missing {config.QCT_TABLES}. Run: "
            f"cd {ROOT} && python -m sgrna.build_guide_index")
    vmap = config.INTERIM / "v_column_mapping.json"
    if not vmap.exists():
        raise SystemExit(
            f"Missing {vmap}. Run: cd {ROOT} && python -m sgrna.decode_v_columns")

    results = {}
    for tag in ("s6_good", "s4_hq"):
        results[tag] = run_set(tag)

    out = OUT / "intergenic_orientation_delta_ci.json"
    out.write_text(json.dumps(results, indent=2))
    # flat CSV for the table
    rows = []
    for tag, r in results.items():
        rows.append(dict(
            set=tag,
            asymmetry_opposite=r["asymmetry_opposite"],
            asymmetry_same=r["asymmetry_same"],
            delta=r["delta"],
            delta_lo=r["delta_lo"],
            delta_hi=r["delta_hi"],
            delta_se=r["delta_se"],
            p_delta_le_0=r["p_delta_le_0"],
            n_opposite=r["n_opposite"],
            n_same=r["n_same"],
            n_boot=r["n_boot"],
        ))
    pd.DataFrame(rows).to_csv(OUT / "intergenic_orientation_delta_ci.csv",
                              index=False)
    print("\n=== delta CI summary ===", flush=True)
    print(pd.DataFrame(rows).to_string(index=False), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
