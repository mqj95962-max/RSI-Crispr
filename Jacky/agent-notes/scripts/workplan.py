#!/usr/bin/env python3
"""Jacky work-plan runner — lives under Jacky/agent-notes/scripts/ (not in sgrna/).

    export RSI09_ROOT=/workspace/Runjia
    python Jacky/agent-notes/scripts/workplan.py --all
"""

from __future__ import annotations

import argparse
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

from sgrna import config  # noqa: E402
from sgrna import evaluate as ev  # noqa: E402


OUT = Path("/workspace/Jacky/agent-notes/results")
OUT.mkdir(parents=True, exist_ok=True)

CROD_GENOME = ROOT / "external_data" / "C_rodentium" / "NC_013716.1.fasta"


def _spearman(a, b) -> float:
    if len(a) < 2:
        return float("nan")
    rho, _ = spearmanr(a, b)
    return 0.0 if np.isnan(rho) else float(rho)


def _genome_for(name: str):
    if name == "TevSpCas9":
        if not CROD_GENOME.exists():
            raise FileNotFoundError(CROD_GENOME)
        return str(CROD_GENOME)
    return None


def guo_ceiling(r2_replicate: float = 0.78, r2_tiling: float = 0.771) -> dict:
    r_rep = float(np.sqrt(r2_replicate))
    r_til = float(np.sqrt(r2_tiling))
    sb_rep = 2 * r_rep / (1 + r_rep)
    out = dict(
        source="Guo et al. 2018 NAR Fig 2b/2c (verify against PDF)",
        r2_replicate_fig2b=r2_replicate,
        r2_tiling_fig2c=r2_tiling,
        pearson_replicate=r_rep,
        pearson_tiling=r_til,
        ceiling_one_replicate_proxy=float(np.sqrt(r_rep)),
        ceiling_two_replicate_mean_sb=float(np.sqrt(sb_rep)),
        ceiling_tiling_independent_library=float(np.sqrt(r_til)),
        slicer_rho_range="0.707–0.721",
        fraction_of_sb_ceiling=0.721 / float(np.sqrt(sb_rep)),
        note=("Attenuation is Pearson additive-noise; √Spearman is an approximation. "
              "within_gene_variance_share from Runjia = 0.8451."),
    )
    (OUT / "guo_reliability_ceiling.json").write_text(json.dumps(out, indent=2))
    return out


def flank_mappability_stratum(seed: int = 41, uniq_quantile: float = 0.75,
                              verbose: bool = True) -> pd.DataFrame:
    from sgrna.build_matrix import load_dataset
    from sgrna.make_diagnostic_blocks import build_mapgc
    from sgrna.run_ablation import cross_validate, load_groups

    X, y, names, ids, fam = load_dataset(["a_flank"], verbose=verbose)
    block = build_mapgc(verbose=verbose)
    uniq_cols = [c for c in block.columns if "eng.mapgc.uniq." in c]
    u = block.set_index(config.ID_COL).reindex(ids)[uniq_cols].mean(axis=1).to_numpy()
    thr = float(np.nanquantile(u[np.isfinite(u)], uniq_quantile))
    keep_hi = np.isfinite(u) & (u >= thr)
    keep_all = np.isfinite(y)

    if verbose:
        print(f"  high-uniqueness: {keep_hi.sum():,} / {len(keep_hi):,} "
              f"(mean uniq >= {thr:.4f})")

    base = np.where(np.asarray(fam) == "base")[0]
    groups = load_groups(ids, by="gene")
    rows = []
    for stratum, mask in (("all", keep_all), ("high_uniq", keep_hi)):
        for label, cols in (("baseline", base), ("+a_flank", np.arange(X.shape[1]))):
            res = cross_validate(
                X[mask][:, cols], y[mask], list(np.asarray(names)[cols]),
                np.asarray(fam)[cols], seeds=(seed,), groups=groups[mask],
                tag="workplan_flank_stratum",
                experiment=f"{stratum}|{label}", verbose=verbose,
            )
            s = res["summary"]
            rows.append(dict(stratum=stratum, arm=label, n=int(mask.sum()),
                             spearman=s["spearman"]["mean"],
                             r2=s["r2"]["mean"]))
    df = pd.DataFrame(rows)
    # deltas within stratum
    deltas = []
    for stratum, sub in df.groupby("stratum"):
        b = sub.loc[sub.arm == "baseline", "spearman"].iloc[0]
        a = sub.loc[sub.arm == "+a_flank", "spearman"].iloc[0]
        deltas.append(dict(stratum=stratum, delta_spearman=a - b,
                           baseline=b, with_flank=a))
    ddf = pd.DataFrame(deltas)
    df.to_csv(OUT / "flank_high_uniqueness_stratum.csv", index=False)
    ddf.to_csv(OUT / "flank_stratum_deltas.csv", index=False)
    if verbose:
        print(ddf.to_string(index=False))
    return df


def replichore_asymmetry(seed: int = 41, verbose: bool = True) -> pd.DataFrame:
    from sgrna import genome
    from sgrna.asymmetry import _side_of
    from sgrna.build_matrix import load_dataset
    from sgrna.io_utils import load_guide_index
    from sgrna.run_ablation import cross_validate

    X, y, names, ids, fam = load_dataset(["a_flank"], verbose=verbose)
    names = np.asarray(names)
    fam = np.asarray(fam)
    base = np.where(fam == "base")[0]
    flank = np.where(fam != "base")[0]
    side = np.array([_side_of(n) for n in names[flank]])
    idx = load_guide_index().set_index(config.ID_COL).reindex(ids)
    repl = idx["left"].map(
        lambda p: genome.replichore(int(p)) if pd.notna(p) and np.isfinite(p) else None)

    rng = np.random.default_rng(seed)
    up_idx = flank[side == "up"]
    dn_idx = flank[side == "dn"]
    dn_matched = rng.choice(dn_idx, size=len(up_idx), replace=False)

    arms = {
        "baseline": base,
        "upstream_only": np.sort(np.concatenate([base, up_idx])),
        "downstream_only": np.sort(np.concatenate([base, dn_idx])),
        "downstream_matched": np.sort(np.concatenate([base, dn_matched])),
        "both": np.arange(X.shape[1]),
    }
    rows = []
    for rep in ("left", "right"):
        keep = (repl == rep).to_numpy()
        if keep.sum() < 500:
            continue
        for arm, cols in arms.items():
            res = cross_validate(
                X[keep][:, cols], y[keep], list(names[cols]), fam[cols],
                seeds=(seed,), tag=f"workplan_replichore_{rep}",
                experiment=arm, verbose=False,
            )
            s = res["summary"]
            rows.append(dict(replichore=rep, arm=arm, n=int(keep.sum()),
                             spearman=s["spearman"]["mean"]))
            if verbose:
                print(f"    {rep:5s} {arm:18s} rho {s['spearman']['mean']:.4f}")
    df = pd.DataFrame(rows)
    # asymmetry summary
    piv = df.pivot_table(index="replichore", columns="arm", values="spearman")
    for arm in ("upstream_only", "downstream_only", "downstream_matched", "both"):
        if arm in piv.columns and "baseline" in piv.columns:
            piv[f"gain_{arm}"] = piv[arm] - piv["baseline"]
    if {"gain_downstream_only", "gain_upstream_only"} <= set(piv.columns):
        piv["asymmetry"] = piv["gain_downstream_only"] - piv["gain_upstream_only"]
    if {"gain_downstream_matched", "gain_upstream_only"} <= set(piv.columns):
        piv["asymmetry_matched"] = (
            piv["gain_downstream_matched"] - piv["gain_upstream_only"])
    df.to_csv(OUT / "asymmetry_replichore.csv", index=False)
    piv.round(4).to_csv(OUT / "asymmetry_replichore_summary.csv")
    if verbose:
        print(piv.round(4).to_string())
    return df


def within_gene_metrics(seeds=(41,), verbose: bool = True) -> dict:
    from sgrna.build_matrix import load_dataset
    from sgrna.io_utils import load_guide_index
    from sgrna.run_ablation import cross_validate, load_groups

    X, y, names, ids, fam = load_dataset([], verbose=verbose)
    gene = (load_guide_index().set_index(config.ID_COL)
            .reindex(ids)["gene_name"].to_numpy(dtype=object))

    res_base = cross_validate(
        X, y, names, fam, seeds=seeds, groups=load_groups(ids, "gene"),
        tag="workplan_within_gene", experiment="baseline",
        return_oof=True, verbose=verbose,
    )
    pred, y_oof = res_base["oof_pred"], res_base["oof_y"]
    ok = np.isfinite(pred) & np.isfinite(y_oof)
    global_r = _spearman(y_oof[ok], pred[ok])
    wg = ev.within_gene_spearman(y_oof[ok], pred[ok], gene[ok], min_guides=10)
    pick = ev.expected_pick_percentile(y_oof[ok], pred[ok])
    pick_g = ev.expected_pick_percentile_within_gene(
        y_oof[ok], pred[ok], gene[ok])

    Xf, yf, nf, idf, ff = load_dataset(["a_flank"], verbose=False)
    gene_f = (load_guide_index().set_index(config.ID_COL)
              .reindex(idf)["gene_name"].to_numpy(dtype=object))
    res_fl = cross_validate(
        Xf, yf, nf, ff, seeds=seeds, groups=load_groups(idf, "gene"),
        tag="workplan_within_gene", experiment="+a_flank",
        return_oof=True, verbose=verbose,
    )
    pf = res_fl["oof_pred"]
    okf = np.isfinite(pf) & np.isfinite(yf)
    global_fl = _spearman(yf[okf], pf[okf])
    wg_fl = ev.within_gene_spearman(yf[okf], pf[okf], gene_f[okf], min_guides=10)

    ceiling_path = Path("/workspace/Jacky/agent-notes/results/ceiling_from_runjia.json")
    label_share = None
    if ceiling_path.exists():
        label_share = json.loads(ceiling_path.read_text())["within_gene_variance_share"]

    out = dict(
        global_spearman_oof_baseline=global_r,
        global_spearman_oof_flank=global_fl,
        within_gene_baseline={k: wg[k] for k in
                              ("n_genes", "mean_rho", "sd_rho", "median_rho")},
        within_gene_with_flanks={k: wg_fl[k] for k in
                                 ("n_genes", "mean_rho", "sd_rho", "median_rho")},
        delta_within_gene_spearman=(
            wg_fl["mean_rho"] - wg["mean_rho"]
            if np.isfinite(wg_fl.get("mean_rho", np.nan))
            and np.isfinite(wg.get("mean_rho", np.nan)) else None),
        pick_percentile_library=pick,
        pick_percentile_within_gene=pick_g,
        label_within_gene_variance_share=label_share,
    )
    (OUT / "within_gene_metrics.json").write_text(json.dumps(out, indent=2))
    if verbose:
        print(json.dumps(out, indent=2))
    return out


def bootstrap_transfer(seed: int = 41, verbose: bool = True) -> pd.DataFrame:
    import xgboost as xgb
    from sgrna.run_ablation import _impute, make_selector
    from sgrna.transfer import build

    pairs = (
        ("WT-SpCas9", "WT-SpCas9"),
        ("WT-SpCas9", "TevSpCas9"),
        ("TevSpCas9", "WT-SpCas9"),
        ("TevSpCas9", "TevSpCas9"),
    )
    rows = []
    for source, target in pairs:
        if verbose:
            print(f"\n== {source} -> {target} ==")
        Xs, ys, names = build(source, genome_fasta=_genome_for(source), verbose=verbose)
        Xt, yt, _ = build(target, genome_fasta=_genome_for(target), verbose=verbose)
        meta_t = json.loads(
            (config.INTERIM / f"transfer_{target}" / "meta.json").read_text())
        meta_s = json.loads(
            (config.INTERIM / f"transfer_{source}" / "meta.json").read_text())
        # group by 100 kb genomic bins from meta left when available
        left = np.asarray(meta_t.get("left") or [np.nan] * len(yt), float)
        if np.isfinite(left).sum() > 100:
            groups = (left // 100_000).astype(object)
            groups = np.where(np.isfinite(left), groups, "unk")
        else:
            groups = np.arange(len(yt))  # row fallback

        A, B = _impute(np.asarray(Xs, float), np.asarray(Xt, float))
        sel = make_selector(seed)
        sel.fit(A, ys)
        top = np.argsort(sel.feature_importances_)[::-1][: config.N_FEATURES]
        m = xgb.XGBRegressor(
            n_estimators=400, learning_rate=0.05, max_depth=4,
            subsample=0.8, colsample_bytree=0.8, random_state=seed,
            n_jobs=-1, tree_method="hist")
        m.fit(A[:, top], ys)
        pred = np.asarray(m.predict(B[:, top])).ravel()
        ci = ev.bootstrap_interval_grouped(
            yt, pred, groups, _spearman, n_boot=2000, seed=seed)
        src_p = set(meta_s.get("protospacers") or [])
        tgt_p = set(meta_t.get("protospacers") or [])
        overlap = (len(src_p & tgt_p) / len(tgt_p)) if tgt_p else float("nan")
        rows.append(dict(
            source=source, target=target, seed=seed,
            spearman=ci["point"], ci_lo=ci["lo"], ci_hi=ci["hi"],
            ci_se=ci.get("se"), guide_overlap=overlap, n_test=len(yt),
            n_groups=ci.get("n_groups"),
        ))
        if verbose:
            print(f"  rho {ci['point']:.4f} [{ci['lo']:.4f}, {ci['hi']:.4f}] "
                  f"overlap={overlap:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "transfer_bootstrap.csv", index=False)
    return df


def bootstrap_cross_kingdom(seed: int = 41, n_boot: int = 2000,
                            verbose: bool = True) -> dict:
    """Shared-column transfer E. coli ↔ human with bootstrap CIs + GC quadratic."""
    from sklearn.linear_model import LinearRegression
    import xgboost as xgb
    from sgrna.io_utils import load_base_matrix
    from sgrna.run_ablation import _impute

    human_path = config.RAW / "human_feature_matrix.csv"
    if not human_path.exists():
        msg = dict(skipped=True, reason=f"missing {human_path}")
        (OUT / "cross_kingdom_bootstrap.json").write_text(json.dumps(msg, indent=2))
        return msg

    if verbose:
        print("  loading E. coli and human matrices (large) ...")
    eco = load_base_matrix()
    hum = pd.read_csv(human_path, header=config.BASE_MATRIX_HEADER_ROW, low_memory=False)

    # drop non-features
    drop = {config.TARGET_COL, config.ID_COL, "sequence", "Sequence", "sgRNA",
            "sgrna", "guide", "Guide"}
    y_e = pd.to_numeric(eco[config.TARGET_COL], errors="coerce")
    # human target column name may also be cut.score
    y_h_col = config.TARGET_COL if config.TARGET_COL in hum.columns else "cut.score"
    y_h = pd.to_numeric(hum[y_h_col], errors="coerce")

    feat_e = [c for c in eco.columns if c not in drop]
    feat_h = [c for c in hum.columns if c not in drop]
    shared = sorted(set(feat_e) & set(feat_h))
    if verbose:
        print(f"  shared feature columns: {len(shared):,}")

    Xe = eco[shared].apply(pd.to_numeric, errors="coerce").to_numpy(np.float32)
    Xh = hum[shared].apply(pd.to_numeric, errors="coerce").to_numpy(np.float32)
    ye = y_e.to_numpy(float)
    yh = y_h.to_numpy(float)
    ok_e = np.isfinite(ye) & np.isfinite(Xe).any(axis=1)
    ok_h = np.isfinite(yh) & np.isfinite(Xh).any(axis=1)
    Xe, ye = Xe[ok_e], ye[ok_e]
    Xh, yh = Xh[ok_h], yh[ok_h]

    def fit_predict(Xtr, ytr, Xte, seed):
        A, B = _impute(Xtr, Xte)
        # subsample features if huge
        n_feat = min(300, A.shape[1])
        # cheap selection via variance then xgb on subset if needed
        m = xgb.XGBRegressor(
            n_estimators=200, learning_rate=0.05, max_depth=4,
            subsample=0.8, colsample_bytree=0.3, random_state=seed,
            n_jobs=-1, tree_method="hist")
        # select top 300 by a quick fit on random column sample for speed
        rng = np.random.default_rng(seed)
        cand = rng.choice(A.shape[1], size=min(2000, A.shape[1]), replace=False)
        m.fit(A[:, cand], ytr)
        top_local = np.argsort(m.feature_importances_)[::-1][:n_feat]
        top = cand[top_local]
        m2 = xgb.XGBRegressor(
            n_estimators=400, learning_rate=0.05, max_depth=4,
            subsample=0.8, colsample_bytree=0.8, random_state=seed,
            n_jobs=-1, tree_method="hist")
        m2.fit(A[:, top], ytr)
        return np.asarray(m2.predict(B[:, top])).ravel()

    if verbose:
        print("  fitting E.coli -> human ...")
    pred_eh = fit_predict(Xe, ye, Xh, seed)
    if verbose:
        print("  fitting human -> E.coli ...")
    pred_he = fit_predict(Xh, yh, Xe, seed)
    # same-species controls (random split)
    rng = np.random.default_rng(seed)
    idx_e = rng.permutation(len(ye))
    tr_e, te_e = idx_e[: int(0.8 * len(ye))], idx_e[int(0.8 * len(ye)):]
    pred_ee = fit_predict(Xe[tr_e], ye[tr_e], Xe[te_e], seed)
    idx_h = rng.permutation(len(yh))
    tr_h, te_h = idx_h[: int(0.8 * len(yh))], idx_h[int(0.8 * len(yh)):]
    pred_hh = fit_predict(Xh[tr_h], yh[tr_h], Xh[te_h], seed)

    def pack(y, p, name):
        ci = ev.bootstrap_interval(y, p, _spearman, n_boot=n_boot, seed=seed)
        return dict(name=name, **ci)

    transfers = [
        pack(yh, pred_eh, "ecoli_to_human"),
        pack(ye, pred_he, "human_to_ecoli"),
        pack(ye[te_e], pred_ee, "ecoli_to_ecoli"),
        pack(yh[te_h], pred_hh, "human_to_human"),
    ]

    # GC quadratic on both
    def gc_curve(X, y, shared_names, label):
        # find a GC-like column if present
        gc_cols = [c for c in shared_names if "gc" in c.lower() or "GC" in c]
        # Noshay often has composition elsewhere; use mean of one-hot G+C if needed
        # Fall back: estimate GC from monomer one-hots if available
        if not gc_cols:
            # approximate from positions if electron/onehot not tractable — skip
            return dict(label=label, skipped=True, reason="no GC column in shared set")
        # use first gc-like column as proxy; better: average if many
        g = np.nanmean(X[:, [shared_names.index(c) for c in gc_cols[:5]]], axis=1)
        ok = np.isfinite(g) & np.isfinite(y)
        g, yy = g[ok], y[ok]
        # linear + quadratic
        Z1 = g.reshape(-1, 1)
        Z2 = np.c_[g, g ** 2]
        r1 = LinearRegression().fit(Z1, yy)
        r2 = LinearRegression().fit(Z2, yy)
        pred1, pred2 = r1.predict(Z1), r2.predict(Z2)
        # turning point for ax^2+bx+c: -b/(2a)
        a, b = r2.coef_[1], r2.coef_[0]
        turn = float(-b / (2 * a)) if a != 0 else float("nan")
        # bootstrap turning point
        rng2 = np.random.default_rng(seed)
        turns = []
        for _ in range(1000):
            ix = rng2.integers(0, len(g), size=len(g))
            rr = LinearRegression().fit(np.c_[g[ix], g[ix] ** 2], yy[ix])
            aa, bb = rr.coef_[1], rr.coef_[0]
            if aa != 0:
                turns.append(-bb / (2 * aa))
        turns = np.asarray(turns)
        return dict(
            label=label, n=int(len(g)),
            linear_rho=_spearman(yy, pred1),
            quadratic_rho=_spearman(yy, pred2),
            linear_r2=float(1 - np.mean((yy - pred1) ** 2) / np.var(yy)),
            quadratic_r2=float(1 - np.mean((yy - pred2) ** 2) / np.var(yy)),
            turning_point=turn,
            turning_point_ci=[float(np.quantile(turns, 0.025)),
                              float(np.quantile(turns, 0.975))] if len(turns) else None,
            n_gc_cols_used=len(gc_cols[:5]),
        )

    gc_e = gc_curve(Xe, ye, shared, "ecoli")
    gc_h = gc_curve(Xh, yh, shared, "human")

    out = dict(
        n_shared_features=len(shared),
        n_ecoli=int(len(ye)), n_human=int(len(yh)),
        transfers=transfers,
        gc_quadratic=dict(ecoli=gc_e, human=gc_h),
    )
    (OUT / "cross_kingdom_bootstrap.json").write_text(json.dumps(out, indent=2))
    if verbose:
        for t in transfers:
            print(f"  {t['name']:20s} rho {t['point']:.4f} "
                  f"[{t['lo']:.4f}, {t['hi']:.4f}]")
    return out


def label_provenance_check(verbose: bool = True) -> dict:
    from sgrna.diagnose import ceiling
    out = {
        "cut_score_chain": (
            "Noshay cut.score ← Guo 2018 Cas9/dCas9 depletion; crisprHAL "
            "33,567 read-count-filtered. Absolute Z-score after geometric mean "
            "of two replicates. See CUT_SCORE_PROVENANCE.md"
        ),
        "within_gene_variance_share_from_runjia": 0.8451212917093656,
    }
    try:
        out["diagnose_ceiling"] = ceiling(verbose=verbose)
    except Exception as exc:
        out["diagnose_ceiling_error"] = str(exc)
    (OUT / "label_provenance.json").write_text(json.dumps(out, indent=2, default=str))
    return out


def inventory_guo_tables(verbose: bool = True) -> dict:
    """Inspect Guo supplementary xlsx if present."""
    path = Path("/workspace/Jacky/agent-notes/guo_tables/Guo2018_Supp_tables_1-10.xlsx")
    out = dict(path=str(path), exists=path.exists(), sheets={})
    if not path.exists():
        (OUT / "guo_tables_inventory.json").write_text(json.dumps(out, indent=2))
        return out
    try:
        xl = pd.ExcelFile(path)
    except Exception as exc:
        out["error"] = str(exc)
        (OUT / "guo_tables_inventory.json").write_text(json.dumps(out, indent=2))
        return out
    for sheet in xl.sheet_names:
        try:
            df = pd.read_excel(path, sheet_name=sheet, nrows=5)
            cols = [str(c) for c in df.columns]
            # also peek full for replicate-like names
            df0 = pd.read_excel(path, sheet_name=sheet, nrows=0)
            all_cols = [str(c) for c in df0.columns]
            rep_hits = [c for c in all_cols
                        if any(k in c.lower() for k in
                               ("replicat", "rep1", "rep2", "rep_1", "rep_2",
                                "count", "read"))]
            out["sheets"][sheet] = dict(
                n_columns=len(all_cols),
                columns_head=all_cols[:40],
                replicate_like_columns=rep_hits,
            )
            if verbose:
                print(f"  sheet {sheet!r}: {len(all_cols)} cols; "
                      f"replicate-like={rep_hits[:10]}")
        except Exception as exc:
            out["sheets"][sheet] = dict(error=str(exc))
    (OUT / "guo_tables_inventory.json").write_text(json.dumps(out, indent=2))
    return out


def ensure_cache(verbose: bool = True):
    from sgrna.build_matrix import build_cache
    build_cache(force=False, verbose=verbose)


def run_all(verbose: bool = True):
    print("=== Guo ceiling ===")
    print(json.dumps(guo_ceiling(), indent=2))
    print("=== Guo table inventory ===")
    inventory_guo_tables(verbose=verbose)
    print("=== ensure base cache ===")
    ensure_cache(verbose=verbose)
    print("=== label provenance ===")
    label_provenance_check(verbose=verbose)
    print("=== within-gene ===")
    within_gene_metrics(verbose=verbose)
    print("=== flank mappability stratum ===")
    flank_mappability_stratum(verbose=verbose)
    print("=== replichore asymmetry ===")
    replichore_asymmetry(verbose=verbose)
    print("=== transfer bootstrap ===")
    bootstrap_transfer(verbose=verbose)
    print("=== cross-kingdom bootstrap ===")
    bootstrap_cross_kingdom(verbose=verbose)
    print("=== DONE ===")
    print("results in", OUT)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--describe", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--guo-ceiling", action="store_true")
    ap.add_argument("--guo-inventory", action="store_true")
    ap.add_argument("--cache", action="store_true")
    ap.add_argument("--flank-stratum", action="store_true")
    ap.add_argument("--replichore", action="store_true")
    ap.add_argument("--within-gene", action="store_true")
    ap.add_argument("--bootstrap-transfer", action="store_true")
    ap.add_argument("--bootstrap-cross-kingdom", action="store_true")
    ap.add_argument("--label-provenance", action="store_true")
    args = ap.parse_args(argv)
    if args.describe:
        print(config.describe())
        print("outputs ->", OUT)
        return 0
    if args.all:
        run_all()
        return 0
    if args.guo_ceiling:
        print(json.dumps(guo_ceiling(), indent=2)); return 0
    if args.guo_inventory:
        print(json.dumps(inventory_guo_tables(), indent=2)); return 0
    if args.cache:
        ensure_cache(); return 0
    if args.flank_stratum:
        print(flank_mappability_stratum()); return 0
    if args.replichore:
        print(replichore_asymmetry()); return 0
    if args.within_gene:
        print(json.dumps(within_gene_metrics(), indent=2)); return 0
    if args.bootstrap_transfer:
        print(bootstrap_transfer()); return 0
    if args.bootstrap_cross_kingdom:
        print(json.dumps(bootstrap_cross_kingdom(), indent=2)); return 0
    if args.label_provenance:
        print(json.dumps(label_provenance_check(), indent=2, default=str)); return 0
    ap.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
