#!/usr/bin/env python3
"""Rebuild base+a_flank for Guo intergenic guides; run orientation-conditioned asymmetry.

Primary set: Data S6 Quality=Good promoter+RBS guides (~9k). Data S4 HQ is almost
entirely RBS (1 promoter / 3644 RBS), so it is a sensitivity check only.

Orientation = guide strand vs nearest gene's strand (same_as_gene / opposite_gene).
If transcription drove the PAM-side flank advantage, asymmetry should differ
(or flip) between those strata; if R-loop / PAM geometry drove it, asymmetry
should be similar in both.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("RSI09_ROOT", "/workspace/Runjia")).resolve()
os.environ["RSI09_ROOT"] = str(ROOT)
sys.path.insert(0, str(ROOT))

from sgrna import config, featurise, genome  # noqa: E402
from sgrna.asymmetry import _side_of  # noqa: E402
from sgrna.features import a_flank  # noqa: E402
from sgrna.run_ablation import cross_validate  # noqa: E402

OUT = Path("/workspace/Jacky/agent-notes/results")
GUO = Path("/workspace/Jacky/agent-notes/guo_tables")
CACHE = config.INTERIM / "intergenic_orientation"
OUT.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)
FLANK_NT = 1000


def load_intergenic_catalog() -> pd.DataFrame:
    prom = pd.read_excel(GUO / "Data S1.xlsx", sheet_name="promoter").assign(kind="promoter")
    rbs = pd.read_excel(GUO / "Data S1.xlsx", sheet_name="RBS").assign(kind="RBS")
    lib = pd.concat([prom, rbs], ignore_index=True).rename(
        columns={"sgRNAID": "gid", "nucleotide sequence": "seq"})
    lib["seq"] = lib["seq"].astype(str).str.upper()
    return lib[["gid", "seq", "kind"]]


def load_sets(lib: pd.DataFrame) -> dict[str, pd.DataFrame]:
    s4 = pd.read_excel(GUO / "Data S4.xlsx", sheet_name="Cas9").rename(
        columns={"sgRNAID": "gid", "score": "y"})
    s6 = pd.read_excel(GUO / "Data S6.xlsx", sheet_name="Cas9").rename(
        columns={"sgRNA": "gid", "Log2_normalized_change": "log2",
                 "Quality": "quality", "gene (promoter)": "locus"})
    s6_good = s6[s6["quality"].astype(str) == "Good"].copy()
    s6_good["y"] = s6_good["log2"].abs()

    sets = {}
    for name, scores in (("s6_good", s6_good), ("s4_hq", s4)):
        cols = ["gid", "y"] + (["locus"] if "locus" in scores.columns else [])
        m = lib.merge(scores[cols], on="gid", how="inner")
        m = m.drop_duplicates(subset=["seq"]).reset_index(drop=True)
        sets[name] = m
    return sets


def locate_and_orient(df: pd.DataFrame) -> pd.DataFrame:
    gseq = genome.load_genome()
    genes = pd.read_csv(config.GENE_TABLE, sep="\t")
    gene_mid = ((genes["left"] + genes["right"]) / 2).to_numpy(float)
    gene_ori = genes["ori"].astype(str).to_numpy()
    gene_name = genes["name"].astype(str).to_numpy()
    gene_left = genes["left"].to_numpy(float)
    gene_right = genes["right"].to_numpy(float)

    rows = []
    for _, r in df.iterrows():
        hits = genome.locate_protospacer(r["seq"], genome=gseq)
        if not hits:
            continue
        # prefer unique; else take first
        h = hits[0]
        if h.get("n_genome_copies", 1) > 1:
            # still usable for orientation audit but flanks ambiguous — skip
            continue
        ctx = genome.guide_context(h["left"], h["right"], h["strand"],
                                   FLANK_NT, genome=gseq)
        # nearest gene by distance to guide midpoint
        mid = (h["left"] + h["right"]) / 2
        # distance to gene body (0 if overlapping)
        dist = np.where(
            (mid >= gene_left) & (mid <= gene_right),
            0.0,
            np.minimum(np.abs(mid - gene_left), np.abs(mid - gene_right)),
        )
        j = int(np.argmin(dist))
        g_ori = gene_ori[j]
        relative = ("same_as_gene" if str(h["strand"]) == g_ori
                    else "opposite_gene")
        rows.append(dict(
            sgRNAID=r["gid"],
            gid=r["gid"],
            kind=r["kind"],
            y=float(r["y"]),
            protospacer=ctx["protospacer"],
            pam=ctx["pam"],
            upstream=ctx["upstream"],
            downstream=ctx["downstream"],
            strand=h["strand"],
            left=int(h["left"]),
            right=int(h["right"]),
            located=True,
            near_gene=gene_name[j],
            gene_ori=g_ori,
            gene_dist=float(dist[j]),
            relative=relative,
            locus=r.get("locus", gene_name[j]),
        ))
    out = pd.DataFrame(rows)
    return out


def build_matrix(idx: pd.DataFrame, tag: str, force: bool = False):
    x_path = CACHE / f"X_{tag}.npy"
    meta_path = CACHE / f"meta_{tag}.json"
    flank_path = CACHE / f"a_flank_{tag}.csv"
    idx_path = CACHE / f"index_{tag}.csv"

    if x_path.exists() and meta_path.exists() and not force:
        meta = json.loads(meta_path.read_text())
        X = np.load(x_path)
        print(f"  cached {tag}: {X.shape}", flush=True)
        return X, np.asarray(meta["y"], float), list(meta["names"]), meta

    print(f"  featurising base for {len(idx):,} guides ({tag})", flush=True)
    base, base_names = featurise.matrix(idx["protospacer"].tolist(), verbose=True)
    base_names = list(base_names)

    print(f"  building a_flank ({tag})", flush=True)
    block = a_flank.build(idx)
    block.to_csv(flank_path, index=False)
    block = block.set_index(config.ID_COL).reindex(idx[config.ID_COL])
    vals = block.to_numpy(dtype=np.float32)
    good = ~np.all(np.isnan(vals), axis=0)
    X = np.hstack([np.asarray(base, dtype=np.float32), vals[:, good]])
    names = base_names + [c for c, k in zip(block.columns, good) if k]
    family_of = ["base"] * len(base_names) + ["a_flank"] * (len(names) - len(base_names))

    meta = dict(
        y=idx["y"].astype(float).tolist(),
        names=names,
        family_of=family_of,
        n_base=len(base_names),
        ids=idx[config.ID_COL].astype(str).tolist(),
        left=idx["left"].astype(int).tolist(),
        relative=idx["relative"].astype(str).tolist(),
        strand=idx["strand"].astype(str).tolist(),
        kind=idx["kind"].astype(str).tolist(),
        near_gene=idx["near_gene"].astype(str).tolist(),
        tag=tag,
        n=int(len(idx)),
    )
    np.save(x_path, X)
    meta_path.write_text(json.dumps(meta))
    idx.to_csv(idx_path, index=False)
    print(f"  built {tag}: {X.shape[0]:,} x {X.shape[1]:,} "
          f"({len(base_names):,} base + {X.shape[1]-len(base_names)} flank)",
          flush=True)
    return X, np.asarray(meta["y"], float), names, meta


def bin100k_groups(lefts: np.ndarray) -> np.ndarray:
    return (np.asarray(lefts, float) // 100_000).astype(int).astype(str)


def run_asymmetry(X, y, names, family_of, groups, stratum_labels, tag: str,
                  seed: int = 41) -> pd.DataFrame:
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

    arms = {
        "baseline": base,
        "upstream_only": np.sort(np.concatenate([base, up_idx])),
        "downstream_only": np.sort(np.concatenate([base, dn_idx])),
        "downstream_matched": np.sort(np.concatenate([base, dn_matched])),
        "both": np.arange(X.shape[1]),
    }
    print(f"  flank columns: up={len(up_idx)} dn={len(dn_idx)} "
          f"matched={len(dn_matched)}", flush=True)

    rows = []
    strata = {"all": np.ones(len(y), dtype=bool)}
    for lab in sorted(set(stratum_labels)):
        strata[lab] = np.asarray(stratum_labels) == lab

    for sname, mask in strata.items():
        if mask.sum() < 200:
            print(f"  skip {sname}: n={mask.sum()}", flush=True)
            continue
        for arm, cols in arms.items():
            print(f"  CV {tag}/{sname}/{arm} n={mask.sum()}", flush=True)
            res = cross_validate(
                X[mask][:, cols], y[mask], list(names[cols]), family_of[cols],
                seeds=(seed,), groups=groups[mask],
                tag=f"intergenic_{tag}_{sname}_{arm}", experiment=arm,
                verbose=True)
            s = res["summary"]
            rows.append(dict(
                set=tag, stratum=sname, arm=arm, n=int(mask.sum()),
                spearman=s["spearman"]["mean"], r2=s["r2"]["mean"]))
    return pd.DataFrame(rows)


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (aset, stratum), sub in df.groupby(["set", "stratum"]):
        piv = sub.set_index("arm")["spearman"]
        b = float(piv.get("baseline", np.nan))
        row = dict(set=aset, stratum=stratum, n=int(sub["n"].iloc[0]),
                   baseline=b)
        for arm in ("upstream_only", "downstream_only",
                    "downstream_matched", "both"):
            if arm in piv.index:
                row[f"gain_{arm}"] = float(piv[arm] - b)
        if "gain_downstream_only" in row and "gain_upstream_only" in row:
            row["asymmetry"] = row["gain_downstream_only"] - row["gain_upstream_only"]
        if "gain_downstream_matched" in row and "gain_upstream_only" in row:
            row["asymmetry_matched"] = (
                row["gain_downstream_matched"] - row["gain_upstream_only"])
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    print("=== load catalogs ===", flush=True)
    lib = load_intergenic_catalog()
    sets = load_sets(lib)
    for k, v in sets.items():
        print(f"  {k}: {len(v)} guides "
              f"(promoter={(v.kind=='promoter').sum()}, "
              f"RBS={(v.kind=='RBS').sum()})", flush=True)

    # Primary = S6 Good; sensitivity = S4 HQ
    all_rows = []
    audit = {}
    for tag in ("s6_good", "s4_hq"):
        print(f"\n=== locate + orient ({tag}) ===", flush=True)
        idx = locate_and_orient(sets[tag])
        print(f"  located unique: {len(idx)}", flush=True)
        print(f"  strand: {idx['strand'].value_counts().to_dict()}", flush=True)
        print(f"  relative: {idx['relative'].value_counts().to_dict()}", flush=True)
        print(f"  kind: {idx['kind'].value_counts().to_dict()}", flush=True)
        audit[tag] = dict(
            n_input=int(len(sets[tag])),
            n_located=int(len(idx)),
            strand=idx["strand"].value_counts().astype(int).to_dict(),
            relative=idx["relative"].value_counts().astype(int).to_dict(),
            kind=idx["kind"].value_counts().astype(int).to_dict(),
            note=("S4 HQ intergenic is almost entirely RBS; S6 Good is the "
                  "primary orientation test set."),
        )
        if len(idx) < 200:
            print(f"  too few guides for {tag}; skip CV", flush=True)
            continue

        X, y, names, meta = build_matrix(idx, tag)
        family_of = meta["family_of"]
        groups = bin100k_groups(np.asarray(meta["left"]))
        print(f"  groups (100kb): {len(set(groups))}", flush=True)

        df = run_asymmetry(X, y, names, family_of, groups,
                           meta["relative"], tag=tag)
        all_rows.append(df)

    if not all_rows:
        raise SystemExit("no sets produced results")
    raw = pd.concat(all_rows, ignore_index=True)
    summ = summarise(raw)
    raw.to_csv(OUT / "intergenic_orientation_asymmetry.csv", index=False)
    summ.to_csv(OUT / "intergenic_orientation_summary.csv", index=False)
    (OUT / "intergenic_orientation_rebuild_audit.json").write_text(
        json.dumps(audit, indent=2))
    print("\n=== summary ===", flush=True)
    print(summ.to_string(index=False), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
