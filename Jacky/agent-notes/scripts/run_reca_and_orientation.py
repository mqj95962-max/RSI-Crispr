#!/usr/bin/env python3
"""ΔrecA flank ablation + intergenic strand orientation audit."""
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

from sgrna import config, genome  # noqa: E402
from sgrna.build_matrix import load_dataset  # noqa: E402
from sgrna.io_utils import load_guide_index  # noqa: E402
from sgrna.run_ablation import cross_validate, load_groups  # noqa: E402

OUT = Path("/workspace/Jacky/agent-notes/results")
GUO = Path("/workspace/Jacky/agent-notes/guo_tables")
OUT.mkdir(parents=True, exist_ok=True)


def load_joins() -> pd.DataFrame:
    idx = load_guide_index().copy()
    idx["gid"] = idx["sgRNAID"].astype(str).str.replace(r"_Cas9$", "", regex=True)
    s4_wt = pd.read_excel(GUO / "Data S4.xlsx", sheet_name="Cas9").rename(
        columns={"sgRNAID": "gid", "score": "s4_cas9"})
    s4_reca = pd.read_excel(GUO / "Data S4.xlsx", sheet_name="ΔrecA Cas9").rename(
        columns={"sgRNAID": "gid", "score": "s4_reca"})
    s6_wt = pd.read_excel(GUO / "Data S6.xlsx", sheet_name="Cas9").rename(
        columns={"sgRNA": "gid", "Log2_normalized_change": "s6_cas9",
                 "Quality": "quality"})
    s6_reca = pd.read_excel(GUO / "Data S6.xlsx", sheet_name="recA(-) Cas9").rename(
        columns={"sgRNA": "gid", "Log2_normalized_change": "s6_reca"})
    m = idx.merge(s4_wt[["gid", "s4_cas9"]], on="gid", how="left")
    m = m.merge(s4_reca[["gid", "s4_reca"]], on="gid", how="left")
    m = m.merge(s6_wt[["gid", "s6_cas9", "quality", "gene (promoter)"]],
                on="gid", how="left")
    m = m.merge(s6_reca[["gid", "s6_reca"]], on="gid", how="left")
    return m


def run_reca(m: pd.DataFrame) -> pd.DataFrame:
    print("=== loading matrix ===", flush=True)
    X, y, names, ids, fam = load_dataset(["a_flank"], verbose=True)
    m2 = m.set_index(config.ID_COL).reindex(ids)
    base = np.where(np.asarray(fam) == "base")[0]
    groups = load_groups(ids, by="gene")

    def one(label, yvec):
        mask = np.isfinite(yvec)
        rows = []
        for arm, cols in (("baseline", base), ("+a_flank", np.arange(X.shape[1]))):
            print(f"  CV {label} {arm} n={mask.sum()}", flush=True)
            res = cross_validate(
                X[mask][:, cols], yvec[mask], list(np.asarray(names)[cols]),
                np.asarray(fam)[cols], seeds=(41,), groups=groups[mask],
                tag=f"workplan_reca_{label}", experiment=arm, verbose=True)
            s = res["summary"]
            rows.append(dict(label=label, arm=arm, n=int(mask.sum()),
                             spearman=s["spearman"]["mean"], r2=s["r2"]["mean"]))
        return rows

    rows = []
    rows += one("published_cut_score", y)
    rows += one("guo_S4_Cas9", m2["s4_cas9"].to_numpy(float))
    rows += one("guo_S4_deltaRecA", m2["s4_reca"].to_numpy(float))
    rows += one("guo_S6_Cas9_abs", np.abs(m2["s6_cas9"].to_numpy(float)))
    rows += one("guo_S6_deltaRecA_abs", np.abs(m2["s6_reca"].to_numpy(float)))
    df = pd.DataFrame(rows)
    summ = []
    for lab, sub in df.groupby("label"):
        b = sub.loc[sub.arm == "baseline", "spearman"].iloc[0]
        a = sub.loc[sub.arm == "+a_flank", "spearman"].iloc[0]
        summ.append(dict(label=lab, baseline=b, with_flank=a, delta=a - b,
                         delta_over_baseline=(a - b) / b if b else np.nan,
                         n=int(sub["n"].iloc[0])))
    sdf = pd.DataFrame(summ)
    df.to_csv(OUT / "reca_flank_ablation.csv", index=False)
    sdf.to_csv(OUT / "reca_flank_summary.csv", index=False)
    print(sdf.to_string(index=False), flush=True)
    return sdf


def orientation_audit() -> dict:
    """Can intergenic guides support a strand-orientation test?"""
    print("=== intergenic orientation audit ===", flush=True)
    prom = pd.read_excel(GUO / "Data S1.xlsx", sheet_name="promoter")
    rbs = pd.read_excel(GUO / "Data S1.xlsx", sheet_name="RBS")
    coding = pd.read_excel(GUO / "Data S1.xlsx", sheet_name="CRISPRi (gene-targeting)")
    s6 = pd.read_excel(GUO / "Data S6.xlsx", sheet_name="Cas9").rename(
        columns={"sgRNA": "gid", "Log2_normalized_change": "score", "Quality": "quality"})
    s4 = pd.read_excel(GUO / "Data S4.xlsx", sheet_name="Cas9").rename(
        columns={"sgRNAID": "gid"})

    def prep(df, kind):
        d = df.rename(columns={"sgRNAID": "gid",
                               "nucleotide sequence": "seq"}).copy()
        d["seq"] = d["seq"].astype(str).str.upper()
        d["kind"] = kind
        return d[["gid", "seq", "kind"]]

    lib = pd.concat([prep(prom, "promoter"), prep(rbs, "RBS"),
                     prep(coding, "coding")], ignore_index=True)
    lib = lib.merge(s6[["gid", "score", "quality", "gene (promoter)"]],
                    on="gid", how="left")
    lib["in_S4"] = lib["gid"].isin(set(s4["gid"].astype(str)))

    # locate on genome (use locate_protospacer; falls back to find)
    g = genome.load_genome()

    def locate(seq: str):
        if len(seq) != 20 or any(c not in "ACGT" for c in seq):
            return None, None
        hits = genome.locate_protospacer(seq, genome=g)
        if hits:
            h = hits[0]
            return int(h["left"]), str(h["strand"])
        i = g.find(seq)
        if i >= 0:
            return i + 1, "+"
        r = genome.revcomp(seq)
        i = g.find(r)
        if i >= 0:
            return i + 1, "-"
        return None, None

    inter = lib[lib.kind.isin(["promoter", "RBS"])].copy()
    locs = [locate(s) for s in inter["seq"]]
    inter["left"] = [a for a, _ in locs]
    inter["strand"] = [b for _, b in locs]
    located = inter[inter["left"].notna()].copy()

    genes = pd.read_csv(config.GENE_TABLE, sep="\t")
    # gene table column names vary; normalise
    rename = {}
    for a, b in (("left", "left"), ("start", "left"), ("begin", "left"),
                 ("right", "right"), ("end", "right"),
                 ("ori", "ori"), ("strand", "ori"), ("orientation", "ori"),
                 ("name", "name"), ("gene", "name"), ("locus_tag", "name")):
        if a in genes.columns and b not in rename.values():
            rename[a] = b
    genes = genes.rename(columns=rename)
    gene_mid = ((genes["left"].astype(float) + genes["right"].astype(float)) / 2).to_numpy()
    gene_ori = genes["ori"].astype(str).to_numpy()
    gene_name = genes["name"].astype(str).to_numpy()

    def nearest(pos):
        j = int(np.argmin(np.abs(gene_mid - pos)))
        return gene_name[j], gene_ori[j]

    nearests = [nearest(p) for p in located["left"]]
    located["near_gene"] = [a for a, _ in nearests]
    located["gene_ori"] = [b for _, b in nearests]
    located["relative"] = np.where(
        located["strand"].astype(str) == located["gene_ori"].astype(str),
        "same_as_gene", "opposite_gene")

    # among those with S6 scores: strand balance + score-by-strand (sanity)
    scored = located[located["score"].notna()].copy()
    strand_score = {}
    if len(scored):
        for st, sub in scored.groupby("strand"):
            strand_score[str(st)] = dict(
                n=int(len(sub)),
                mean_abs_score=float(np.nanmean(np.abs(sub["score"]))),
                median_abs_score=float(np.nanmedian(np.abs(sub["score"]))),
            )

    coding_seqs = set(coding["nucleotide sequence"].astype(str).str.upper())
    ours = set(load_guide_index()["protospacer"].astype(str).str.upper())
    inter_seqs = set(inter["seq"].astype(str).str.upper())

    summary = dict(
        n_intergenic_designed=int(len(inter)),
        n_intergenic_in_S6=int(inter["score"].notna().sum()),
        n_intergenic_in_S4_HQ=int(inter["in_S4"].sum()),
        n_intergenic_located_on_genome=int(len(located)),
        strand_counts=located["strand"].value_counts(dropna=False).astype(int).to_dict(),
        relative_to_nearest_gene=located["relative"].value_counts().astype(int).to_dict(),
        by_kind_strand={f"{k[0]}|{k[1]}": int(v)
                        for k, v in located.groupby(["kind", "strand"]).size().items()},
        scored_strand_abs_activity=strand_score,
        coding_in_our_matrix=int(len(coding_seqs & ours)),
        intergenic_seq_overlap_with_13880=int(len(inter_seqs & ours)),
        note=(
            "Our Noshay 13,880 matrix is almost entirely the coding/CRISPRi "
            "sub-library. Intergenic promoter/RBS guides exist in Guo Data "
            "S1/S4/S6 with BOTH strands, but they are not in the published "
            "Noshay matrix, so the asymmetry ablation cannot be re-run on them "
            "without rebuilding features for those guides."
        ),
    )
    # JSON-safe keys
    summary["strand_counts"] = {str(k): int(v) for k, v in summary["strand_counts"].items()}
    summary["relative_to_nearest_gene"] = {
        str(k): int(v) for k, v in summary["relative_to_nearest_gene"].items()}
    located.head(5000).to_csv(OUT / "intergenic_located_sample.csv", index=False)
    (OUT / "intergenic_orientation_audit.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main(argv: list[str] | None = None):
    argv = list(sys.argv[1:] if argv is None else argv)
    only_orient = "--orient-only" in argv
    skip_reca = only_orient or "--skip-reca" in argv or (OUT / "reca_flank_summary.csv").exists()
    if not only_orient:
        m = load_joins()
        diag = dict(
            n_guides=len(m),
            n_s4_cas9=int(m["s4_cas9"].notna().sum()),
            n_s4_reca=int(m["s4_reca"].notna().sum()),
            spearman_cut_vs_s4=float(
                spearmanr(m["cut_score"], m["s4_cas9"], nan_policy="omit").statistic),
            spearman_cut_vs_s6abs=float(
                spearmanr(m["cut_score"], m["s6_cas9"].abs(), nan_policy="omit").statistic),
            spearman_s4_cas9_vs_reca=float(
                spearmanr(m["s4_cas9"], m["s4_reca"], nan_policy="omit").statistic),
            no_per_replicate_columns=True,
            published_data_files="Data S1–S7 + Supplementary Material PDF",
            replicate_evidence=(
                "PDF Tables S4/S5 list R1/R2 library sizes; Fig S12 shows replicate "
                "agreement. Per-guide replicate values were averaged before release."
            ),
        )
        (OUT / "guo_label_join.json").write_text(json.dumps(diag, indent=2))
        print(json.dumps(diag, indent=2), flush=True)
        if skip_reca:
            print("skipping ΔrecA CV (existing results or --skip-reca)", flush=True)
        else:
            run_reca(m)
    orientation_audit()
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
