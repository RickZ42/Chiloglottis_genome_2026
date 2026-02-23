#!/usr/bin/env python3
"""High-confidence inversion interior-gene divergence vs matched colinear controls.

Uses precomputed H1-H2 one-to-one anchor divergence metrics and matches each foreground
gene to a non-inversion colinear control on the same chromosome with closest CDS length.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

try:
    from scipy.stats import mannwhitneyu, wilcoxon
except Exception:  # pragma: no cover
    mannwhitneyu = None
    wilcoxon = None


DEFAULT_HC = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/syri/syri_asm10/H1_vs_H2syri.highconfINV10kb_coordinates.tsv"
)
DEFAULT_OVERLAPS = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/H1/breaker/H1_INV_gene_overlaps.tsv"
)
DEFAULT_DIV = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/H1H2_inversion_divergence/H1H2_anchor_unique_pairs.inversion_vs_colinear.divergence.tsv"
)
DEFAULT_H1BED = Path(
    "/g/data/xf3/zz3507/compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H1t1.bed"
)
DEFAULT_OUTDIR = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/H1H2_inversion_divergence/highconf_interior_matched"
)

METRICS = [
    "dS_NG86",
    "dN_NG86",
    "dN_dS_NG86",
    "cds_nt_pdistance",
    "aa_pdistance",
    "cds_nt_identity",
    "aa_identity",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--highconf", type=Path, default=DEFAULT_HC)
    p.add_argument("--overlaps", type=Path, default=DEFAULT_OVERLAPS)
    p.add_argument("--divergence", type=Path, default=DEFAULT_DIV)
    p.add_argument("--h1-bed", type=Path, default=DEFAULT_H1BED)
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument(
        "--priority-ids",
        default="INV1936,INV1915",
        help="Comma-separated inversion IDs to report separately",
    )
    return p.parse_args()


def finite_series(df: pd.DataFrame, col: str) -> pd.Series:
    s = pd.to_numeric(df[col], errors="coerce")
    return s[s.map(lambda x: pd.notna(x) and math.isfinite(float(x)))]


def summary_rows(
    fg: pd.DataFrame, bg: pd.DataFrame, label: str, test_prefix: str = "mw"
) -> List[dict]:
    rows = []
    for m in METRICS:
        x = finite_series(fg, m)
        y = finite_series(bg, m)
        if len(x) == 0 or len(y) == 0:
            rows.append(
                dict(
                    comparison=label,
                    metric=m,
                    fg_n=int(len(x)),
                    fg_median=math.nan,
                    fg_mean=math.nan,
                    bg_n=int(len(y)),
                    bg_median=math.nan,
                    bg_mean=math.nan,
                    test=f"{test_prefix}_two_sided",
                    statistic=math.nan,
                    pvalue=math.nan,
                )
            )
            continue
        stat = p = math.nan
        if mannwhitneyu is not None:
            try:
                res = mannwhitneyu(x, y, alternative="two-sided")
                stat, p = float(res.statistic), float(res.pvalue)
            except Exception:
                pass
        rows.append(
            dict(
                comparison=label,
                metric=m,
                fg_n=int(len(x)),
                fg_median=float(x.median()),
                fg_mean=float(x.mean()),
                bg_n=int(len(y)),
                bg_median=float(y.median()),
                bg_mean=float(y.mean()),
                test=f"{test_prefix}_two_sided",
                statistic=stat,
                pvalue=p,
            )
        )
    return rows


def paired_summary_rows(pairs: pd.DataFrame, label: str) -> List[dict]:
    rows = []
    for m in METRICS:
        x = pd.to_numeric(pairs[f"fg_{m}"], errors="coerce")
        y = pd.to_numeric(pairs[f"ctrl_{m}"], errors="coerce")
        valid = x.notna() & y.notna() & x.map(math.isfinite) & y.map(math.isfinite)
        x = x[valid]
        y = y[valid]
        if len(x) == 0:
            rows.append(
                dict(
                    comparison=label,
                    metric=m,
                    pair_n=0,
                    median_delta_fg_minus_ctrl=math.nan,
                    mean_delta_fg_minus_ctrl=math.nan,
                    test="wilcoxon_two_sided",
                    statistic=math.nan,
                    pvalue=math.nan,
                )
            )
            continue
        d = x - y
        stat = p = math.nan
        if wilcoxon is not None:
            try:
                # zero_method='wilcox' drops zero differences
                res = wilcoxon(d, alternative="two-sided", zero_method="wilcox")
                stat, p = float(res.statistic), float(res.pvalue)
            except Exception:
                pass
        rows.append(
            dict(
                comparison=label,
                metric=m,
                pair_n=int(len(d)),
                median_delta_fg_minus_ctrl=float(d.median()),
                mean_delta_fg_minus_ctrl=float(d.mean()),
                test="wilcoxon_two_sided",
                statistic=stat,
                pvalue=p,
            )
        )
    return rows


def greedy_match_same_chr_by_cds_len(
    fg: pd.DataFrame, ctrl_pool: pd.DataFrame
) -> pd.DataFrame:
    """Deterministic 1:1 matching without replacement.

    Matching key priority:
      1) same chromosome
      2) minimum |log2((ctrl_cds+1)/(fg_cds+1))|
      3) minimum absolute CDS length difference
      4) lexical H1_gene
    """

    remaining: Dict[str, pd.DataFrame] = {
        chr_: sub.sort_values(["H1_gene"]).copy()
        for chr_, sub in ctrl_pool.groupby("H1_chr", dropna=False)
    }

    match_rows = []
    fg_order = fg.copy().sort_values(["H1_chr", "H1_cds_len_bp", "H1_gene"], ascending=[True, False, True])

    for _, r in fg_order.iterrows():
        chr_ = r["H1_chr"]
        cand = remaining.get(chr_, pd.DataFrame()).copy()
        if cand.empty:
            # fallback to all remaining controls if same chromosome depleted
            cand = pd.concat([x for x in remaining.values() if not x.empty], ignore_index=True)
            same_chr = False
        else:
            same_chr = True
        if cand.empty:
            continue
        fg_len = float(r["H1_cds_len_bp"])
        cand["len_absdiff"] = (pd.to_numeric(cand["H1_cds_len_bp"], errors="coerce") - fg_len).abs()
        cand["len_log2diff"] = (
            (pd.to_numeric(cand["H1_cds_len_bp"], errors="coerce") + 1.0) / (fg_len + 1.0)
        ).map(lambda x: abs(math.log2(x)) if x > 0 else math.inf)
        cand = cand.sort_values(["len_log2diff", "len_absdiff", "H1_gene"])
        best = cand.iloc[0]

        # remove matched control from its source chromosome bucket
        src_chr = best["H1_chr"]
        src_bucket = remaining.get(src_chr, pd.DataFrame()).copy()
        if not src_bucket.empty:
            remaining[src_chr] = src_bucket[src_bucket["H1_gene"] != best["H1_gene"]].copy()

        row = {}
        for c in fg.columns:
            row[f"fg_{c}"] = r[c]
        for c in ctrl_pool.columns:
            row[f"ctrl_{c}"] = best[c]
        row["matched_same_chr"] = bool(same_chr and (best["H1_chr"] == chr_))
        row["match_len_absdiff"] = float(best["len_absdiff"])
        row["match_len_log2diff"] = float(best["len_log2diff"])
        match_rows.append(row)

    return pd.DataFrame(match_rows)


def main() -> None:
    args = parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    priority_ids = [x.strip() for x in args.priority_ids.split(",") if x.strip()]

    hc = pd.read_csv(args.highconf, sep="\t")
    hc = hc.rename(columns={"ID": "INV_ID"})
    hc_ids = set(hc["INV_ID"].astype(str))

    overlaps = pd.read_csv(args.overlaps, sep="\t")
    overlaps = overlaps[overlaps["INV_ID"].isin(hc_ids)].copy()
    overlaps["is_interior"] = (overlaps["GeneStart"] >= overlaps["RefStart"]) & (overlaps["GeneEnd"] <= overlaps["RefEnd"])
    interior = overlaps[overlaps["is_interior"]].copy()
    interior = interior.drop_duplicates(subset=["INV_ID", "GeneID"])

    div = pd.read_csv(args.divergence, sep="\t")
    div = div[div["status"] == "ok"].copy()

    # H1 t1 coordinates for chromosome-aware matching
    bed_cols = ["H1_chr", "H1_start0", "H1_end", "H1_gene", "score", "strand"]
    bed = pd.read_csv(args.h1_bed, sep="\t", header=None, names=bed_cols)
    bed = bed[bed["H1_gene"].astype(str).str.endswith(".t1")].copy()
    bed["H1_root"] = bed["H1_gene"].astype(str).str.replace(r"\.t\d+$", "", regex=True)
    bed = bed.drop_duplicates(subset=["H1_gene"])

    div = div.merge(
        bed[["H1_gene", "H1_chr", "H1_start0", "H1_end", "strand"]],
        on="H1_gene",
        how="left",
        validate="one_to_one",
    )

    interior_gene_set = set(interior["GeneID"].astype(str))
    div["is_highconf_interior_gene"] = div["H1_root"].astype(str).isin(interior_gene_set)

    fg_all = div[div["is_highconf_interior_gene"]].copy()
    ctrl_pool = div[~div["is_inversion_gene"].astype(bool)].copy()

    # Attach inversion IDs to foreground rows
    fg_all = fg_all.merge(
        interior[["INV_ID", "GeneID", "RefChr", "RefStart", "RefEnd", "QryChr", "QryStart", "QryEnd"]].rename(columns={"GeneID": "H1_root"}),
        on="H1_root",
        how="left",
        validate="many_to_one",
    )

    # Counts and missing coverage by inversion
    anchor_roots = set(div["H1_root"].astype(str))
    interior["in_anchor_divergence_table"] = interior["GeneID"].astype(str).isin(anchor_roots)
    counts = []
    for inv_id, sub in interior.groupby("INV_ID"):
        missing = sorted(set(sub["GeneID"].astype(str)) - anchor_roots)
        used = sorted(set(sub["GeneID"].astype(str)) & anchor_roots)
        counts.append(
            {
                "INV_ID": inv_id,
                "interior_genes_total": int(sub["GeneID"].nunique()),
                "anchor_paired_genes_used": int(len(used)),
                "missing_from_anchor_pair_set": int(len(missing)),
                "used_genes": ",".join(used),
                "missing_genes": ",".join(missing),
            }
        )
    counts_df = pd.DataFrame(counts).sort_values("INV_ID")

    # 1:1 matched controls for all foreground genes
    matched_pairs = greedy_match_same_chr_by_cds_len(fg_all, ctrl_pool)

    # Save tables
    interior.to_csv(outdir / "highconf10kb_interior_genes.anchor_presence.tsv", sep="\t", index=False)
    counts_df.to_csv(outdir / "highconf10kb_interior_gene_counts.tsv", sep="\t", index=False)
    fg_all.sort_values(["INV_ID", "H1_root"]).to_csv(outdir / "highconf10kb_interior.anchor_divergence_fg.tsv", sep="\t", index=False)
    matched_pairs.to_csv(outdir / "highconf10kb_interior_vs_matched_colinear.pairs.tsv", sep="\t", index=False)

    # Summaries
    summary_rows_all = []
    summary_rows_all.extend(summary_rows(fg_all, ctrl_pool, label="highconf10kb_interior_vs_all_colinear"))
    if not matched_pairs.empty:
        fg_match = matched_pairs[[c for c in matched_pairs.columns if c.startswith("fg_")]].copy()
        fg_match.columns = [c[3:] for c in fg_match.columns]
        ctrl_match = matched_pairs[[c for c in matched_pairs.columns if c.startswith("ctrl_")]].copy()
        ctrl_match.columns = [c[5:] for c in ctrl_match.columns]
        summary_rows_all.extend(summary_rows(fg_match, ctrl_match, label="highconf10kb_interior_vs_matched_colinear", test_prefix="mw"))
        paired_rows = paired_summary_rows(matched_pairs, label="highconf10kb_interior_vs_matched_colinear_paired")
    else:
        paired_rows = []

    # Priority inversions (vs all colinear and matched subset if possible)
    for inv_id in priority_ids:
        fg_sub = fg_all[fg_all["INV_ID"] == inv_id].copy()
        if fg_sub.empty:
            continue
        summary_rows_all.extend(summary_rows(fg_sub, ctrl_pool, label=f"{inv_id}_interior_vs_all_colinear"))
        mp_sub = matched_pairs[matched_pairs["fg_INV_ID"] == inv_id].copy() if not matched_pairs.empty else pd.DataFrame()
        if not mp_sub.empty:
            fg_s = mp_sub[[c for c in mp_sub.columns if c.startswith("fg_")]].copy()
            fg_s.columns = [c[3:] for c in fg_s.columns]
            ctrl_s = mp_sub[[c for c in mp_sub.columns if c.startswith("ctrl_")]].copy()
            ctrl_s.columns = [c[5:] for c in ctrl_s.columns]
            summary_rows_all.extend(summary_rows(fg_s, ctrl_s, label=f"{inv_id}_interior_vs_matched_colinear"))
            paired_rows.extend(paired_summary_rows(mp_sub, label=f"{inv_id}_interior_vs_matched_colinear_paired"))

    metric_summary = pd.DataFrame(summary_rows_all)
    metric_summary.to_csv(outdir / "highconf10kb_interior.metric_summary.tsv", sep="\t", index=False)
    paired_summary = pd.DataFrame(paired_rows)
    paired_summary.to_csv(outdir / "highconf10kb_interior.matched_pairs_paired_tests.tsv", sep="\t", index=False)

    # Text summary
    lines = []
    lines.append("High-confidence inversion (10 kb set) interior-gene divergence summary")
    lines.append("=" * 74)
    lines.append(f"High-confidence set: {args.highconf}")
    lines.append("Interior gene definition: GeneStart >= RefStart and GeneEnd <= RefEnd in H1 overlap table")
    lines.append("Foreground analysis universe: precomputed H1-H2 one-to-one anchor divergence rows with status=ok")
    lines.append("Matched control definition: non-inversion colinear genes (is_inversion_gene=False), 1:1 greedy same-chromosome nearest CDS length match")
    lines.append("")
    lines.append(f"High-confidence inversions in set: {', '.join(sorted(hc_ids))}")
    lines.append(f"High-confidence interior genes total (H1 overlaps): {interior['GeneID'].nunique()}")
    lines.append(f"High-confidence interior genes in anchor divergence table: {fg_all['H1_root'].nunique()}")
    lines.append(f"Missing (no anchor divergence row): {interior['GeneID'].nunique() - fg_all['H1_root'].nunique()}")
    lines.append(f"All-colinear background genes: {ctrl_pool.shape[0]}")
    lines.append(f"Matched control pairs built: {matched_pairs.shape[0]}")
    if not matched_pairs.empty:
        same_chr_n = int(pd.to_numeric(matched_pairs['matched_same_chr'], errors='coerce').fillna(False).sum())
        lines.append(f"Matched pairs on same chromosome: {same_chr_n}/{matched_pairs.shape[0]}")
        lines.append(
            f"Matched CDS length absdiff median: {matched_pairs['match_len_absdiff'].median():.1f} bp "
            f"(mean {matched_pairs['match_len_absdiff'].mean():.1f})"
        )
        lines.append(
            f"Matched CDS length |log2 ratio| median: {matched_pairs['match_len_log2diff'].median():.4f}"
        )
    lines.append("")
    lines.append("Per-inversion anchor coverage (interior genes):")
    for _, r in counts_df.iterrows():
        lines.append(
            f"  {r['INV_ID']}: total={r['interior_genes_total']}, used={r['anchor_paired_genes_used']}, missing={r['missing_from_anchor_pair_set']}"
        )
    lines.append("")
    for label in [
        "highconf10kb_interior_vs_all_colinear",
        "highconf10kb_interior_vs_matched_colinear",
        "highconf10kb_interior_vs_matched_colinear_paired",
    ] + [f"{x}_interior_vs_all_colinear" for x in priority_ids] + [f"{x}_interior_vs_matched_colinear" for x in priority_ids]:
        if label.endswith("_paired"):
            sub = paired_summary[paired_summary["comparison"] == label] if not paired_summary.empty else pd.DataFrame()
            if sub.empty:
                continue
            lines.append(f"[{label}]")
            for metric in ["dS_NG86", "dN_NG86", "dN_dS_NG86", "cds_nt_pdistance", "aa_pdistance"]:
                rr = sub[sub["metric"] == metric]
                if rr.empty:
                    continue
                rr = rr.iloc[0]
                lines.append(
                    f"  {metric}: pairs={int(rr['pair_n'])}, median_delta(fg-ctrl)={rr['median_delta_fg_minus_ctrl']:.6g}, "
                    f"p={rr['pvalue']:.6g}" if pd.notna(rr["pvalue"]) else
                    f"  {metric}: pairs={int(rr['pair_n'])}, median_delta(fg-ctrl)={rr['median_delta_fg_minus_ctrl']:.6g}, p=NA"
                )
            lines.append("")
        else:
            sub = metric_summary[metric_summary["comparison"] == label]
            if sub.empty:
                continue
            lines.append(f"[{label}]")
            for metric in ["dS_NG86", "dN_NG86", "dN_dS_NG86", "cds_nt_pdistance", "aa_pdistance"]:
                rr = sub[sub["metric"] == metric]
                if rr.empty:
                    continue
                rr = rr.iloc[0]
                lines.append(
                    f"  {metric}: fg n={int(rr['fg_n'])}, median={rr['fg_median']:.6g}; bg n={int(rr['bg_n'])}, "
                    f"median={rr['bg_median']:.6g}; p={rr['pvalue']:.6g}" if pd.notna(rr["pvalue"]) else
                    f"  {metric}: fg n={int(rr['fg_n'])}, median={rr['fg_median']:.6g}; bg n={int(rr['bg_n'])}, median={rr['bg_median']:.6g}; p=NA"
                )
            lines.append("")

    (outdir / "highconf10kb_interior.divergence_summary.txt").write_text("\n".join(lines) + "\n")

    print(f"[INFO] wrote output directory: {outdir}")
    for fn in [
        "highconf10kb_interior.divergence_summary.txt",
        "highconf10kb_interior_gene_counts.tsv",
        "highconf10kb_interior.anchor_divergence_fg.tsv",
        "highconf10kb_interior_vs_matched_colinear.pairs.tsv",
        "highconf10kb_interior.metric_summary.tsv",
        "highconf10kb_interior.matched_pairs_paired_tests.tsv",
    ]:
        print(f"[INFO]   {outdir / fn}")


if __name__ == "__main__":
    main()

