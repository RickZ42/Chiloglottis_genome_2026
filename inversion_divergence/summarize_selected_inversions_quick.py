#!/usr/bin/env python3
"""Quick divergence summary for selected SyRI inversions (H1 genes vs colinear background).

This script reuses the previously computed H1-H2 one-to-one anchor divergence table and
the H1 inversion gene-overlap table, so it does not recompute alignments/dN/dS.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import pandas as pd

try:
    from scipy.stats import mannwhitneyu
except Exception:  # pragma: no cover
    mannwhitneyu = None


DEFAULT_OVERLAPS = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/H1/breaker/H1_INV_gene_overlaps.tsv"
)
DEFAULT_DIVERGENCE = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/H1H2_inversion_divergence/"
    "H1H2_anchor_unique_pairs.inversion_vs_colinear.divergence.tsv"
)
DEFAULT_OUTDIR = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/H1H2_inversion_divergence/selected_inversions_quick"
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


def finite_numeric(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    return s[s.map(lambda x: pd.notna(x) and math.isfinite(float(x)))]


def compare_metrics(fg: pd.DataFrame, bg: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for metric in METRICS:
        x = finite_numeric(fg.get(metric, pd.Series(dtype=float)))
        y = finite_numeric(bg.get(metric, pd.Series(dtype=float)))
        if len(x) == 0 or len(y) == 0:
            rows.append(
                {
                    "comparison": label,
                    "metric": metric,
                    "fg_n": len(x),
                    "fg_median": math.nan,
                    "fg_mean": math.nan,
                    "bg_n": len(y),
                    "bg_median": math.nan,
                    "bg_mean": math.nan,
                    "mannwhitney_u": math.nan,
                    "pvalue_two_sided": math.nan,
                }
            )
            continue
        if mannwhitneyu is None:
            u = p = math.nan
        else:
            try:
                res = mannwhitneyu(x, y, alternative="two-sided")
                u, p = float(res.statistic), float(res.pvalue)
            except Exception:
                u = p = math.nan
        rows.append(
            {
                "comparison": label,
                "metric": metric,
                "fg_n": int(len(x)),
                "fg_median": float(x.median()),
                "fg_mean": float(x.mean()),
                "bg_n": int(len(y)),
                "bg_median": float(y.median()),
                "bg_mean": float(y.mean()),
                "mannwhitney_u": u,
                "pvalue_two_sided": p,
            }
        )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--inv-ids",
        default="INV1936,INV1915",
        help="Comma-separated inversion IDs (default: INV1936,INV1915)",
    )
    p.add_argument("--overlaps", type=Path, default=DEFAULT_OVERLAPS)
    p.add_argument("--divergence", type=Path, default=DEFAULT_DIVERGENCE)
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    inv_ids = [x.strip() for x in args.inv_ids.split(",") if x.strip()]
    if not inv_ids:
        raise SystemExit("No inversion IDs provided")

    outdir = args.outdir / ("_".join(inv_ids))
    outdir.mkdir(parents=True, exist_ok=True)

    overlaps = pd.read_csv(args.overlaps, sep="\t", dtype={"INV_ID": str, "GeneID": str})
    div = pd.read_csv(args.divergence, sep="\t")

    # Keep only valid one-to-one anchor divergence rows (same table generated from unique pairs).
    div_ok = div[div["status"] == "ok"].copy()
    div_ok["H1_root"] = div_ok["H1_root"].astype(str)

    # Background = all non-inversion genes in the same anchor divergence universe.
    background = div_ok[~div_ok["is_inversion_gene"].astype(bool)].copy()

    selected_overlaps = overlaps[overlaps["INV_ID"].isin(inv_ids)].copy()
    selected_overlaps = selected_overlaps.drop_duplicates(subset=["INV_ID", "GeneID"])
    selected_gene_set = set(selected_overlaps["GeneID"].astype(str))

    selected_pairs = div_ok[div_ok["H1_root"].isin(selected_gene_set)].copy()
    selected_pairs = selected_pairs.merge(
        selected_overlaps[["INV_ID", "GeneID", "RefChr", "RefStart", "RefEnd", "QryChr", "QryStart", "QryEnd"]]
        .drop_duplicates(),
        left_on="H1_root",
        right_on="GeneID",
        how="left",
        validate="many_to_one",
    )

    # Save gene list with anchor presence/absence.
    anchor_roots = set(div_ok["H1_root"].astype(str))
    selected_overlaps["in_anchor_divergence_table"] = selected_overlaps["GeneID"].astype(str).isin(anchor_roots)
    selected_overlaps.to_csv(outdir / "selected_inversion_overlap_genes.with_anchor_presence.tsv", sep="\t", index=False)

    # Save selected pair table.
    cols_front = [
        "INV_ID",
        "RefChr",
        "RefStart",
        "RefEnd",
        "QryChr",
        "QryStart",
        "QryEnd",
        "H1_gene",
        "H1_root",
        "H2_gene",
        "H2_root",
    ]
    ordered_cols = cols_front + [c for c in selected_pairs.columns if c not in cols_front]
    selected_pairs[ordered_cols].sort_values(["INV_ID", "H1_root"]).to_csv(
        outdir / "selected_inversions.anchor_divergence_pairs.tsv", sep="\t", index=False
    )

    # Metric comparisons: combined and per inversion.
    metric_tables = []
    combined = compare_metrics(selected_pairs, background, label=f"{'+'.join(inv_ids)}_vs_colinear")
    metric_tables.append(combined)
    for inv_id in inv_ids:
        inv_df = selected_pairs[selected_pairs["INV_ID"] == inv_id].copy()
        metric_tables.append(compare_metrics(inv_df, background, label=f"{inv_id}_vs_colinear"))
    metric_summary = pd.concat(metric_tables, ignore_index=True)
    metric_summary.to_csv(outdir / "selected_inversions_vs_colinear.metric_summary.tsv", sep="\t", index=False)

    # Compact counts table per inversion.
    counts_rows = []
    for inv_id in inv_ids:
        ov = selected_overlaps[selected_overlaps["INV_ID"] == inv_id]
        sp = selected_pairs[selected_pairs["INV_ID"] == inv_id]
        missing = sorted(set(ov["GeneID"].astype(str)) - set(sp["H1_root"].astype(str)))
        counts_rows.append(
            {
                "INV_ID": inv_id,
                "overlap_genes_total": int(ov.shape[0]),
                "anchor_paired_genes_used": int(sp["H1_root"].nunique()),
                "missing_from_anchor_pair_set": int(len(missing)),
                "missing_geneIDs": ",".join(missing),
            }
        )
    counts_df = pd.DataFrame(counts_rows)
    counts_df.to_csv(outdir / "selected_inversions.gene_counts.tsv", sep="\t", index=False)

    # Text summary.
    lines = []
    lines.append("Quick divergence summary for selected inversions (reusing H1-H2 anchor divergence table)")
    lines.append("=" * 78)
    lines.append(f"Selected inversion IDs: {', '.join(inv_ids)}")
    lines.append(f"Overlaps table: {args.overlaps}")
    lines.append(f"Divergence table (precomputed one-to-one anchors): {args.divergence}")
    lines.append("")
    lines.append(f"Background definition: status=ok AND is_inversion_gene=False in divergence table")
    lines.append(f"Background genes (anchor pairs): {background.shape[0]}")
    lines.append("")
    lines.append("Per-inversion gene availability:")
    for _, r in counts_df.iterrows():
        lines.append(
            f"  {r['INV_ID']}: overlap_genes_total={r['overlap_genes_total']}, "
            f"anchor_paired_genes_used={r['anchor_paired_genes_used']}, "
            f"missing_from_anchor_pair_set={r['missing_from_anchor_pair_set']}"
        )
        if r["missing_geneIDs"]:
            lines.append(f"    missing_geneIDs: {r['missing_geneIDs']}")
    lines.append("")
    lines.append(f"Combined selected inversion genes used (unique H1 roots): {selected_pairs['H1_root'].nunique()}")
    lines.append("")

    # Add a compact metric summary for dS/dN/dNdS and distances.
    for comp in metric_summary["comparison"].drop_duplicates():
        lines.append(f"[{comp}]")
        sub = metric_summary[metric_summary["comparison"] == comp]
        for metric in ["dS_NG86", "dN_NG86", "dN_dS_NG86", "cds_nt_pdistance", "aa_pdistance"]:
            r = sub[sub["metric"] == metric]
            if r.empty:
                continue
            r = r.iloc[0]
            lines.append(
                f"  {metric}: fg n={int(r['fg_n'])}, median={r['fg_median']:.6g}; "
                f"bg n={int(r['bg_n'])}, median={r['bg_median']:.6g}; "
                f"p={r['pvalue_two_sided']:.6g}" if pd.notna(r["pvalue_two_sided"]) else
                f"  {metric}: fg n={int(r['fg_n'])}, median={r['fg_median']:.6g}; "
                f"bg n={int(r['bg_n'])}, median={r['bg_median']:.6g}; p=NA"
            )
        lines.append("")

    summary_path = outdir / "selected_inversions_quick.divergence_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")

    print(f"[INFO] wrote: {outdir}")
    for fn in [
        "selected_inversion_overlap_genes.with_anchor_presence.tsv",
        "selected_inversions.anchor_divergence_pairs.tsv",
        "selected_inversions.gene_counts.tsv",
        "selected_inversions_vs_colinear.metric_summary.tsv",
        "selected_inversions_quick.divergence_summary.txt",
    ]:
        print(f"[INFO]   {outdir / fn}")


if __name__ == "__main__":
    main()

