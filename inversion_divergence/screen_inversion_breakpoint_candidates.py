#!/usr/bin/env python3
"""Prioritize inversion breakpoint cases for follow-up structural/RNA inspection.

This is a coarse haplotype-side screen intended to shortlist candidate inversions for:
  - breakpoint-overlapping genes
  - promoter-proximal breakpoints
  - exon / intron-boundary-adjacent breakpoints
  - expressed nearby genes worth manual IGV follow-up

It does not attempt formal differential-expression inference. When expression tables
are supplied, it combines genomic breakpoint context with all-library RNA abundance
summaries so the top ranked candidates can be reviewed manually in IGV/H1-H2
follow-up.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


DEFAULT_INV = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/syri/syri_asm10/H1_vs_H2syri.highconfINV.tsv"
)
DEFAULT_T1_BED = Path(
    "/g/data/xf3/zz3507/compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H1t1.bed"
)
DEFAULT_GTF = Path("/g/data/xf3/zz3507/Output/20260127Genome/H1/breaker/braker.gtf")
DEFAULT_GENE_COUNTS = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/H1/inversion_rna_minimap2/rna_gene_counts.tsv"
)
DEFAULT_INV_OVERVIEW = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/H1/inversion_rna_minimap2/"
    "inv_rna_expr_minimap2_all.overview.tsv"
)
DEFAULT_HIC = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/hic_inv_context_joint/hic_inv_context.summary.tsv"
)
DEFAULT_OUTDIR = Path(
    "/g/data/xf3/zz3507/Output/20260127Genome/other_syri_inversion/priority_breakpoint_screen"
)


@dataclass(frozen=True)
class Transcript:
    gene_id: str
    transcript_id: str
    chrom: str
    start: int
    end: int
    strand: str
    promoter_start: int
    promoter_end: int
    mean_expr: float
    max_expr: float


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inversions", type=Path, default=DEFAULT_INV)
    p.add_argument("--t1-bed", "--h1-t1-bed", dest="t1_bed", type=Path, default=DEFAULT_T1_BED)
    p.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    p.add_argument("--gene-counts", type=Path, default=DEFAULT_GENE_COUNTS)
    p.add_argument("--inv-overview", type=Path, default=DEFAULT_INV_OVERVIEW)
    p.add_argument("--hic-summary", type=Path, default=DEFAULT_HIC)
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument("--coord-prefix", choices=["Ref", "Qry"], default="Ref")
    p.add_argument("--label", default="H1", help="Short label for the screened haplotype/reference.")
    p.add_argument("--promoter-bp", type=int, default=2000)
    p.add_argument("--flank-bp", type=int, default=10000)
    p.add_argument("--boundary-bp", type=int, default=200)
    p.add_argument("--top-n", type=int, default=10)
    p.add_argument(
        "--skip-expression",
        action="store_true",
        help="Ignore gene-count / overview tables and run a structure-only breakpoint screen.",
    )
    p.add_argument(
        "--strict-only",
        action="store_true",
        help="Restrict the screen to rows with HighConfidence=True in the input table.",
    )
    return p.parse_args()


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return False
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def parse_gene_counts(path: Path) -> Dict[str, Tuple[float, float]]:
    if not path or not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", comment="#")
    sample_cols = [c for c in df.columns if c not in {"Geneid", "Chr", "Start", "End", "Strand", "Length"}]
    counts = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    means = counts.mean(axis=1)
    maxima = counts.max(axis=1)
    return {
        str(gene_id): (float(mean_val), float(max_val))
        for gene_id, mean_val, max_val in zip(df["Geneid"].astype(str), means, maxima)
    }


def load_transcripts(
    bed_path: Path,
    expr_by_gene: Dict[str, Tuple[float, float]],
    promoter_bp: int,
) -> Tuple[Dict[str, Transcript], Dict[str, List[Transcript]]]:
    tx_by_gene: Dict[str, Transcript] = {}
    tx_by_chr: Dict[str, List[Transcript]] = {}
    with bed_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, start0, end0, tx_id, _, strand = line.rstrip("\n").split("\t")[:6]
            if not tx_id.endswith(".t1"):
                continue
            gene_id = tx_id.rsplit(".", 1)[0]
            start = int(start0) + 1
            end = int(end0)
            mean_expr, max_expr = expr_by_gene.get(gene_id, (math.nan, math.nan))
            if strand == "+":
                promoter_start = max(1, start - promoter_bp)
                promoter_end = max(1, start - 1)
            else:
                promoter_start = end + 1
                promoter_end = end + promoter_bp
            tx = Transcript(
                gene_id=gene_id,
                transcript_id=tx_id,
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                promoter_start=promoter_start,
                promoter_end=promoter_end,
                mean_expr=mean_expr,
                max_expr=max_expr,
            )
            tx_by_gene[gene_id] = tx
            tx_by_chr.setdefault(chrom, []).append(tx)
    for chrom in tx_by_chr:
        tx_by_chr[chrom].sort(key=lambda x: (x.start, x.end, x.gene_id))
    return tx_by_gene, tx_by_chr


def parse_gtf_attributes(attr_field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_field.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " in part:
            key, value = part.split(" ", 1)
            attrs[key] = value.strip().strip('"')
        else:
            attrs[part] = ""
    return attrs


def load_feature_intervals(
    gtf_path: Path,
) -> Tuple[Dict[str, List[Tuple[int, int]]], Dict[str, List[Tuple[int, int]]]]:
    exons: Dict[str, List[Tuple[int, int]]] = {}
    introns: Dict[str, List[Tuple[int, int]]] = {}
    with gtf_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            feature = fields[2]
            if feature not in {"exon", "intron"}:
                continue
            attrs = parse_gtf_attributes(fields[8])
            tx_id = attrs.get("transcript_id", "")
            gene_id = attrs.get("gene_id", "")
            if not tx_id.endswith(".t1") or not gene_id:
                continue
            start = int(fields[3])
            end = int(fields[4])
            if feature == "exon":
                exons.setdefault(gene_id, []).append((start, end))
            else:
                introns.setdefault(gene_id, []).append((start, end))
    for d in (exons, introns):
        for gene_id in d:
            d[gene_id].sort()
    return exons, introns


def interval_contains(interval: Tuple[int, int], pos: int) -> bool:
    return interval[0] <= pos <= interval[1]


def min_distance_to_boundaries(exons: Sequence[Tuple[int, int]], pos: int) -> Optional[int]:
    if not exons:
        return None
    return min(min(abs(pos - start), abs(pos - end)) for start, end in exons)


def distance_to_gene(start: int, end: int, pos: int) -> int:
    if start <= pos <= end:
        return 0
    if pos < start:
        return start - pos
    return pos - end


def classify_gene_body_relation(
    gene_id: str,
    pos: int,
    exon_map: Dict[str, List[Tuple[int, int]]],
    intron_map: Dict[str, List[Tuple[int, int]]],
    boundary_bp: int,
) -> Tuple[str, Optional[int], bool]:
    exons = exon_map.get(gene_id, [])
    introns = intron_map.get(gene_id, [])
    in_exon = any(interval_contains(iv, pos) for iv in exons)
    in_intron = any(interval_contains(iv, pos) for iv in introns)
    boundary_dist = min_distance_to_boundaries(exons, pos)
    boundary_near = boundary_dist is not None and boundary_dist <= boundary_bp
    if in_exon:
        relation = "exon"
    elif in_intron:
        relation = "intron"
    else:
        relation = "gene_body_other"
    return relation, boundary_dist, boundary_near


def format_expr(value: float) -> str:
    if pd.isna(value):
        return "NA"
    return f"{value:.1f}"


def choose_best_event(events: Sequence[dict]) -> Optional[dict]:
    if not events:
        return None

    def event_rank(event: dict) -> Tuple[int, float]:
        relation = event["relation"]
        base = {
            "exon": 5,
            "promoter": 4,
            "intron": 2,
            "flank": 1,
            "gene_body_other": 1,
        }.get(relation, 0)
        if event.get("boundary_near"):
            base += 2
        return (base, float(event.get("mean_expr", float("-inf"))))

    return max(events, key=event_rank)


def event_short_label(event: dict) -> str:
    label = f"{event['gene_id']}:{event['relation']}"
    if event.get("boundary_near"):
        label += f":bdry{event['boundary_dist']}bp"
    if not pd.isna(event.get("mean_expr", math.nan)):
        label += f":expr{event['mean_expr']:.1f}"
    return label


def analyze_breakpoint(
    chrom: str,
    pos: int,
    side: str,
    transcripts: Sequence[Transcript],
    exon_map: Dict[str, List[Tuple[int, int]]],
    intron_map: Dict[str, List[Tuple[int, int]]],
    boundary_bp: int,
    flank_bp: int,
) -> Tuple[dict, List[dict]]:
    events: List[dict] = []
    nearest_left: Optional[dict] = None
    nearest_right: Optional[dict] = None
    nearest_flank: Optional[dict] = None

    for tx in transcripts:
        gene_dist = distance_to_gene(tx.start, tx.end, pos)
        if tx.end < pos:
            left_row = {
                "gene_id": tx.gene_id,
                "distance_to_gene_bp": gene_dist,
                "mean_expr": tx.mean_expr,
            }
            if nearest_left is None or left_row["distance_to_gene_bp"] < nearest_left["distance_to_gene_bp"]:
                nearest_left = left_row
        elif tx.start > pos:
            right_row = {
                "gene_id": tx.gene_id,
                "distance_to_gene_bp": gene_dist,
                "mean_expr": tx.mean_expr,
            }
            if nearest_right is None or right_row["distance_to_gene_bp"] < nearest_right["distance_to_gene_bp"]:
                nearest_right = right_row

        if tx.start <= pos <= tx.end:
            relation, boundary_dist, boundary_near = classify_gene_body_relation(
                tx.gene_id, pos, exon_map, intron_map, boundary_bp
            )
            events.append(
                {
                    "boundary_side": side,
                    "breakpoint_chr": chrom,
                    "breakpoint_pos": pos,
                    "gene_id": tx.gene_id,
                    "relation": relation,
                    "distance_to_gene_bp": 0,
                    "boundary_dist": boundary_dist,
                    "boundary_near": boundary_near,
                    "mean_expr": tx.mean_expr,
                    "max_expr": tx.max_expr,
                    "strand": tx.strand,
                }
            )
        elif tx.promoter_start <= pos <= tx.promoter_end:
            promoter_dist = min(abs(pos - tx.start), abs(pos - tx.end))
            events.append(
                {
                    "boundary_side": side,
                    "breakpoint_chr": chrom,
                    "breakpoint_pos": pos,
                    "gene_id": tx.gene_id,
                    "relation": "promoter",
                    "distance_to_gene_bp": promoter_dist,
                    "boundary_dist": None,
                    "boundary_near": False,
                    "mean_expr": tx.mean_expr,
                    "max_expr": tx.max_expr,
                    "strand": tx.strand,
                }
            )
        elif gene_dist <= flank_bp:
            flank_event = {
                "boundary_side": side,
                "breakpoint_chr": chrom,
                "breakpoint_pos": pos,
                "gene_id": tx.gene_id,
                "relation": "flank",
                "distance_to_gene_bp": gene_dist,
                "boundary_dist": None,
                "boundary_near": False,
                "mean_expr": tx.mean_expr,
                "max_expr": tx.max_expr,
                "strand": tx.strand,
            }
            if nearest_flank is None or flank_event["distance_to_gene_bp"] < nearest_flank["distance_to_gene_bp"]:
                nearest_flank = flank_event

    if nearest_flank is not None and nearest_flank["gene_id"] not in {e["gene_id"] for e in events}:
        events.append(nearest_flank)

    overlap_events = [e for e in events if e["relation"] in {"exon", "intron", "gene_body_other"}]
    promoter_events = [e for e in events if e["relation"] == "promoter"]
    flank_events = [e for e in events if e["relation"] == "flank"]
    best = choose_best_event(events)

    summary = {
        "boundary_side": side,
        "breakpoint_chr": chrom,
        "breakpoint_pos": pos,
        "n_overlap_genes": len(overlap_events),
        "n_exon_hits": sum(e["relation"] == "exon" for e in overlap_events),
        "n_intron_hits": sum(e["relation"] == "intron" for e in overlap_events),
        "n_promoter_hits": len(promoter_events),
        "n_boundary_near_hits": sum(bool(e["boundary_near"]) for e in overlap_events),
        "best_event": event_short_label(best) if best else "NA",
        "nearest_left_gene": nearest_left["gene_id"] if nearest_left else "NA",
        "nearest_left_dist_bp": nearest_left["distance_to_gene_bp"] if nearest_left else math.nan,
        "nearest_left_mean_expr": nearest_left["mean_expr"] if nearest_left else math.nan,
        "nearest_right_gene": nearest_right["gene_id"] if nearest_right else "NA",
        "nearest_right_dist_bp": nearest_right["distance_to_gene_bp"] if nearest_right else math.nan,
        "nearest_right_mean_expr": nearest_right["mean_expr"] if nearest_right else math.nan,
        "has_expressed_boundary_gene": any(
            (not pd.isna(e["mean_expr"])) and e["mean_expr"] >= 100 for e in events if e["relation"] != "flank"
        ),
        "has_expressed_flank_gene": any(
            (not pd.isna(e["mean_expr"])) and e["mean_expr"] >= 100 and e["distance_to_gene_bp"] <= flank_bp
            for e in flank_events
        ),
        "max_boundary_related_mean_expr": max(
            [e["mean_expr"] for e in events if not pd.isna(e["mean_expr"])] or [math.nan]
        ),
    }
    return summary, events


def compute_priority_score(row: pd.Series) -> int:
    score = 0
    if row["strict_highconfidence"]:
        score += 1
    if row["has_hic_panel"]:
        score += 2
    if row["any_exon_hit"]:
        score += 5
    if row["any_promoter_hit"]:
        score += 4
    if row["any_boundary_near_hit"]:
        score += 4
    if row["any_intron_hit"]:
        score += 2
    if row["any_overlap_gene"]:
        score += 1
    if row["max_boundary_related_mean_expr"] >= 500:
        score += 2
    elif row["max_boundary_related_mean_expr"] >= 100:
        score += 1
    if row["nearest_gene_within_2kb"]:
        score += 1
    return score


def build_markdown_summary(
    inv_df: pd.DataFrame,
    top_df: pd.DataFrame,
    events_df: pd.DataFrame,
    args: argparse.Namespace,
) -> str:
    lines = []
    lines.append("# Inversion breakpoint priority screen")
    lines.append("")
    lines.append(f"- Screen label: `{args.label}`")
    lines.append(f"- Coordinate prefix used: `{args.coord_prefix}`")
    lines.append(f"- Input inversion table: `{args.inversions}`")
    lines.append(f"- Transcript BED: `{args.t1_bed}`")
    lines.append(f"- GTF: `{args.gtf}`")
    lines.append(f"- Total inversions screened: `{inv_df.shape[0]}`")
    lines.append(f"- Promoter window: `{args.promoter_bp}` bp")
    lines.append(f"- Flanking-expression window: `{args.flank_bp}` bp")
    lines.append(f"- Exon-boundary proximity threshold: `{args.boundary_bp}` bp")
    lines.append("")
    lines.append("## Category counts")
    lines.append("")
    for label, value in [
        ("Strict high-confidence inversions", int(inv_df["strict_highconfidence"].sum())),
        ("Any breakpoint-overlapping gene", int(inv_df["any_overlap_gene"].sum())),
        ("Any exon-overlapping breakpoint", int(inv_df["any_exon_hit"].sum())),
        ("Any promoter-overlapping breakpoint", int(inv_df["any_promoter_hit"].sum())),
        ("Any intron-overlapping breakpoint", int(inv_df["any_intron_hit"].sum())),
        ("Any exon-boundary-near breakpoint", int(inv_df["any_boundary_near_hit"].sum())),
        ("Any boundary-related gene with mean expression >= 100", int(inv_df["any_expressed_boundary_gene"].sum())),
    ]:
        lines.append(f"- {label}: `{value}`")
    lines.append("")
    lines.append(f"## Top {args.top_n} priority inversions")
    lines.append("")
    for _, row in top_df.iterrows():
        lines.append(
            f"- `{row['ID']}` score={row['priority_score']} "
            f"[strictHC={row['strict_highconfidence']}, HIC={row['has_hic_panel']}] "
            f"L={row['left_best_event']} | R={row['right_best_event']}"
        )
    lines.append("")
    top_ids = set(top_df["ID"].astype(str))
    top_events = events_df[events_df["INV_ID"].astype(str).isin(top_ids)].copy()
    if not top_events.empty:
        lines.append("## Top-candidate breakpoint events")
        lines.append("")
        for _, row in top_events.sort_values(
            ["INV_ID", "boundary_side", "relation", "distance_to_gene_bp", "gene_id"]
        ).iterrows():
            lines.append(
                f"- `{row['INV_ID']}` {row['boundary_side']} `{row['gene_id']}` "
                f"relation={row['relation']} dist={row['distance_to_gene_bp']}bp "
                f"expr_mean={format_expr(row['mean_expr'])}"
                + (
                    f" exon_boundary_dist={int(row['boundary_dist'])}bp"
                    if pd.notna(row["boundary_dist"])
                    else ""
                )
            )
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    expr_by_gene = {} if args.skip_expression else parse_gene_counts(args.gene_counts)
    tx_by_gene, tx_by_chr = load_transcripts(args.t1_bed, expr_by_gene, args.promoter_bp)
    exon_map, intron_map = load_feature_intervals(args.gtf)

    inv_df = pd.read_csv(args.inversions, sep="\t")
    if "HighConfidence" in inv_df.columns:
        inv_df["strict_highconfidence"] = inv_df["HighConfidence"].map(parse_bool)
    else:
        inv_df["strict_highconfidence"] = False
    if args.strict_only:
        inv_df = inv_df[inv_df["strict_highconfidence"]].copy()
    chrom_col = f"{args.coord_prefix}Chr"
    start_col = f"{args.coord_prefix}Start"
    end_col = f"{args.coord_prefix}End"
    inv_df = inv_df.sort_values([chrom_col, start_col, end_col, "ID"]).reset_index(drop=True)

    hic_df = pd.read_csv(args.hic_summary, sep="\t") if args.hic_summary.exists() else pd.DataFrame()
    hic_map = {
        str(row["ID"]): {
            "hic_inter_contacts": float(row["Inter_contacts"]),
            "hic_bin_size": int(row["bin_size"]),
            "hic_flank_bp": int(row["flank_bp"]),
        }
        for _, row in hic_df.iterrows()
    }

    overview_df = (
        pd.read_csv(args.inv_overview, sep="\t")
        if (not args.skip_expression and args.inv_overview and args.inv_overview.exists())
        else pd.DataFrame()
    )
    overview_map = {
        str(row["INV_ID"]): {
            "interior_genes_with_expression": int(row["GenesWithExpression"]),
            "interior_mean_expr_alllibs": float(row["MeanExpr_AcrossSamples"]),
            "interior_max_expr_alllibs": float(row["MaxExpr_AcrossSamples"]),
        }
        for _, row in overview_df.iterrows()
    }

    screen_rows: List[dict] = []
    event_rows: List[dict] = []

    for _, inv in inv_df.iterrows():
        inv_id = str(inv["ID"])
        chrom = str(inv[chrom_col])
        left_bp = int(inv[start_col])
        right_bp = int(inv[end_col])
        transcripts = tx_by_chr.get(chrom, [])

        left_summary, left_events = analyze_breakpoint(
            chrom=chrom,
            pos=left_bp,
            side="left",
            transcripts=transcripts,
            exon_map=exon_map,
            intron_map=intron_map,
            boundary_bp=args.boundary_bp,
            flank_bp=args.flank_bp,
        )
        right_summary, right_events = analyze_breakpoint(
            chrom=chrom,
            pos=right_bp,
            side="right",
            transcripts=transcripts,
            exon_map=exon_map,
            intron_map=intron_map,
            boundary_bp=args.boundary_bp,
            flank_bp=args.flank_bp,
        )

        for event in left_events + right_events:
            event_rows.append({"INV_ID": inv_id, **event})

        overview = overview_map.get(inv_id, {})
        hic = hic_map.get(inv_id, {})

        row = {
            "ID": inv_id,
            "ScreenLabel": args.label,
            "CoordPrefix": args.coord_prefix,
            "ScreenChr": chrom,
            "ScreenStart": left_bp,
            "ScreenEnd": right_bp,
            "RefChr": str(inv["RefChr"]),
            "RefStart": int(inv["RefStart"]),
            "RefEnd": int(inv["RefEnd"]),
            "RefLen": int(inv["RefLen"]),
            "QryChr": str(inv["QryChr"]),
            "QryStart": int(inv["QryStart"]),
            "QryEnd": int(inv["QryEnd"]),
            "QryLen": int(inv["QryLen"]),
            "MinLen": int(inv["MinLen"]),
            "RefRepeatFrac": float(inv["RefRepeatFrac"]),
            "QryRepeatFrac": float(inv["QryRepeatFrac"]),
            "BPMedianRepeatFrac": float(inv["BPMedianRepeatFrac"]),
            "strict_highconfidence": bool(inv["strict_highconfidence"]),
            "has_hic_panel": inv_id in hic_map,
            "hic_inter_contacts": hic.get("hic_inter_contacts", math.nan),
            "left_best_event": left_summary["best_event"],
            "right_best_event": right_summary["best_event"],
            "left_n_overlap_genes": left_summary["n_overlap_genes"],
            "right_n_overlap_genes": right_summary["n_overlap_genes"],
            "left_n_promoter_hits": left_summary["n_promoter_hits"],
            "right_n_promoter_hits": right_summary["n_promoter_hits"],
            "left_n_exon_hits": left_summary["n_exon_hits"],
            "right_n_exon_hits": right_summary["n_exon_hits"],
            "left_n_intron_hits": left_summary["n_intron_hits"],
            "right_n_intron_hits": right_summary["n_intron_hits"],
            "left_n_boundary_near_hits": left_summary["n_boundary_near_hits"],
            "right_n_boundary_near_hits": right_summary["n_boundary_near_hits"],
            "left_nearest_left_gene": left_summary["nearest_left_gene"],
            "left_nearest_left_dist_bp": left_summary["nearest_left_dist_bp"],
            "left_nearest_right_gene": left_summary["nearest_right_gene"],
            "left_nearest_right_dist_bp": left_summary["nearest_right_dist_bp"],
            "right_nearest_left_gene": right_summary["nearest_left_gene"],
            "right_nearest_left_dist_bp": right_summary["nearest_left_dist_bp"],
            "right_nearest_right_gene": right_summary["nearest_right_gene"],
            "right_nearest_right_dist_bp": right_summary["nearest_right_dist_bp"],
            "interior_genes_with_expression": overview.get("interior_genes_with_expression", math.nan),
            "interior_mean_expr_alllibs": overview.get("interior_mean_expr_alllibs", math.nan),
            "interior_max_expr_alllibs": overview.get("interior_max_expr_alllibs", math.nan),
        }
        row["any_overlap_gene"] = (row["left_n_overlap_genes"] + row["right_n_overlap_genes"]) > 0
        row["any_exon_hit"] = (row["left_n_exon_hits"] + row["right_n_exon_hits"]) > 0
        row["any_intron_hit"] = (row["left_n_intron_hits"] + row["right_n_intron_hits"]) > 0
        row["any_promoter_hit"] = (row["left_n_promoter_hits"] + row["right_n_promoter_hits"]) > 0
        row["any_boundary_near_hit"] = (
            row["left_n_boundary_near_hits"] + row["right_n_boundary_near_hits"]
        ) > 0
        row["any_expressed_boundary_gene"] = bool(
            left_summary["has_expressed_boundary_gene"] or right_summary["has_expressed_boundary_gene"]
        )
        row["any_expressed_flank_gene"] = bool(
            left_summary["has_expressed_flank_gene"] or right_summary["has_expressed_flank_gene"]
        )
        row["max_boundary_related_mean_expr"] = max(
            left_summary["max_boundary_related_mean_expr"],
            right_summary["max_boundary_related_mean_expr"],
        )
        row["nearest_gene_within_2kb"] = any(
            pd.notna(x) and float(x) <= 2000
            for x in [
                row["left_nearest_left_dist_bp"],
                row["left_nearest_right_dist_bp"],
                row["right_nearest_left_dist_bp"],
                row["right_nearest_right_dist_bp"],
            ]
        )
        row["priority_score"] = 0
        screen_rows.append(row)

    screen_df = pd.DataFrame(screen_rows)
    screen_df["priority_score"] = screen_df.apply(compute_priority_score, axis=1)
    screen_df = screen_df.sort_values(
        ["priority_score", "strict_highconfidence", "has_hic_panel", "hic_inter_contacts", "MinLen", "ID"],
        ascending=[False, False, False, False, False, True],
    ).reset_index(drop=True)
    top_df = screen_df.head(args.top_n).copy()
    events_df = pd.DataFrame(event_rows)
    if not events_df.empty:
        events_df = events_df.sort_values(
            ["INV_ID", "boundary_side", "relation", "distance_to_gene_bp", "gene_id"]
        ).reset_index(drop=True)

    all_tsv = args.outdir / "all_inversion_breakpoint_priority_screen.tsv"
    top_tsv = args.outdir / f"top{args.top_n}_inversion_breakpoint_priority_screen.tsv"
    events_tsv = args.outdir / "all_breakpoint_candidate_events.tsv"
    md_path = args.outdir / "breakpoint_priority_screen.summary.md"

    screen_df.to_csv(all_tsv, sep="\t", index=False)
    top_df.to_csv(top_tsv, sep="\t", index=False)
    events_df.to_csv(events_tsv, sep="\t", index=False)
    md_path.write_text(build_markdown_summary(screen_df, top_df, events_df, args))

    print(f"Wrote: {all_tsv}")
    print(f"Wrote: {top_tsv}")
    print(f"Wrote: {events_tsv}")
    print(f"Wrote: {md_path}")
    print("")
    print("Top candidates:")
    display_cols = [
        "ID",
        "priority_score",
        "strict_highconfidence",
        "has_hic_panel",
        "any_exon_hit",
        "any_promoter_hit",
        "any_intron_hit",
        "any_boundary_near_hit",
        "left_best_event",
        "right_best_event",
    ]
    print(top_df[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
