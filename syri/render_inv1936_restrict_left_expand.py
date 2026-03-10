#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, FancyArrowPatch, PathPatch
from matplotlib.path import Path as MplPath


TRACK_COLORS = {"H1": "#4C78A8", "H2": "#E45756"}
ORIENTATION_COLORS = {"+": "b", "-": "g"}
NON_INV_LINK_COLOR = "#BFBFBF"
INV_INTERVAL_COLOR = "#F57C00"


@dataclass
class Window:
    chrom: str
    start: int
    end: int


@dataclass
class Feature:
    chrom: str
    start: int
    end: int
    name: str
    score: str
    strand: str

    @property
    def mid(self) -> float:
        return (self.start + self.end) / 2.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Expand the INV1936 restrict window to the left and redraw the microsynteny figure."
    )
    p.add_argument("--inv-id", default="INV1936")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--restrict-summary", type=Path, required=True)
    p.add_argument("--full-h1-bed", type=Path, required=True)
    p.add_argument("--full-h2-bed", type=Path, required=True)
    p.add_argument("--blocks", type=Path, required=True)
    p.add_argument("--extra-bed", type=Path, required=True)
    p.add_argument("--left-expand-bp", type=int, required=True)
    p.add_argument("--output-tag", default=None)
    p.add_argument("--track-right-label", default="left flank | inversion core | right flank")
    return p.parse_args()


def read_summary(path: Path) -> dict[str, str]:
    with path.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return {row["metric"]: row["value"] for row in rows}


def parse_window(text: str) -> Window:
    chrom, coords = text.split(":")
    start_text, end_text = coords.split("-")
    return Window(chrom=chrom, start=int(start_text), end=int(end_text))


def read_bed(path: Path) -> list[Feature]:
    feats: list[Feature] = []
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name, score, strand = line.rstrip("\n").split("\t")[:6]
            feats.append(
                Feature(
                    chrom=chrom,
                    start=int(start),
                    end=int(end),
                    name=name,
                    score=score,
                    strand=strand,
                )
            )
    return feats


def overlaps(feature: Feature, window: Window) -> bool:
    return feature.chrom == window.chrom and feature.end >= window.start and feature.start <= window.end


def filter_bed(features: list[Feature], window: Window) -> list[Feature]:
    return [feat for feat in features if overlaps(feat, window)]


def write_bed(path: Path, features: list[Feature], prefix: str | None = None) -> None:
    with path.open("w") as out:
        for feat in features:
            chrom = f"{prefix}_{feat.chrom}" if prefix else feat.chrom
            out.write(
                "\t".join(
                    [
                        chrom,
                        str(feat.start),
                        str(feat.end),
                        feat.name,
                        feat.score,
                        feat.strand,
                    ]
                )
                + "\n"
            )


def mb_text(value: int) -> str:
    return f"{value / 1_000_000:.2f}"


def write_layout(path: Path, h1_window: Window, h2_window: Window, label_suffix: str) -> None:
    with path.open("w") as out:
        out.write("# x, y, rotation, ha, va, color, ratio, label\n")
        out.write(
            "0.5, 0.70, 0, left, center, {color}, 1, H1 {chrom} ({start}-{end} Mb): {suffix}\n".format(
                color=TRACK_COLORS["H1"],
                chrom=h1_window.chrom,
                start=mb_text(h1_window.start),
                end=mb_text(h1_window.end),
                suffix=label_suffix,
            )
        )
        out.write(
            "0.5, 0.30, 0, left, center, {color}, 1, H2 {chrom} ({start}-{end} Mb): {suffix}\n".format(
                color=TRACK_COLORS["H2"],
                chrom=h2_window.chrom,
                start=mb_text(h2_window.start),
                end=mb_text(h2_window.end),
                suffix=label_suffix,
            )
        )
        out.write("# edges\n")
        out.write("e, 0, 1\n")


def write_summary(path: Path, h1_window: Window, h2_window: Window, h1_count: int, h2_count: int) -> None:
    rows = [
        ("H1_window", f"{h1_window.chrom}:{h1_window.start}-{h1_window.end}"),
        ("H2_window", f"{h2_window.chrom}:{h2_window.start}-{h2_window.end}"),
        ("H1_genes_in_window", str(h1_count)),
        ("H2_genes_in_window", str(h2_count)),
    ]
    with path.open("w") as out:
        out.write("metric\tvalue\n")
        for metric, value in rows:
            out.write(f"{metric}\t{value}\n")


def parse_blocks(path: Path) -> list[tuple[str, str, str | None]]:
    pairs: list[tuple[str, str, str | None]] = []
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            left, right = line.rstrip("\n").split("\t")[:2]
            color = None
            gene_left = left
            if "*" in left:
                maybe_color, maybe_gene = left.split("*", 1)
                if maybe_gene:
                    color = maybe_color
                    gene_left = maybe_gene
            pairs.append((gene_left, right, color))
    return pairs


def x_from_bp(window: Window, bp: float, x0: float, x1: float) -> float:
    span = max(1, window.end - window.start)
    return x0 + (bp - window.start) / span * (x1 - x0)


def add_ribbon(
    ax,
    top_start: float,
    top_end: float,
    bottom_start: float,
    bottom_end: float,
    y_top: float,
    y_bottom: float,
    color: str,
    alpha: float,
    lw: float,
    zorder: int,
) -> None:
    ymid_top = (y_top + y_bottom) / 2 + 0.08
    ymid_bottom = (y_top + y_bottom) / 2 - 0.02
    verts = [
        (top_start, y_top),
        (top_start, ymid_top),
        (bottom_start, ymid_top),
        (bottom_start, y_bottom),
        (bottom_end, y_bottom),
        (bottom_end, ymid_bottom),
        (top_end, ymid_bottom),
        (top_end, y_top),
        (top_start, y_top),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    ax.add_patch(
        PathPatch(
            MplPath(verts, codes),
            facecolor=color,
            edgecolor=color,
            alpha=alpha,
            linewidth=lw,
            zorder=zorder,
        )
    )


def draw_gene_track(
    ax,
    features: list[Feature],
    window: Window,
    y: float,
    x0: float,
    x1: float,
    extra_regions: list[Feature],
) -> None:
    ax.plot([x0, x1], [y, y], color="#666666", linewidth=2.0, zorder=3)
    for feat in extra_regions:
        if feat.chrom != f"H1_{window.chrom}" and feat.chrom != f"H2_{window.chrom}":
            continue
        start = x_from_bp(window, feat.start, x0, x1)
        end = x_from_bp(window, feat.end, x0, x1)
        ax.plot([start, end], [y, y], color=INV_INTERVAL_COLOR, linewidth=6.0, zorder=4, solid_capstyle="butt")

    for feat in features:
        start = x_from_bp(window, feat.start, x0, x1)
        end = x_from_bp(window, feat.end, x0, x1)
        dx = end - start
        color = ORIENTATION_COLORS.get(feat.strand, "b")
        if abs(dx) < 0.004:
            dx = 0.004 if feat.strand == "+" else -0.004
        head_length = min(max(abs(dx) * 0.35, 0.004), 0.012)
        ax.add_patch(
            FancyArrow(
                start,
                y,
                dx,
                0,
                width=0.0045,
                head_width=0.018,
                head_length=head_length,
                length_includes_head=True,
                color=color,
                linewidth=0,
                zorder=5,
            )
        )


def draw_track_label(ax, label: str, y: float, color: str) -> None:
    main, _, tail = label.partition(": ")
    ax.text(0.015, y + 0.028, main, ha="left", va="bottom", fontsize=10, fontweight="bold", color=color)
    if tail:
        ax.text(0.015, y - 0.006, tail, ha="left", va="top", fontsize=8.5, color=color)


def render_main_panel(
    ax,
    h1_features: list[Feature],
    h2_features: list[Feature],
    h1_window: Window,
    h2_window: Window,
    h1_label: str,
    h2_label: str,
    blocks: list[tuple[str, str, str | None]],
    extra_regions: list[Feature],
) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    x0, x1 = 0.28, 0.94
    y_top, y_bottom = 0.70, 0.30
    draw_track_label(ax, h1_label, y_top, TRACK_COLORS["H1"])
    draw_track_label(ax, h2_label, y_bottom, TRACK_COLORS["H2"])

    draw_gene_track(ax, h1_features, h1_window, y_top, x0, x1, extra_regions)
    draw_gene_track(ax, h2_features, h2_window, y_bottom, x0, x1, extra_regions)

    h1_map = {feat.name: feat for feat in h1_features}
    h2_map = {feat.name: feat for feat in h2_features}

    for left_name, right_name, highlight_color in blocks:
        h1_feat = h1_map.get(left_name)
        h2_feat = h2_map.get(right_name)
        if not h1_feat or not h2_feat:
            continue
        top_start = x_from_bp(h1_window, h1_feat.start, x0, x1)
        top_end = x_from_bp(h1_window, h1_feat.end, x0, x1)
        bottom_start = x_from_bp(h2_window, h2_feat.start, x0, x1)
        bottom_end = x_from_bp(h2_window, h2_feat.end, x0, x1)
        color = highlight_color or NON_INV_LINK_COLOR
        alpha = 0.9 if highlight_color else 0.45
        lw = 1.5 if highlight_color else 1.0
        zorder = 2 if highlight_color else 1
        add_ribbon(ax, top_start, top_end, bottom_start, bottom_end, y_top, y_bottom, color, alpha, lw, zorder)


def draw_legend(legend_ax, h1_window: Window, h2_window: Window, left_expand_bp: int, highlight_color: str) -> None:
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.axis("off")

    legend_ax.text(0.04, 0.92, "Legend", fontsize=18, fontweight="bold", ha="left", va="center")

    legend_ax.plot([0.08, 0.32], [0.78, 0.78], color=highlight_color, linewidth=4.5, solid_capstyle="round")
    legend_ax.text(0.36, 0.78, "Inversion links\n(core + 1 left gene)", fontsize=11.5, ha="left", va="center")

    legend_ax.plot([0.08, 0.32], [0.62, 0.62], color="#9E9E9E", linewidth=4.0, solid_capstyle="round")
    legend_ax.text(0.36, 0.62, "Non-inversion\ncollinear links", fontsize=11.5, ha="left", va="center")

    legend_ax.plot([0.08, 0.32], [0.46, 0.46], color=INV_INTERVAL_COLOR, linewidth=10.0, solid_capstyle="butt")
    legend_ax.text(0.36, 0.46, "Inversion interval\n(track highlight)", fontsize=11.5, ha="left", va="center")

    legend_ax.add_patch(
        FancyArrow(
            0.10,
            0.28,
            0.20,
            0,
            width=0.01,
            head_width=0.05,
            head_length=0.04,
            length_includes_head=True,
            color="b",
            linewidth=0,
        )
    )
    legend_ax.text(0.36, 0.28, "Gene arrow direction\n= feature orientation", fontsize=11.5, ha="left", va="center")

    mb_value = left_expand_bp / 1_000_000.0
    if mb_value.is_integer():
        mb_label = f"{int(mb_value)} Mb"
    else:
        mb_label = f"{mb_value:.2f} Mb"
    window_note = "\n".join(
        [
            f"Display window (left expanded by {mb_label}):",
            f"H1 {h1_window.chrom}: {h1_window.start:,}-{h1_window.end:,}",
            f"H2 {h2_window.chrom}: {h2_window.start:,}-{h2_window.end:,}",
        ]
    )
    legend_ax.text(
        0.04,
        0.08,
        window_note,
        fontsize=10.5,
        ha="left",
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#F7F7F7", edgecolor="#777777", linewidth=1.0),
    )


def render_figure(
    outpath: Path,
    h1_features: list[Feature],
    h2_features: list[Feature],
    h1_window: Window,
    h2_window: Window,
    h1_label: str,
    h2_label: str,
    blocks: list[tuple[str, str, str | None]],
    extra_regions: list[Feature],
    left_expand_bp: int,
    with_legend: bool,
) -> None:
    highlight_color = next((color for _, _, color in blocks if color), "lightskyblue")

    if with_legend:
        fig = plt.figure(figsize=(12.3, 7.6))
        main_ax = fig.add_axes((0.03, 0.08, 0.66, 0.84))
        legend_ax = fig.add_axes((0.73, 0.08, 0.25, 0.84))
        draw_legend(legend_ax, h1_window, h2_window, left_expand_bp, highlight_color)
    else:
        fig = plt.figure(figsize=(9.0, 7.0))
        main_ax = fig.add_axes((0.03, 0.06, 0.94, 0.88))

    render_main_panel(
        main_ax,
        h1_features,
        h2_features,
        h1_window,
        h2_window,
        h1_label,
        h2_label,
        blocks,
        extra_regions,
    )
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def maybe_copy_extra(extra_bed: Path, out_extra: Path) -> None:
    if extra_bed.resolve() != out_extra.resolve():
        copyfile(extra_bed, out_extra)


def default_output_tag(left_expand_bp: int) -> str:
    if left_expand_bp % 1_000_000 == 0:
        return f"restrict_left{left_expand_bp // 1_000_000}Mb"
    return f"restrict_left{left_expand_bp}bp"


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    output_tag = args.output_tag or default_output_tag(args.left_expand_bp)
    summary = read_summary(args.restrict_summary)
    h1_base = parse_window(summary["H1_window"])
    h2_base = parse_window(summary["H2_window"])
    h1_window = Window(h1_base.chrom, max(0, h1_base.start - args.left_expand_bp), h1_base.end)
    h2_window = Window(h2_base.chrom, max(0, h2_base.start - args.left_expand_bp), h2_base.end)

    h1_full = read_bed(args.full_h1_bed)
    h2_full = read_bed(args.full_h2_bed)
    h1_features = filter_bed(h1_full, h1_window)
    h2_features = filter_bed(h2_full, h2_window)

    h1_bed = args.outdir / f"H1_{args.inv_id}_{output_tag}.t1.bed"
    h2_bed = args.outdir / f"H2_{args.inv_id}_{output_tag}.t1.bed"
    h1_prefixed_bed = args.outdir / f"H1_{args.inv_id}_{output_tag}.prefixed.t1.bed"
    h2_prefixed_bed = args.outdir / f"H2_{args.inv_id}_{output_tag}.prefixed.t1.bed"
    merged_bed = args.outdir / f"H1_H2_{args.inv_id}_{output_tag}.prefixed.merged.bed"
    layout_path = args.outdir / f"H1_H2_{args.inv_id}_{output_tag}.layout"
    summary_path = args.outdir / f"{args.inv_id}_{output_tag}.summary.tsv"
    extra_path = args.outdir / f"H1_H2_{args.inv_id}_{output_tag}.prefixed.inversion_core.extra.bed"

    write_bed(h1_bed, h1_features)
    write_bed(h2_bed, h2_features)
    write_bed(h1_prefixed_bed, h1_features, prefix="H1")
    write_bed(h2_prefixed_bed, h2_features, prefix="H2")
    with merged_bed.open("w") as out:
        write_bed(h1_prefixed_bed, h1_features, prefix="H1")
    # Re-open so the merged BED preserves H1 then H2 order.
    with merged_bed.open("w") as out:
        for feat in h1_features:
            out.write(f"H1_{feat.chrom}\t{feat.start}\t{feat.end}\t{feat.name}\t{feat.score}\t{feat.strand}\n")
        for feat in h2_features:
            out.write(f"H2_{feat.chrom}\t{feat.start}\t{feat.end}\t{feat.name}\t{feat.score}\t{feat.strand}\n")

    h1_label = (
        f"H1 {h1_window.chrom} ({mb_text(h1_window.start)}-{mb_text(h1_window.end)} Mb): "
        f"{args.track_right_label}"
    )
    h2_label = (
        f"H2 {h2_window.chrom} ({mb_text(h2_window.start)}-{mb_text(h2_window.end)} Mb): "
        f"{args.track_right_label}"
    )
    write_layout(layout_path, h1_window, h2_window, args.track_right_label)
    write_summary(summary_path, h1_window, h2_window, len(h1_features), len(h2_features))
    maybe_copy_extra(args.extra_bed, extra_path)

    blocks = parse_blocks(args.blocks)
    extra_regions = read_bed(extra_path)
    base_prefix = args.outdir / (
        f"H1_vs_H2_{args.inv_id}_jcvi_microsynteny.{output_tag}.strictAnchor.highlightINVplusLeft"
    )
    render_figure(
        outpath=Path(str(base_prefix) + ".nolabels.png"),
        h1_features=h1_features,
        h2_features=h2_features,
        h1_window=h1_window,
        h2_window=h2_window,
        h1_label=h1_label,
        h2_label=h2_label,
        blocks=blocks,
        extra_regions=extra_regions,
        left_expand_bp=args.left_expand_bp,
        with_legend=False,
    )
    render_figure(
        outpath=Path(str(base_prefix) + ".WITH_LEGEND.png"),
        h1_features=h1_features,
        h2_features=h2_features,
        h1_window=h1_window,
        h2_window=h2_window,
        h1_label=h1_label,
        h2_label=h2_label,
        blocks=blocks,
        extra_regions=extra_regions,
        left_expand_bp=args.left_expand_bp,
        with_legend=True,
    )
    render_figure(
        outpath=Path(str(base_prefix) + ".nolabels.pdf"),
        h1_features=h1_features,
        h2_features=h2_features,
        h1_window=h1_window,
        h2_window=h2_window,
        h1_label=h1_label,
        h2_label=h2_label,
        blocks=blocks,
        extra_regions=extra_regions,
        left_expand_bp=args.left_expand_bp,
        with_legend=False,
    )
    render_figure(
        outpath=Path(str(base_prefix) + ".WITH_LEGEND.pdf"),
        h1_features=h1_features,
        h2_features=h2_features,
        h1_window=h1_window,
        h2_window=h2_window,
        h1_label=h1_label,
        h2_label=h2_label,
        blocks=blocks,
        extra_regions=extra_regions,
        left_expand_bp=args.left_expand_bp,
        with_legend=True,
    )


if __name__ == "__main__":
    main()
