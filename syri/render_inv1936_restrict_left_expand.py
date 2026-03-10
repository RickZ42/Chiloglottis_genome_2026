#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile

from PIL import Image, ImageColor, ImageDraw, ImageFont


TRACK_COLORS = {"H1": "#4C78A8", "H2": "#E45756"}
ORIENTATION_COLORS = {"+": "#4C78A8", "-": "#54A24B"}
NON_INV_LINK_COLOR = "#C7C7C7"
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


def maybe_copy_extra(extra_bed: Path, out_extra: Path) -> None:
    if extra_bed.resolve() != out_extra.resolve():
        copyfile(extra_bed, out_extra)


def default_output_tag(left_expand_bp: int) -> str:
    if left_expand_bp % 1_000_000 == 0:
        return f"restrict_left{left_expand_bp // 1_000_000}Mb"
    return f"restrict_left{left_expand_bp}bp"


def load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    candidates = []
    if bold:
        candidates.extend(
            [
                "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
            ]
        )
    candidates.extend(
        [
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/dejavu/DejaVuSans.ttf",
        ]
    )
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def rgba(color: str, alpha: int = 255) -> tuple[int, int, int, int]:
    r, g, b = ImageColor.getrgb(color)
    return (r, g, b, alpha)


def x_from_bp(window: Window, bp: float, x0: float, x1: float) -> float:
    span = max(1, window.end - window.start)
    return x0 + (bp - window.start) / span * (x1 - x0)


def cubic_bezier(p0, p1, p2, p3, n=40):
    pts = []
    for i in range(n + 1):
        t = i / n
        x = ((1 - t) ** 3) * p0[0] + 3 * ((1 - t) ** 2) * t * p1[0] + 3 * (1 - t) * (t**2) * p2[0] + (t**3) * p3[0]
        y = ((1 - t) ** 3) * p0[1] + 3 * ((1 - t) ** 2) * t * p1[1] + 3 * (1 - t) * (t**2) * p2[1] + (t**3) * p3[1]
        pts.append((x, y))
    return pts


def draw_ribbon(draw: ImageDraw.ImageDraw, top_start, top_end, bottom_start, bottom_end, y_top, y_bottom, color, alpha):
    ymid_top = (y_top + y_bottom) / 2 + 120
    ymid_bottom = (y_top + y_bottom) / 2 - 30
    curve1 = cubic_bezier((top_start, y_top), (top_start, ymid_top), (bottom_start, ymid_top), (bottom_start, y_bottom))
    curve2 = cubic_bezier((bottom_end, y_bottom), (bottom_end, ymid_bottom), (top_end, ymid_bottom), (top_end, y_top))
    polygon = curve1 + [(bottom_end, y_bottom)] + curve2 + [(top_start, y_top)]
    draw.polygon([(int(x), int(y)) for x, y in polygon], fill=rgba(color, alpha))


def draw_arrow(draw: ImageDraw.ImageDraw, start: float, end: float, y: float, color: str):
    thickness = 16
    dx = end - start
    if abs(dx) < 12:
        dx = 12 if dx >= 0 else -12
        end = start + dx
    head = min(max(abs(dx) * 0.35, 10), 22)
    if dx >= 0:
        body_end = end - head
        pts = [
            (start, y - thickness // 3),
            (body_end, y - thickness // 3),
            (body_end, y - thickness // 2),
            (end, y),
            (body_end, y + thickness // 2),
            (body_end, y + thickness // 3),
            (start, y + thickness // 3),
        ]
    else:
        body_start = end + head
        pts = [
            (start, y - thickness // 3),
            (body_start, y - thickness // 3),
            (body_start, y - thickness // 2),
            (end, y),
            (body_start, y + thickness // 2),
            (body_start, y + thickness // 3),
            (start, y + thickness // 3),
        ]
    draw.polygon([(int(x), int(y)) for x, y in pts], fill=rgba(color))


def draw_track_label(draw: ImageDraw.ImageDraw, label: str, y: int, color: str, bold_font, small_font):
    main, _, tail = label.partition(": ")
    draw.text((20, y - 44), main, fill=rgba(color), font=bold_font)
    if tail:
        draw.text((20, y - 12), tail, fill=rgba(color), font=small_font)


def draw_main_panel(
    img: Image.Image,
    h1_features: list[Feature],
    h2_features: list[Feature],
    h1_window: Window,
    h2_window: Window,
    h1_label: str,
    h2_label: str,
    blocks: list[tuple[str, str, str | None]],
    extra_regions: list[Feature],
    with_legend: bool,
) -> None:
    x0 = 300
    x1 = 1760 if with_legend else 2320
    y_top = 600 if with_legend else 720
    y_bottom = 1230 if with_legend else 1380

    base_draw = ImageDraw.Draw(img)
    overlay = Image.new("RGBA", img.size, (255, 255, 255, 0))
    overlay_draw = ImageDraw.Draw(overlay)
    label_font = load_font(34, bold=True)
    small_font = load_font(28)

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
        alpha = 190 if highlight_color else 110
        draw_ribbon(overlay_draw, top_start, top_end, bottom_start, bottom_end, y_top, y_bottom, color, alpha)

    img.alpha_composite(overlay)
    base_draw = ImageDraw.Draw(img)
    draw_track_label(base_draw, h1_label, y_top, TRACK_COLORS["H1"], label_font, small_font)
    draw_track_label(base_draw, h2_label, y_bottom, TRACK_COLORS["H2"], label_font, small_font)

    for y in (y_top, y_bottom):
        base_draw.line([(x0, y), (x1, y)], fill=rgba("#666666"), width=6)

    for feat in extra_regions:
        if feat.chrom == f"H1_{h1_window.chrom}":
            start = x_from_bp(h1_window, feat.start, x0, x1)
            end = x_from_bp(h1_window, feat.end, x0, x1)
            base_draw.line([(start, y_top), (end, y_top)], fill=rgba(INV_INTERVAL_COLOR), width=18)
        elif feat.chrom == f"H2_{h2_window.chrom}":
            start = x_from_bp(h2_window, feat.start, x0, x1)
            end = x_from_bp(h2_window, feat.end, x0, x1)
            base_draw.line([(start, y_bottom), (end, y_bottom)], fill=rgba(INV_INTERVAL_COLOR), width=18)

    for feat in h1_features:
        start = x_from_bp(h1_window, feat.start, x0, x1)
        end = x_from_bp(h1_window, feat.end, x0, x1)
        draw_arrow(base_draw, start, end, y_top, ORIENTATION_COLORS.get(feat.strand, TRACK_COLORS["H1"]))
    for feat in h2_features:
        start = x_from_bp(h2_window, feat.start, x0, x1)
        end = x_from_bp(h2_window, feat.end, x0, x1)
        draw_arrow(base_draw, start, end, y_bottom, ORIENTATION_COLORS.get(feat.strand, TRACK_COLORS["H2"]))


def draw_legend(img: Image.Image, h1_window: Window, h2_window: Window, left_expand_bp: int, highlight_color: str):
    draw = ImageDraw.Draw(img)
    title_font = load_font(50, bold=True)
    item_font = load_font(34)
    note_font = load_font(30)

    x0 = 2200
    draw.text((x0, 120), "Legend", fill=rgba("#111111"), font=title_font)

    draw.line([(x0 + 20, 300), (x0 + 170, 300)], fill=rgba(highlight_color), width=12)
    draw.text((x0 + 210, 260), "Inversion links\n(core + 1 left gene)", fill=rgba("#222222"), font=item_font, spacing=6)

    draw.line([(x0 + 20, 540), (x0 + 170, 540)], fill=rgba("#9E9E9E"), width=12)
    draw.text((x0 + 210, 500), "Non-inversion\ncollinear links", fill=rgba("#222222"), font=item_font, spacing=6)

    draw.line([(x0 + 20, 780), (x0 + 170, 780)], fill=rgba(INV_INTERVAL_COLOR), width=28)
    draw.text((x0 + 210, 740), "Inversion interval\n(track highlight)", fill=rgba("#222222"), font=item_font, spacing=6)

    draw_arrow(draw, x0 + 20, x0 + 170, 1030, ORIENTATION_COLORS["+"])
    draw.text((x0 + 210, 990), "Gene arrow direction\n= feature orientation", fill=rgba("#222222"), font=item_font, spacing=6)

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
    box = (x0 - 10, 1450, 3050, 1840)
    draw.rounded_rectangle(box, radius=20, fill=rgba("#F5F5F5"), outline=rgba("#7A7A7A"), width=3)
    draw.multiline_text((x0 + 10, 1475), window_note, fill=rgba("#222222"), font=note_font, spacing=10)


def render_png(
    png_path: Path,
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
    size = (3075, 1906) if with_legend else (2700, 2100)
    img = Image.new("RGBA", size, (255, 255, 255, 255))
    draw_main_panel(
        img,
        h1_features,
        h2_features,
        h1_window,
        h2_window,
        h1_label,
        h2_label,
        blocks,
        extra_regions,
        with_legend,
    )
    if with_legend:
        highlight_color = next((color for _, _, color in blocks if color), "lightskyblue")
        draw_legend(img, h1_window, h2_window, left_expand_bp, highlight_color)

    img.save(png_path)


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
        for feat in h1_features:
            out.write(f"H1_{feat.chrom}\t{feat.start}\t{feat.end}\t{feat.name}\t{feat.score}\t{feat.strand}\n")
        for feat in h2_features:
            out.write(f"H2_{feat.chrom}\t{feat.start}\t{feat.end}\t{feat.name}\t{feat.score}\t{feat.strand}\n")
    write_layout(layout_path, h1_window, h2_window, args.track_right_label)
    write_summary(summary_path, h1_window, h2_window, len(h1_features), len(h2_features))
    maybe_copy_extra(args.extra_bed, extra_path)

    blocks = parse_blocks(args.blocks)
    extra_regions = read_bed(extra_path)
    h1_label = (
        f"H1 {h1_window.chrom} ({mb_text(h1_window.start)}-{mb_text(h1_window.end)} Mb): "
        f"{args.track_right_label}"
    )
    h2_label = (
        f"H2 {h2_window.chrom} ({mb_text(h2_window.start)}-{mb_text(h2_window.end)} Mb): "
        f"{args.track_right_label}"
    )

    base_prefix = args.outdir / (
        f"H1_vs_H2_{args.inv_id}_jcvi_microsynteny.{output_tag}.strictAnchor.highlightINVplusLeft"
    )
    render_png(
        png_path=Path(str(base_prefix) + ".nolabels.png"),
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
    render_png(
        png_path=Path(str(base_prefix) + ".WITH_LEGEND.png"),
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
