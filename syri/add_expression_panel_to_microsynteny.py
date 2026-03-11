#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from PIL import Image, ImageChops, ImageColor, ImageDraw, ImageFont


META_COLUMNS = {
    "INV_ID",
    "RefChr",
    "RefStart",
    "RefEnd",
    "QryChr",
    "QryStart",
    "QryEnd",
    "GeneID",
    "GeneChr",
    "GeneStart",
    "GeneEnd",
}

FEATURECOUNTS_META_COLUMNS = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}


@dataclass
class Feature:
    chrom: str
    start: int
    end: int
    name: str
    strand: str


@dataclass
class ExpressionGene:
    gene_id: str
    mean_expr: float


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Add an H1 RNA expression panel above an existing microsynteny figure."
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--background-png", type=Path, help="Existing microsynteny PNG.")
    src.add_argument("--background-pdf", type=Path, help="Microsynteny PDF; page 1 will be rasterized with ImageMagick.")
    p.add_argument("--pdf-density", type=int, default=200, help="Rasterization density for --background-pdf.")
    p.add_argument("--h1-bed", type=Path, required=True, help="Top-track H1 BED used by the figure.")
    p.add_argument("--expr-table", type=Path, required=True, help="Per-gene expression table.")
    p.add_argument(
        "--expr-metric",
        choices=["counts", "tpm"],
        default="counts",
        help="Expression metric to plot. TPM is available for featureCounts-style tables with Length column.",
    )
    p.add_argument("--inv-id", default="INV1936", help="Inversion ID present in the expression table.")
    p.add_argument("--window-info", type=Path, default=None, help="Optional inversion window summary TSV with H1_inv_start/H1_inv_end.")
    p.add_argument("--output-prefix", type=Path, required=True, help="Output prefix without suffix.")
    p.add_argument("--track-x0", type=int, default=None, help="Optional manual override for top track left pixel.")
    p.add_argument("--track-x1", type=int, default=None, help="Optional manual override for top track right pixel.")
    p.add_argument("--panel-height", type=int, default=330, help="Height of the new expression panel.")
    p.add_argument("--panel-gap", type=int, default=36, help="Gap between the expression panel and the original figure.")
    p.add_argument("--ymax", type=float, default=None, help="Optional fixed y-axis maximum for the expression panel.")
    p.add_argument(
        "--crop-background",
        action="store_true",
        help="Crop whitespace around the original figure before adding the expression panel.",
    )
    p.add_argument(
        "--crop-margin",
        type=int,
        default=24,
        help="Margin in pixels to keep around the cropped background content.",
    )
    p.add_argument("--title", default="H1 RNA expression", help="Panel title.")
    p.add_argument("--bar-color", default="#4C78A8", help="Expression bar fill color.")
    p.add_argument("--bar-outline", default="#274E76", help="Expression bar border color.")
    p.add_argument("--inv-shade-color", default="#F57C00", help="Inversion interval highlight color.")
    return p.parse_args()


def load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    candidates: list[str] = []
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


def normalize_inv_id(token: str) -> str:
    text = token.strip().upper()
    if not text:
        return ""
    if text.startswith("INV"):
        suffix = re.sub(r"^INV[_-]?", "", text)
        if suffix.isdigit():
            return f"INV{int(suffix)}"
        return text
    if text.isdigit():
        return f"INV{int(text)}"
    return text


def trim_transcript(name: str) -> str:
    return name.split(".", 1)[0]


def read_bed(path: Path) -> list[Feature]:
    features: list[Feature] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name, _score, strand = line.rstrip("\n").split("\t")[:6]
            features.append(
                Feature(
                    chrom=chrom,
                    start=int(start),
                    end=int(end),
                    name=name,
                    strand=strand,
                )
            )
    return features


def parse_float(text: str) -> float | None:
    value = text.strip()
    if not value or value in {"NA", "NaN"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_expression(
    path: Path,
    inv_id: str,
    expr_metric: str,
) -> tuple[dict[str, ExpressionGene], int, tuple[int, int] | None, str]:
    genes: dict[str, ExpressionGene] = {}
    inv_interval: tuple[int, int] | None = None
    sample_count = 0
    wanted = normalize_inv_id(inv_id)
    expr_label = "counts"
    with path.open() as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
        reader = csv.DictReader(lines, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Failed to parse header: {path}")
        fieldnames = set(reader.fieldnames)
        if {"INV_ID", "GeneID"}.issubset(fieldnames):
            if expr_metric == "tpm":
                raise ValueError("TPM requires a featureCounts-style table with a Length column.")
            sample_cols = [col for col in reader.fieldnames if col not in META_COLUMNS]
            sample_count = len(sample_cols)
            for row in reader:
                if normalize_inv_id(row.get("INV_ID", "")) != wanted:
                    continue
                gene_id = row["GeneID"].strip()
                values = [v for v in (parse_float(row.get(col, "")) for col in sample_cols) if v is not None]
                if not values:
                    continue
                genes[gene_id] = ExpressionGene(gene_id=gene_id, mean_expr=sum(values) / len(values))
                if inv_interval is None:
                    inv_interval = (int(row["RefStart"]), int(row["RefEnd"]))
        elif "Geneid" in fieldnames:
            sample_cols = [col for col in reader.fieldnames if col not in FEATURECOUNTS_META_COLUMNS]
            sample_count = len(sample_cols)
            rows = [row for row in reader]
            if expr_metric == "tpm":
                expr_label = "TPM"
                parsed_rows: list[tuple[str, float, list[float]]] = []
                rpk_sums = [0.0] * sample_count
                for row in rows:
                    gene_id = row["Geneid"].strip()
                    length_bp = parse_float(row.get("Length", ""))
                    if not gene_id or length_bp is None or length_bp <= 0:
                        continue
                    length_kb = length_bp / 1000.0
                    counts: list[float] = []
                    for col in sample_cols:
                        value = parse_float(row.get(col, ""))
                        counts.append(0.0 if value is None else max(0.0, value))
                    parsed_rows.append((gene_id, length_kb, counts))
                    for i, count in enumerate(counts):
                        rpk_sums[i] += count / length_kb

                scaling_factors = [
                    (rpk_sum / 1_000_000.0) if rpk_sum > 0 else None
                    for rpk_sum in rpk_sums
                ]
                for gene_id, length_kb, counts in parsed_rows:
                    tpms: list[float] = []
                    for i, count in enumerate(counts):
                        factor = scaling_factors[i]
                        if factor is None or factor <= 0:
                            continue
                        tpms.append((count / length_kb) / factor)
                    if not tpms:
                        continue
                    genes[gene_id] = ExpressionGene(gene_id=gene_id, mean_expr=sum(tpms) / len(tpms))
            else:
                expr_label = "counts"
                for row in rows:
                    gene_id = row["Geneid"].strip()
                    values = [v for v in (parse_float(row.get(col, "")) for col in sample_cols) if v is not None]
                    if not values:
                        continue
                    genes[gene_id] = ExpressionGene(gene_id=gene_id, mean_expr=sum(values) / len(values))
        else:
            raise ValueError(f"Unsupported expression table format: {path}")
    return genes, sample_count, inv_interval, expr_label


def load_inv_interval_from_window_info(path: Path, inv_id: str) -> tuple[int, int] | None:
    wanted = normalize_inv_id(inv_id)
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if normalize_inv_id(row.get("ID", "")) != wanted:
                continue
            return int(row["H1_inv_start"]), int(row["H1_inv_end"])
    return None


def rasterize_pdf(pdf_path: Path, density: int, outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    png_path = outdir / f"{pdf_path.stem}.page1.png"
    if shutil.which("convert"):
        cmd = [
            "convert",
            "-density",
            str(density),
            f"{pdf_path}[0]",
            str(png_path),
        ]
    elif shutil.which("gs"):
        cmd = [
            "gs",
            "-q",
            "-dNOPAUSE",
            "-dBATCH",
            "-sDEVICE=pngalpha",
            f"-r{density}",
            f"-sOutputFile={png_path}",
            str(pdf_path),
        ]
    else:
        raise RuntimeError("Neither 'convert' nor 'gs' is available for PDF rasterization.")
    subprocess.run(cmd, check=True)
    return png_path


def load_background(args: argparse.Namespace) -> Image.Image:
    if args.background_png is not None:
        return Image.open(args.background_png).convert("RGBA")
    with tempfile.TemporaryDirectory(prefix="microsynteny_bg_") as tmpdir:
        png_path = rasterize_pdf(args.background_pdf, args.pdf_density, Path(tmpdir))
        return Image.open(png_path).convert("RGBA")


def detect_top_track_span(img: Image.Image) -> tuple[int, int, int]:
    rgb = img.convert("RGB")
    width, height = rgb.size
    candidates: list[tuple[int, int, int, int]] = []
    for y in range(max(0, int(height * 0.02)), min(height, int(height * 0.98))):
        xs: list[int] = []
        for x in range(width):
            r, g, b = rgb.getpixel((x, y))
            if (r + g + b) / 3 < 170:
                xs.append(x)
        if not xs:
            continue
        x0 = min(xs)
        x1 = max(xs)
        count = len(xs)
        if x1 - x0 < int(width * 0.45):
            continue
        candidates.append((count, x0, x1, y))
    if not candidates:
        raise ValueError("Failed to detect the top track span from the background figure.")

    max_span = max(x1 - x0 for _count, x0, x1, _y in candidates)
    near_top_track = [cand for cand in candidates if (cand[2] - cand[1]) >= max_span * 0.9]
    best = min(near_top_track, key=lambda item: item[3])
    return best[1], best[2], best[3]


def crop_background_to_content(
    img: Image.Image,
    bg_color: tuple[int, int, int, int],
    margin: int,
) -> Image.Image:
    bg_rgb = Image.new("RGB", img.size, bg_color[:3])
    diff = ImageChops.difference(img.convert("RGB"), bg_rgb)
    bbox = diff.getbbox()
    if bbox is None:
        return img

    left = max(0, bbox[0] - margin)
    top = max(0, bbox[1] - margin)
    right = min(img.size[0], bbox[2] + margin)
    bottom = min(img.size[1], bbox[3] + margin)
    return img.crop((left, top, right, bottom))


def x_from_bp(bp: float, genome_start: int, genome_end: int, x0: int, x1: int) -> float:
    span = max(1, genome_end - genome_start)
    return x0 + (bp - genome_start) / span * (x1 - x0)


def round_up(value: float) -> int:
    if value <= 0:
        return 1
    magnitude = 10 ** int(math.floor(math.log10(value)))
    step = magnitude
    if value / magnitude < 2:
        step = magnitude / 2
    return int(math.ceil(value / step) * step)


def format_mb(bp: int) -> str:
    return f"{bp / 1_000_000:.2f} Mb"


def text_size(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont) -> tuple[int, int]:
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0], bbox[3] - bbox[1]


def write_summary_table(
    path: Path,
    features: list[Feature],
    expr_map: dict[str, ExpressionGene],
    expr_label: str,
) -> None:
    label = expr_label.strip().lower() if expr_label else "expression"
    value_col = f"mean_{label}"
    with path.open("w") as handle:
        handle.write(f"gene_id\tchrom\tstart\tend\tstrand\t{value_col}\n")
        for feat in features:
            gene_id = trim_transcript(feat.name)
            expr = expr_map.get(gene_id)
            if expr is None:
                continue
            handle.write(
                "\t".join(
                    [
                        gene_id,
                        feat.chrom,
                        str(feat.start),
                        str(feat.end),
                        feat.strand,
                        f"{expr.mean_expr:.6f}",
                    ]
                )
                + "\n"
            )


def draw_expression_panel(
    canvas: Image.Image,
    bg_color: tuple[int, int, int, int],
    features: list[Feature],
    expr_map: dict[str, ExpressionGene],
    inv_interval: tuple[int, int] | None,
    sample_count: int,
    genome_start: int,
    genome_end: int,
    x0: int,
    x1: int,
    panel_height: int,
    ymax_override: float | None,
    expr_label: str,
    title: str,
    bar_color: str,
    bar_outline: str,
    inv_shade_color: str,
) -> None:
    overlay = Image.new("RGBA", canvas.size, (255, 255, 255, 0))
    draw = ImageDraw.Draw(overlay)
    bold = load_font(28, bold=True)
    regular = load_font(20)
    small = load_font(18)

    panel_top = 0
    panel_bottom = panel_height
    plot_top = 92
    plot_bottom = panel_bottom - 52
    plot_height = plot_bottom - plot_top

    draw.rectangle([(0, panel_top), (canvas.size[0], panel_bottom)], fill=bg_color)

    expressed_genes = sum(1 for feat in features if trim_transcript(feat.name) in expr_map)
    subtitle = f"Mean {expr_label} across {sample_count} RNA-seq libraries ({expressed_genes} genes in this H1 window)"
    title_w, title_h = text_size(draw, title, bold)
    subtitle_w, subtitle_h = text_size(draw, subtitle, regular)
    draw.text(((canvas.size[0] - title_w) / 2, 18), title, fill=rgba("#111111"), font=bold)
    draw.text(((canvas.size[0] - subtitle_w) / 2, 18 + title_h + 6), subtitle, fill=rgba("#555555"), font=regular)

    values = [
        expr_map[trim_transcript(feat.name)].mean_expr
        for feat in features
        if trim_transcript(feat.name) in expr_map
    ]
    if ymax_override is not None:
        ymax = max(1.0, float(ymax_override))
    else:
        ymax = round_up(max(values)) if values else 1

    if inv_interval is not None:
        inv_x0 = x_from_bp(inv_interval[0], genome_start, genome_end, x0, x1)
        inv_x1 = x_from_bp(inv_interval[1], genome_start, genome_end, x0, x1)
        draw.rectangle(
            [(inv_x0, plot_top), (inv_x1, plot_bottom)],
            fill=rgba(inv_shade_color, 35),
        )
        draw.text((inv_x0 + 6, plot_top + 6), "Inversion interval", fill=rgba(inv_shade_color), font=small)

    tick_values = [0, ymax / 2, ymax]
    for tick in tick_values:
        y = plot_bottom - (tick / ymax) * plot_height
        draw.line([(x0, y), (x1, y)], fill=rgba("#D5D5D5"), width=1)
        label = f"{int(tick):,}"
        label_w, label_h = text_size(draw, label, small)
        draw.text((x0 - label_w - 10, y - label_h / 2), label, fill=rgba("#555555"), font=small)

    draw.line([(x0, plot_bottom), (x1, plot_bottom)], fill=rgba("#666666"), width=2)
    draw.line([(x0, plot_top), (x0, plot_bottom)], fill=rgba("#666666"), width=2)

    for feat in features:
        expr = expr_map.get(trim_transcript(feat.name))
        if expr is None:
            continue
        bar_x0 = x_from_bp(feat.start, genome_start, genome_end, x0, x1)
        bar_x1 = x_from_bp(feat.end, genome_start, genome_end, x0, x1)
        if bar_x1 <= bar_x0:
            bar_x1 = bar_x0 + 1
        scaled_expr = min(expr.mean_expr, ymax)
        bar_top = plot_bottom - (scaled_expr / ymax) * plot_height
        bar_top = max(plot_top, bar_top)
        draw.rectangle(
            [(bar_x0, bar_top), (bar_x1, plot_bottom)],
            fill=rgba(bar_color, 210),
            outline=rgba(bar_outline),
            width=1,
        )

    start_label = format_mb(genome_start)
    end_label = format_mb(genome_end)
    start_w, start_h = text_size(draw, start_label, small)
    end_w, end_h = text_size(draw, end_label, small)
    draw.text((x0, plot_bottom + 10), start_label, fill=rgba("#555555"), font=small)
    draw.text((x1 - end_w, plot_bottom + 10), end_label, fill=rgba("#555555"), font=small)
    axis_label = f"{features[0].chrom} coordinates aligned to the original H1 track"
    axis_w, _axis_h = text_size(draw, axis_label, small)
    draw.text(((x0 + x1 - axis_w) / 2, plot_bottom + 10), axis_label, fill=rgba("#555555"), font=small)

    canvas.alpha_composite(overlay)


def main() -> None:
    args = parse_args()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    background = load_background(args)
    if args.crop_background:
        background = crop_background_to_content(
            background,
            background.getpixel((0, 0)),
            args.crop_margin,
        )
    h1_features = read_bed(args.h1_bed)
    expr_map, sample_count, inv_interval, expr_label = load_expression(
        args.expr_table,
        args.inv_id,
        args.expr_metric,
    )
    if inv_interval is None and args.window_info is not None:
        inv_interval = load_inv_interval_from_window_info(args.window_info, args.inv_id)

    if not h1_features:
        raise ValueError(f"No H1 features found in {args.h1_bed}")
    if not expr_map:
        raise ValueError(f"No expression rows found for {args.inv_id} in {args.expr_table}")

    if args.track_x0 is not None and args.track_x1 is not None:
        track_x0, track_x1 = args.track_x0, args.track_x1
    else:
        track_x0, track_x1, _track_y = detect_top_track_span(background)

    genome_start = min(feat.start for feat in h1_features)
    genome_end = max(feat.end for feat in h1_features)
    panel_offset = args.panel_height + args.panel_gap
    bg_color = background.getpixel((0, 0))

    canvas = Image.new(
        "RGBA",
        (background.size[0], background.size[1] + panel_offset),
        bg_color,
    )
    draw_expression_panel(
        canvas=canvas,
        bg_color=bg_color,
        features=h1_features,
        expr_map=expr_map,
        inv_interval=inv_interval,
        sample_count=sample_count,
        genome_start=genome_start,
        genome_end=genome_end,
        x0=track_x0,
        x1=track_x1,
        panel_height=args.panel_height,
        ymax_override=args.ymax,
        expr_label=expr_label,
        title=args.title,
        bar_color=args.bar_color,
        bar_outline=args.bar_outline,
        inv_shade_color=args.inv_shade_color,
    )
    canvas.alpha_composite(background, dest=(0, panel_offset))

    png_path = Path(f"{args.output_prefix}.png")
    pdf_path = Path(f"{args.output_prefix}.pdf")
    table_path = Path(f"{args.output_prefix}.expression.tsv")

    canvas.save(png_path)
    canvas.convert("RGB").save(pdf_path)
    write_summary_table(table_path, h1_features, expr_map, expr_label)


if __name__ == "__main__":
    main()
