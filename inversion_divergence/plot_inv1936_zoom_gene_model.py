#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


OUTDIR = Path("/g/data/xf3/zz3507/Output/20260127Genome/inv1936_truncation_analysis")

W = 2200
H = 1100
MARGIN_X = 140
TRACK_W = 1880


def iv(start: int, end: int) -> tuple[int, int]:
    return (start, end)


H1 = {
    "title": "H1 g13041.t1",
    "subtitle": "Breakpoint-overlapping gene on scaffold_5; repeat context: (AT)n repeat",
    "gene_start": 7223150,
    "gene_end": 7278940,
    "breakpoint": 7277582,
    "segment": iv(7277582, 7278940),
    "repeat_spans": [iv(7277542, 7277606)],
    "note": "Breakpoint in intron 7263637-7278117; 536 bp before the next exon",
    "exons": [
        iv(7223150, 7223315), iv(7230042, 7230142), iv(7240967, 7241065),
        iv(7262738, 7262918), iv(7263002, 7263050), iv(7263141, 7263219),
        iv(7263389, 7263463), iv(7263559, 7263636), iv(7278118, 7278276),
        iv(7278571, 7278666), iv(7278752, 7278940),
    ],
    "cds": [
        iv(7223150, 7223315), iv(7230042, 7230142), iv(7240967, 7241065),
        iv(7262738, 7262918), iv(7263002, 7263050), iv(7263141, 7263219),
        iv(7263389, 7263463), iv(7263559, 7263636), iv(7278118, 7278276),
        iv(7278571, 7278666), iv(7278752, 7278940),
    ],
}

H2 = {
    "title": "H2 g13243.t1",
    "subtitle": "Anchored mirrored counterpart on scaffold_5; repeat context: repeat-family hits",
    "gene_start": 7448673,
    "gene_end": 7501391,
    "breakpoint": 7450031,
    "segment": iv(7448673, 7450031),
    "repeat_spans": [iv(7450008, 7450166), iv(7450017, 7450226)],
    "note": "Breakpoint in intron 7449496-7463092; 536 bp after the previous exon",
    "exons": [
        iv(7448673, 7448861), iv(7448947, 7449042), iv(7449337, 7449495),
        iv(7463093, 7463170), iv(7463266, 7463340), iv(7463510, 7463588),
        iv(7463679, 7463727), iv(7463811, 7463991), iv(7483513, 7483611),
        iv(7494444, 7494544), iv(7501226, 7501391),
    ],
    "cds": [
        iv(7448673, 7448861), iv(7448947, 7449042), iv(7449337, 7449495),
        iv(7463093, 7463170), iv(7463266, 7463340), iv(7463510, 7463588),
        iv(7463679, 7463727), iv(7463811, 7463991), iv(7483513, 7483611),
        iv(7494444, 7494544), iv(7501226, 7501391),
    ],
}


def font(size: int) -> ImageFont.ImageFont:
    try:
        return ImageFont.truetype("DejaVuSans.ttf", size)
    except OSError:
        return ImageFont.load_default()


FONT_TITLE = font(42)
FONT_SUB = font(26)
FONT_LABEL = font(24)
FONT_SMALL = font(22)


def xmap(pos: int, gene_start: int, gene_end: int) -> int:
    return int(MARGIN_X + (pos - gene_start) * TRACK_W / (gene_end - gene_start))


def draw_locus(draw: ImageDraw.ImageDraw, top: int, locus: dict[str, object]) -> None:
    gene_start = int(locus["gene_start"])
    gene_end = int(locus["gene_end"])
    breakpoint = int(locus["breakpoint"])
    seg_start, seg_end = locus["segment"]

    title_y = top
    subtitle_y = top + 50
    track_y = top + 140
    exon_y = top + 120
    cds_y = top + 135
    note_y = top + 260

    draw.text((MARGIN_X, title_y), locus["title"], fill="#111111", font=FONT_TITLE)
    draw.text((MARGIN_X, subtitle_y), locus["subtitle"], fill="#4d4d4d", font=FONT_SUB)

    for start, end in locus["repeat_spans"]:
        draw.rectangle([xmap(start, gene_start, gene_end), top + 88, xmap(end, gene_start, gene_end), top + 248], fill="#f6b0a7")

    draw.rectangle([xmap(seg_start, gene_start, gene_end), top + 96, xmap(seg_end, gene_start, gene_end), top + 240], fill="#fdd49e")
    draw.line([(MARGIN_X, track_y), (MARGIN_X + TRACK_W, track_y)], fill="#666666", width=5)

    for start, end in locus["exons"]:
        draw.rectangle([xmap(start, gene_start, gene_end), exon_y, xmap(end, gene_start, gene_end), exon_y + 36], fill="#cfcfcf", outline="#6b6b6b", width=2)
    for start, end in locus["cds"]:
        draw.rectangle([xmap(start, gene_start, gene_end), cds_y, xmap(end, gene_start, gene_end), cds_y + 18], fill="#4c78a8", outline="#274e76", width=2)

    bp_x = xmap(breakpoint, gene_start, gene_end)
    draw.line([(bp_x, top + 70), (bp_x, top + 250)], fill="#e66101", width=4)
    draw.text((bp_x - 68, top + 18), f"Breakpoint\n{breakpoint:,}", fill="#b35806", font=FONT_LABEL)

    draw.text((MARGIN_X + 580, note_y), "Boundary-separated block: 1,359 bp, 3 exons, 441 coding bp (~147 aa)", fill="#7f3b08", font=FONT_LABEL)
    draw.text((MARGIN_X, note_y), locus["note"], fill="#444444", font=FONT_LABEL)

    axis_y = top + 320
    draw.line([(MARGIN_X, axis_y), (MARGIN_X + TRACK_W, axis_y)], fill="#555555", width=2)
    for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
        x = int(MARGIN_X + TRACK_W * frac)
        pos = int(gene_start + (gene_end - gene_start) * frac)
        draw.line([(x, axis_y), (x, axis_y + 12)], fill="#555555", width=2)
        draw.text((x - 34, axis_y + 16), f"{pos / 1_000_000:.3f}", fill="#555555", font=FONT_SMALL)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    image = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(image)

    draw.text((MARGIN_X, 30), "INV1936 zoom-in gene model", fill="#111111", font=FONT_TITLE)
    draw.text(
        (MARGIN_X, 86),
        "Both haplotypes show an intronic breakpoint that separates the same terminal 3-exon block across the inversion boundary.",
        fill="#444444",
        font=FONT_SUB,
    )

    draw_locus(draw, 170, H1)
    draw_locus(draw, 610, H2)

    png = OUTDIR / "INV1936_zoom_gene_model.png"
    pdf = OUTDIR / "INV1936_zoom_gene_model.pdf"
    image.save(png, dpi=(300, 300))
    image.save(pdf, resolution=300)


if __name__ == "__main__":
    main()
