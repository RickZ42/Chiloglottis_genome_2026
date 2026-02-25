#!/usr/bin/env python3
"""
Plot local gene-order anchor dotplot for an inversion candidate.

This keeps the original gene order on both haplotypes, so an inversion appears
as a negative-slope trend (top-left to bottom-right) across matched anchors.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_bed_order(bedfile: Path):
    genes = []
    with bedfile.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, start, end, gene = parts[:4]
            genes.append((gene, chrom, int(start), int(end)))
    gene_to_idx = {g: i + 1 for i, (g, *_rest) in enumerate(genes)}
    return genes, gene_to_idx


def read_blocks_pairs(blocksfile: Path):
    pairs = []
    with blocksfile.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            g1, g2 = parts[:2]
            if g2 == ".":
                continue
            pairs.append((g1, g2))
    return pairs


def parse_window_info(infofile: Path):
    if not infofile.exists():
        return {}
    with infofile.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            return row
    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h1_bed", type=Path)
    ap.add_argument("h2_bed", type=Path)
    ap.add_argument("blocks", type=Path)
    ap.add_argument("--window-info", type=Path, default=None)
    ap.add_argument("--title", default="Local gene-order anchor dotplot")
    ap.add_argument("--outprefix", type=Path, required=True)
    args = ap.parse_args()

    h1_genes, h1_idx = read_bed_order(args.h1_bed)
    h2_genes, h2_idx = read_bed_order(args.h2_bed)
    pairs = read_blocks_pairs(args.blocks)

    pts = []
    for g1, g2 in pairs:
        if g1 in h1_idx and g2 in h2_idx:
            pts.append((h1_idx[g1], h2_idx[g2], g1, g2))
    pts.sort(key=lambda x: x[0])

    if not pts:
        raise SystemExit("No matched anchor pairs found in blocks file.")

    window_info = parse_window_info(args.window_info) if args.window_info else {}

    fig, ax = plt.subplots(figsize=(7.2, 6.2), dpi=220)

    # Reference inversion trend line (negative slope) for visual cue.
    ax.plot(
        [1, len(h1_genes)],
        [len(h2_genes), 1],
        linestyle="--",
        linewidth=1.2,
        color="#BBBBBB",
        zorder=1,
        label="Expected inversion trend",
    )

    xs = [x for x, *_ in pts]
    ys = [y for _, y, *_ in pts]

    # Connect anchors in H1 order to make monotonic trend obvious.
    ax.plot(xs, ys, color="#D62728", linewidth=1.4, alpha=0.8, zorder=2)
    ax.scatter(xs, ys, s=28, color="#1F77B4", edgecolor="white", linewidth=0.5, zorder=3)

    # Light labels for endpoints.
    first = pts[0]
    last = pts[-1]
    for x, y, g1, g2, ha, va in [
        (*first, "left", "bottom"),
        (*last, "right", "top"),
    ]:
        ax.text(
            x,
            y,
            f"{g1} ↔ {g2}",
            fontsize=7,
            ha=ha,
            va=va,
            color="#333333",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.75),
            zorder=4,
        )

    ax.set_xlim(0.5, len(h1_genes) + 0.5)
    ax.set_ylim(0.5, len(h2_genes) + 0.5)
    ax.invert_yaxis()  # top-left to bottom-right looks like inversion trend
    ax.set_xlabel(f"H1 local gene order (n={len(h1_genes)})")
    ax.set_ylabel(f"H2 local gene order (n={len(h2_genes)})")
    ax.grid(True, linestyle=":", linewidth=0.6, color="#DDDDDD")
    ax.set_axisbelow(True)

    title = args.title
    subtitle = None
    if window_info:
        inv_id = window_info.get("ID", "")
        h1_chr = window_info.get("H1_chr", "")
        h1s = window_info.get("H1_inv_start", "")
        h1e = window_info.get("H1_inv_end", "")
        h2_chr = window_info.get("H2_chr", "")
        h2s = window_info.get("H2_inv_start", "")
        h2e = window_info.get("H2_inv_end", "")
        subtitle = (
            f"{inv_id}: H1 {h1_chr}:{h1s}-{h1e} vs H2 {h2_chr}:{h2s}-{h2e} "
            f"(matched anchors={len(pts)})"
        )
    ax.set_title(title, fontsize=12, weight="bold", pad=10)
    if subtitle:
        fig.text(0.5, 0.955, subtitle, ha="center", va="center", fontsize=8, color="#444444")

    ax.legend(loc="lower left", fontsize=7, frameon=False)

    fig.tight_layout(rect=(0, 0, 1, 0.94 if subtitle else 1))

    outprefix = args.outprefix
    outprefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(outprefix.with_suffix(f".{ext}"))


if __name__ == "__main__":
    main()

