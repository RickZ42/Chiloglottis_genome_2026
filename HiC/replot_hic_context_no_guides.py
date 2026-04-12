#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fmt_mb(bp: float) -> str:
    return f"{bp / 1e6:.2f}".rstrip("0").rstrip(".")


def dense_from_dump(path: Path, h1_start: int, h1_end: int, h2_start: int, h2_end: int, bin_size: int) -> np.ndarray:
    h1_bins = int(math.ceil((h1_end - h1_start + 1) / float(bin_size)))
    h2_bins = int(math.ceil((h2_end - h2_start + 1) / float(bin_size)))
    mat = np.zeros((h1_bins, h2_bins), dtype=np.float64)
    if not path.exists() or path.stat().st_size == 0:
        return mat
    with path.open() as f:
        for ln in f:
            if not ln.strip():
                continue
            c = ln.split()
            if len(c) < 3:
                continue
            try:
                x = int(float(c[0]))
                y = int(float(c[1]))
                v = float(c[2])
            except ValueError:
                continue
            if not math.isfinite(v):
                continue
            i = int((x - h1_start) // bin_size)
            j = int((y - h2_start) // bin_size)
            if i < 0 or j < 0 or i >= h1_bins or j >= h2_bins:
                continue
            mat[i, j] += v
    return mat


def main() -> None:
    ap = argparse.ArgumentParser(description="Replot a Hi-C context panel without exact inversion guide lines.")
    ap.add_argument("--summary-tsv", type=Path, required=True)
    ap.add_argument("--inv-id", required=True)
    ap.add_argument("--out-prefix", type=Path, required=True)
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    row = None
    with args.summary_tsv.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["ID"] == args.inv_id:
                row = r
                break
    if row is None:
        raise SystemExit(f"{args.inv_id} not found in {args.summary_tsv}")

    h1_start = int(row["H1_window_start"])
    h1_end = int(row["H1_window_end"])
    h2_start = int(row["H2_window_start"])
    h2_end = int(row["H2_window_end"])
    bin_size = int(row["bin_size"])
    dump_tsv = Path(row["Inter_dump"])

    mat = dense_from_dump(dump_tsv, h1_start, h1_end, h2_start, h2_end, bin_size)
    vmax = np.percentile(np.log1p(mat[mat > 0]), 99.5) if np.any(mat > 0) else 1.0

    fig, ax = plt.subplots(figsize=(7.4, 6.8), dpi=300, facecolor="white")
    ax.set_facecolor("#fbfbfb")
    im = ax.imshow(
        np.log1p(mat),
        origin="lower",
        cmap="Reds",
        interpolation="nearest",
        aspect="equal",
        vmin=0.0,
        vmax=vmax,
    )

    nbin_y, nbin_x = mat.shape
    xticks = np.linspace(0, max(0, nbin_x - 1), 5)
    yticks = np.linspace(0, max(0, nbin_y - 1), 5)
    ax.set_xticks(xticks)
    ax.set_xticklabels([fmt_mb(x) for x in np.linspace(h2_start, h2_end, 5)], fontsize=10, fontweight="bold")
    ax.set_yticks(yticks)
    ax.set_yticklabels([fmt_mb(y) for y in np.linspace(h1_start, h1_end, 5)], fontsize=10, fontweight="bold")
    ax.set_xlabel("H2 position (Mb)", fontsize=12, fontweight="bold")
    ax.set_ylabel("H1 position (Mb)", fontsize=12, fontweight="bold")
    title = args.title or f"{args.inv_id} | H1 {row['RefChr']} vs H2 {row['QryChr']} | {row['Norm']}-normalized"
    ax.set_title(title, fontsize=14, fontweight="bold", pad=10)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("log1p normalized contact", fontsize=11, fontweight="bold")
    cbar.ax.tick_params(labelsize=9, width=1.0, length=4)

    ax.text(
        0.01,
        -0.10,
        "Guide lines omitted because the inversion is near the effective Hi-C bin limit; synteny/gbrowse define exact boundaries.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        color="#444444",
    )

    out_prefix = args.out_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".png"), facecolor="white", bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")


if __name__ == "__main__":
    main()
