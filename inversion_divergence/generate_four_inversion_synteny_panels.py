#!/usr/bin/env python3
"""
Generate unified local H1/H2 gene-order synteny dotplots for selected inversions.

This uses:
- SyRI inversion coordinates from the breakpoint priority screen table
- Global H1/H2 `.t1` BED gene models
- Global H1/H2 lifted anchors

For each inversion, genes within +/- flank_bp around the inversion interval are
extracted on both haplotypes, lifted anchors are filtered to those local genes,
and a local gene-order anchor dotplot is rendered using the existing plotting
script.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


ROOT = Path("/g/data/xf3/zz3507")
DEFAULT_INV_TABLE = ROOT / "Output/20260127Genome/other_syri_inversion/priority_breakpoint_screen/all_inversion_breakpoint_priority_screen.tsv"
DEFAULT_H1_BED = ROOT / "compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H1t1.bed"
DEFAULT_H2_BED = ROOT / "compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H2t1.bed"
DEFAULT_ANCHORS = ROOT / "compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H1t1.H2t1.lifted.anchors"
DEFAULT_PLOT_SCRIPT = ROOT / "script/Chiloglottis_genome_2026/syri/plot_local_gene_order_anchor_dotplot.py"


def read_inv_rows(path: Path, wanted: set[str]) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            inv_id = row["ID"]
            if inv_id in wanted:
                rows[inv_id] = row
    return rows


def read_bed(path: Path) -> list[tuple[str, int, int, str, str]]:
    recs: list[tuple[str, int, int, str, str]] = []
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, gene, _score, strand = line.rstrip("\n").split("\t")[:6]
            recs.append((chrom, int(start), int(end), gene, strand))
    return recs


def read_anchor_pairs(path: Path) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            a, b, *_rest = line.rstrip("\n").split("\t")
            pairs.append((a, b))
    return pairs


def subset_bed(
    bed: list[tuple[str, int, int, str, str]],
    chrom: str,
    start: int,
    end: int,
) -> list[tuple[str, int, int, str, str]]:
    return [r for r in bed if r[0] == chrom and not (r[2] < start or r[1] > end)]


def write_bed(path: Path, recs: list[tuple[str, int, int, str, str]]) -> None:
    with path.open("w") as f:
        for chrom, start, end, gene, strand in recs:
            f.write(f"{chrom}\t{start}\t{end}\t{gene}\t0\t{strand}\n")


def write_blocks(path: Path, pairs: list[tuple[str, str]]) -> None:
    with path.open("w") as f:
        for a, b in pairs:
            f.write(f"{a}\t{b}\n")


def write_window_info(path: Path, row: dict[str, str]) -> None:
    fieldnames = [
        "ID",
        "H1_chr",
        "H1_inv_start",
        "H1_inv_end",
        "H2_chr",
        "H2_inv_start",
        "H2_inv_end",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        w.writerow(
            {
                "ID": row["ID"],
                "H1_chr": row["RefChr"],
                "H1_inv_start": row["RefStart"],
                "H1_inv_end": row["RefEnd"],
                "H2_chr": row["QryChr"],
                "H2_inv_start": row["QryStart"],
                "H2_inv_end": row["QryEnd"],
            }
        )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inv-ids", nargs="+", required=True)
    ap.add_argument("--flank-bp", type=int, default=2_000_000)
    ap.add_argument("--inv-table", type=Path, default=DEFAULT_INV_TABLE)
    ap.add_argument("--h1-bed", type=Path, default=DEFAULT_H1_BED)
    ap.add_argument("--h2-bed", type=Path, default=DEFAULT_H2_BED)
    ap.add_argument("--anchors", type=Path, default=DEFAULT_ANCHORS)
    ap.add_argument("--plot-script", type=Path, default=DEFAULT_PLOT_SCRIPT)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    wanted = set(args.inv_ids)
    inv_rows = read_inv_rows(args.inv_table, wanted)
    missing = sorted(wanted - set(inv_rows))
    if missing:
        raise SystemExit(f"Missing inversion rows: {', '.join(missing)}")

    h1_bed = read_bed(args.h1_bed)
    h2_bed = read_bed(args.h2_bed)
    all_pairs = read_anchor_pairs(args.anchors)

    outdir = args.outdir
    data_dir = outdir / "data"
    png_dir = outdir / "png"
    pdf_dir = outdir / "pdf"
    data_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(parents=True, exist_ok=True)
    pdf_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, str | int]] = []

    for inv_id in args.inv_ids:
        row = inv_rows[inv_id]
        ref_chr = row["RefChr"]
        ref_start = int(row["RefStart"])
        ref_end = int(row["RefEnd"])
        qry_chr = row["QryChr"]
        qry_start = int(row["QryStart"])
        qry_end = int(row["QryEnd"])

        h1_start = max(0, ref_start - args.flank_bp)
        h1_end = ref_end + args.flank_bp
        h2_start = max(0, qry_start - args.flank_bp)
        h2_end = qry_end + args.flank_bp

        local_h1 = subset_bed(h1_bed, ref_chr, h1_start, h1_end)
        local_h2 = subset_bed(h2_bed, qry_chr, h2_start, h2_end)
        h1_genes = {g for *_rest, g, _strand in [(c, s, e, g, st) for c, s, e, g, st in local_h1]}
        h2_genes = {g for *_rest, g, _strand in [(c, s, e, g, st) for c, s, e, g, st in local_h2]}
        local_pairs = [(a, b) for a, b in all_pairs if a in h1_genes and b in h2_genes]

        inv_dir = data_dir / inv_id
        inv_dir.mkdir(parents=True, exist_ok=True)
        h1_local_bed = inv_dir / f"H1_{inv_id}_flank{args.flank_bp}.t1.bed"
        h2_local_bed = inv_dir / f"H2_{inv_id}_flank{args.flank_bp}.t1.bed"
        blocks = inv_dir / f"H1_H2_{inv_id}_flank{args.flank_bp}.blocks"
        window_info = inv_dir / f"{inv_id}_flank{args.flank_bp}_window_info.tsv"
        outprefix = inv_dir / f"H1_vs_H2_{inv_id}_geneorder_anchor_dotplot_flank{args.flank_bp}"

        write_bed(h1_local_bed, local_h1)
        write_bed(h2_local_bed, local_h2)
        write_blocks(blocks, local_pairs)
        write_window_info(window_info, row)

        title = f"{inv_id} local gene-order synteny dotplot (H1 vs H2)"
        subprocess.run(
            [
                "python",
                str(args.plot_script),
                str(h1_local_bed),
                str(h2_local_bed),
                str(blocks),
                "--window-info",
                str(window_info),
                "--title",
                title,
                "--outprefix",
                str(outprefix),
            ],
            check=True,
        )

        png_path = outprefix.with_suffix(".png")
        pdf_path = outprefix.with_suffix(".pdf")
        png_link = png_dir / f"{inv_id}.geneorder_synteny_flank{args.flank_bp}.png"
        pdf_link = pdf_dir / f"{inv_id}.geneorder_synteny_flank{args.flank_bp}.pdf"
        if png_link.exists() or png_link.is_symlink():
            png_link.unlink()
        if pdf_link.exists() or pdf_link.is_symlink():
            pdf_link.unlink()
        png_link.symlink_to(png_path)
        pdf_link.symlink_to(pdf_path)

        manifest_rows.append(
            {
                "ID": inv_id,
                "RefChr": ref_chr,
                "RefStart": ref_start,
                "RefEnd": ref_end,
                "QryChr": qry_chr,
                "QryStart": qry_start,
                "QryEnd": qry_end,
                "flank_bp": args.flank_bp,
                "h1_gene_count": len(local_h1),
                "h2_gene_count": len(local_h2),
                "anchor_pair_count": len(local_pairs),
                "png": str(png_link),
                "pdf": str(pdf_link),
            }
        )

    manifest_path = outdir / "manifest.tsv"
    with manifest_path.open("w", newline="") as f:
        fieldnames = [
            "ID",
            "RefChr",
            "RefStart",
            "RefEnd",
            "QryChr",
            "QryStart",
            "QryEnd",
            "flank_bp",
            "h1_gene_count",
            "h2_gene_count",
            "anchor_pair_count",
            "png",
            "pdf",
        ]
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        w.writerows(manifest_rows)


if __name__ == "__main__":
    main()
