#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import shutil
import subprocess
from pathlib import Path


ROOT = Path("/g/data/xf3/zz3507")
OUTDIR = ROOT / "Output/20260127Genome/inv1936_H1_H2_IGV_comparison_bundle"
H1_DIR = OUTDIR / "H1"
H2_DIR = OUTDIR / "H2"

SAMTOOLS = Path("/g/data/xf3/miniconda/envs/common-tools/bin/samtools")

H1_REF = ROOT / "Output/20260127Genome/H1/H1_20260127.FINAL.top20.fa"
H1_GTF = ROOT / "Output/20260127Genome/H1/breaker/braker.gtf"
H1_GFF = ROOT / "Output/20260127Genome/H1/RepeatMasker/H1_20260127.FINAL.top20.fa.out.gff"
H1_LOCAL_BED = ROOT / "Output/20260127Genome/syri/syri_asm10/INV1936_jcvi_microsynteny/H1_INV1936_local.t1.bed"
H1_BP_BED = ROOT / "Output/20260127Genome/inv1936_truncation_analysis/INV1936_H1_breakpoint.bed"
H1_GENE_BED = ROOT / "Output/20260127Genome/inv1936_truncation_analysis/INV1936_H1_gene.bed"
H1_BAMS = [
    ROOT / "Output/20260127Genome/H1/inversion_rna_minimap2/bam/141L_S13.sorted.bam",
    ROOT / "Output/20260127Genome/H1/inversion_rna_minimap2/bam/142L_S14.sorted.bam",
    ROOT / "Output/20260127Genome/H1/inversion_rna_minimap2/bam/143L_S15.sorted.bam",
]

H2_REF = ROOT / "Output/20260127Genome/H2/H2_20260127.FINAL.top20.ordered.renamed.fa"
H2_GTF = ROOT / "Output/20260127Genome/H2/breaker_H2_2026_0211/braker.gtf"
H2_GFF = ROOT / "Output/20260127Genome/H2/RepeatMasker/H2_20260127.FINAL.top20.ordered.renamed.fa.out.gff"
H2_LOCAL_BED = ROOT / "Output/20260127Genome/syri/syri_asm10/INV1936_jcvi_microsynteny/H2_INV1936_local.t1.bed"
H2_BP_BED = ROOT / "Output/20260127Genome/inv1936_truncation_analysis/INV1936_H2_breakpoint.bed"
H2_GENE_BED = ROOT / "Output/20260127Genome/inv1936_truncation_analysis/INV1936_H2_gene.bed"
H2_BAMS = [
    ROOT / "Output/20260127Genome/inv1936_targeted_h2_remap/bam/141L_S13.H2_targeted.sorted.bam",
    ROOT / "Output/20260127Genome/inv1936_targeted_h2_remap/bam/142L_S14.H2_targeted.sorted.bam",
    ROOT / "Output/20260127Genome/inv1936_targeted_h2_remap/bam/143L_S15.H2_targeted.sorted.bam",
]

SCAFFOLD = "scaffold_5"

H1_GENE_START = 7_223_150
H1_GENE_END = 7_278_940
H1_BP = 7_277_582
H1_INV_START = 7_277_582
H1_INV_END = 8_517_764
H1_BROAD_START = 7_213_150
H1_BROAD_END = 7_288_940
H1_ZOOM_START = 7_276_382
H1_ZOOM_END = 7_279_782

H2_GENE_START = 7_448_673
H2_GENE_END = 7_501_391
H2_BP = 7_450_031
H2_INV_START = 6_399_204
H2_INV_END = 7_450_031
H2_BROAD_START = 7_438_673
H2_BROAD_END = 7_511_391
H2_ZOOM_START = 7_448_831
H2_ZOOM_END = 7_452_231


def ensure_dirs() -> None:
    H1_DIR.mkdir(parents=True, exist_ok=True)
    H2_DIR.mkdir(parents=True, exist_ok=True)


def extract_scaffold(ref_fa: Path, out_fa: Path, scaffold: str = SCAFFOLD) -> None:
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    found = False
    with ref_fa.open() as src, out_fa.open("w") as dst:
        keep = False
        for line in src:
            if line.startswith(">"):
                name = line[1:].split()[0]
                keep = name == scaffold
                if keep:
                    found = True
                    dst.write(line)
                continue
            if keep:
                dst.write(line)
    if not found:
        raise ValueError(f"{scaffold} not found in {ref_fa}")
    subprocess.run([str(SAMTOOLS), "faidx", str(out_fa)], check=True)


def _iter_lines(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            yield line


def overlap(start0: int, end0: int, qstart0: int, qend0: int) -> bool:
    return start0 < qend0 and end0 > qstart0


def subset_gtf(in_path: Path, out_path: Path, chrom: str, start1: int, end1: int) -> None:
    qstart0 = start1 - 1
    qend0 = end1
    with out_path.open("w") as out:
        for line in _iter_lines(in_path):
            if not line.strip() or line.startswith("#"):
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != chrom:
                continue
            start0 = int(parts[3]) - 1
            end0 = int(parts[4])
            if overlap(start0, end0, qstart0, qend0):
                out.write(line)


def subset_gff(in_path: Path, out_path: Path, chrom: str, start1: int, end1: int) -> None:
    qstart0 = start1 - 1
    qend0 = end1
    with out_path.open("w") as out:
        for line in _iter_lines(in_path):
            if not line.strip() or line.startswith("#"):
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != chrom:
                continue
            start0 = int(parts[3]) - 1
            end0 = int(parts[4])
            if overlap(start0, end0, qstart0, qend0):
                out.write(line)


def subset_bed(in_path: Path, out_path: Path, chrom: str, start1: int, end1: int) -> None:
    qstart0 = start1 - 1
    qend0 = end1
    with in_path.open() as src, out_path.open("w") as out:
        for line in src:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != chrom:
                continue
            start0 = int(parts[1])
            end0 = int(parts[2])
            if overlap(start0, end0, qstart0, qend0):
                out.write(line)


def write_interval_beds(
    outdir: Path,
    prefix: str,
    chrom: str,
    inv_start1: int,
    inv_end1: int,
    gene_start1: int,
    gene_end1: int,
    gene_label: str,
    gene_strand: str,
) -> None:
    interval = outdir / f"{prefix}_interval.bed"
    boundaries = outdir / f"{prefix}_boundaries.bed"
    gene = outdir / f"{prefix}_gene.bed"
    with interval.open("w") as fh:
        fh.write(f"{chrom}\t{inv_start1 - 1}\t{inv_end1}\t{prefix}_interval\t0\t+\n")
    with boundaries.open("w") as fh:
        fh.write(f"{chrom}\t{inv_start1 - 1}\t{inv_start1}\t{prefix}_left_boundary\t0\t+\n")
        fh.write(f"{chrom}\t{inv_end1 - 1}\t{inv_end1}\t{prefix}_right_boundary\t0\t+\n")
    with gene.open("w") as fh:
        fh.write(f"{chrom}\t{gene_start1 - 1}\t{gene_end1}\t{gene_label}\t0\t{gene_strand}\n")


def slice_bam(in_bam: Path, out_bam: Path, region: str) -> None:
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([str(SAMTOOLS), "view", "-b", "-o", str(out_bam), str(in_bam), region], check=True)
    subprocess.run([str(SAMTOOLS), "index", str(out_bam)], check=True)


def write_manifest(out_path: Path) -> None:
    rows = [
        ["haplotype", "sample", "bam"],
    ]
    for bam in H1_BAMS:
        rows.append(["H1", bam.name.replace(".sorted.bam", ""), str(bam)])
    for bam in H2_BAMS:
        rows.append(["H2", bam.name.replace(".sorted.bam", ""), str(bam)])
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerows(rows)


def write_readme(out_path: Path) -> None:
    text = f"""INV1936 H1/H2 IGV-ready comparison bundle

Purpose
- Compare the same sunflower-stage labellum RNA subset across H1 and H2 for INV1936.
- Support side-by-side inspection of gene structure, breakpoint position, repeat context, and RNA read/junction behavior.

Mapping note
- H1 BAMs are slices from the existing full-genome H1 minimap2 RNA alignments.
- H2 BAMs are derived from IGV-focused remapping of the same three libraries to the H2 scaffold_5 local reference, used here to enable direct read-level inspection of the focal locus under local resource limits.

Samples
- 141L_S13
- 142L_S14
- 143L_S15

Bundle layout
- H1/
  - H1_20260127.FINAL.scaffold_5.fa (+ .fai)
  - braker.INV1936_H1_broad_region.gtf
  - H1_20260127.FINAL.top20.fa.out.INV1936_H1_broad_region.gff
  - INV1936_H1_breakpoint.bed
  - INV1936_H1_gene.bed
  - INV1936_H1_interval.bed
  - INV1936_H1_boundaries.bed
  - H1_INV1936_local.t1.bed
  - H1 broad-region BAM slices for 141L_S13, 142L_S14, 143L_S15
- H2/
  - H2_20260127.FINAL.scaffold_5.fa (+ .fai)
  - braker.INV1936_H2_broad_region.gtf
  - H2_20260127.FINAL.top20.ordered.renamed.fa.out.INV1936_H2_broad_region.gff
  - INV1936_H2_breakpoint.bed
  - INV1936_H2_gene.bed
  - INV1936_H2_interval.bed
  - INV1936_H2_boundaries.bed
  - H2_INV1936_local.t1.bed
  - H2 broad-region BAM slices for 141L_S13, 142L_S14, 143L_S15

Suggested IGV load order
1. Load the scaffold_5 FASTA for one haplotype.
2. Load the broad-region BRAKER GTF.
3. Load the broad-region RepeatMasker GFF.
4. Load breakpoint BED.
5. Load gene BED.
6. Load inversion interval and boundaries BED.
7. Load the local t1 BED.
8. Load the three BAMs for that haplotype.

Key inspection windows
- H1 broad: {SCAFFOLD}:{H1_BROAD_START:,}-{H1_BROAD_END:,}
- H1 breakpoint zoom: {SCAFFOLD}:{H1_ZOOM_START:,}-{H1_ZOOM_END:,}
- H2 broad: {SCAFFOLD}:{H2_BROAD_START:,}-{H2_BROAD_END:,}
- H2 breakpoint zoom: {SCAFFOLD}:{H2_ZOOM_START:,}-{H2_ZOOM_END:,}

What to inspect
1. Whether the full g13041 / g13243 gene bodies retain read support across both haplotypes.
2. Whether the terminal exon block across the breakpoint retains coverage in both H1 and H2.
3. Whether splice junctions connecting the terminal exons are present in both haplotypes.
4. Whether one haplotype shows a clear loss of terminal-exon support relative to the other.
5. Whether local flanking coverage or junction behavior changes near the boundary.
"""
    out_path.write_text(text)


def main() -> None:
    ensure_dirs()

    extract_scaffold(H1_REF, H1_DIR / "H1_20260127.FINAL.scaffold_5.fa")
    extract_scaffold(H2_REF, H2_DIR / "H2_20260127.FINAL.scaffold_5.fa")

    subset_gtf(H1_GTF, H1_DIR / "braker.INV1936_H1_broad_region.gtf", SCAFFOLD, H1_BROAD_START, H1_BROAD_END)
    subset_gff(H1_GFF, H1_DIR / "H1_20260127.FINAL.top20.fa.out.INV1936_H1_broad_region.gff", SCAFFOLD, H1_BROAD_START, H1_BROAD_END)
    shutil.copy2(H1_BP_BED, H1_DIR / "INV1936_H1_breakpoint.bed")
    write_interval_beds(H1_DIR, "INV1936_H1", SCAFFOLD, H1_INV_START, H1_INV_END, H1_GENE_START, H1_GENE_END, "g13041.t1_H1", "+")
    subset_bed(H1_LOCAL_BED, H1_DIR / "H1_INV1936_local.t1.bed", SCAFFOLD, H1_BROAD_START, H1_BROAD_END)

    subset_gtf(H2_GTF, H2_DIR / "braker.INV1936_H2_broad_region.gtf", SCAFFOLD, H2_BROAD_START, H2_BROAD_END)
    subset_gff(H2_GFF, H2_DIR / "H2_20260127.FINAL.top20.ordered.renamed.fa.out.INV1936_H2_broad_region.gff", SCAFFOLD, H2_BROAD_START, H2_BROAD_END)
    shutil.copy2(H2_BP_BED, H2_DIR / "INV1936_H2_breakpoint.bed")
    write_interval_beds(H2_DIR, "INV1936_H2", SCAFFOLD, H2_INV_START, H2_INV_END, H2_GENE_START, H2_GENE_END, "g13243.t1_H2", "-")
    subset_bed(H2_LOCAL_BED, H2_DIR / "H2_INV1936_local.t1.bed", SCAFFOLD, H2_BROAD_START, H2_BROAD_END)

    for bam in H1_BAMS:
        out_bam = H1_DIR / f"{bam.stem}.INV1936_broad_region.bam"
        slice_bam(bam, out_bam, f"{SCAFFOLD}:{H1_BROAD_START}-{H1_BROAD_END}")

    for bam in H2_BAMS:
        out_bam = H2_DIR / f"{bam.stem}.INV1936_broad_region.bam"
        slice_bam(bam, out_bam, f"{SCAFFOLD}:{H2_BROAD_START}-{H2_BROAD_END}")

    write_manifest(OUTDIR / "INV1936_H1_H2_bam_manifest.tsv")
    write_readme(OUTDIR / "README.txt")


if __name__ == "__main__":
    main()
