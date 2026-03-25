#!/usr/bin/env python3
from __future__ import annotations

import csv
import subprocess
from itertools import groupby
from pathlib import Path


OUTDIR = Path("/g/data/xf3/zz3507/Output/20260127Genome/inv1936_pygenometracks")
SAMTOOLS = Path("/g/data/xf3/miniconda/envs/common-tools/bin/samtools")

H1_REGION = ("scaffold_5", 7_262_500, 7_279_200)
H2_REGION = ("scaffold_5", 7_448_400, 7_464_200)

H1_GTF = Path("/g/data/xf3/zz3507/Output/20260127Genome/H1/breaker/braker.gtf")
H2_GTF = Path("/g/data/xf3/zz3507/Output/20260127Genome/H2/breaker_H2_2026_0211/braker.gtf")
H1_REPEAT = Path("/g/data/xf3/zz3507/Output/20260127Genome/H1/RepeatMasker/H1_20260127.FINAL.top20.fa.out.gff")
H2_REPEAT = Path("/g/data/xf3/zz3507/Output/20260127Genome/H2/RepeatMasker/H2_20260127.FINAL.top20.ordered.renamed.fa.out.gff")

H1_BREAKPOINT_BED = Path("/g/data/xf3/zz3507/Output/20260127Genome/inv1936_truncation_analysis/INV1936_H1_breakpoint.bed")
H2_BREAKPOINT_BED = Path("/g/data/xf3/zz3507/Output/20260127Genome/inv1936_truncation_analysis/INV1936_H2_breakpoint.bed")
H1_GENE_BED = Path("/g/data/xf3/zz3507/Output/20260127Genome/inv1936_truncation_analysis/INV1936_H1_gene.bed")
H2_GENE_BED = Path("/g/data/xf3/zz3507/Output/20260127Genome/inv1936_truncation_analysis/INV1936_H2_gene.bed")

RNA_BAMS = [
    Path("/g/data/xf3/zz3507/Output/20260127Genome/H1/inversion_rna_minimap2/bam/141L_S13.sorted.bam"),
    Path("/g/data/xf3/zz3507/Output/20260127Genome/H1/inversion_rna_minimap2/bam/142L_S14.sorted.bam"),
    Path("/g/data/xf3/zz3507/Output/20260127Genome/H1/inversion_rna_minimap2/bam/143L_S15.sorted.bam"),
]


def overlaps(chrom: str, start: int, end: int, region: tuple[str, int, int]) -> bool:
    rchrom, rstart, rend = region
    return chrom == rchrom and not (end < rstart or start > rend)


def subset_text_track(src: Path, dst: Path, region: tuple[str, int, int], feature_types: set[str] | None = None) -> None:
    with src.open() as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            ftype = parts[2] if len(parts) > 2 else ""
            start = int(parts[3])
            end = int(parts[4])
            if feature_types is not None and ftype not in feature_types:
                continue
            if overlaps(chrom, start, end, region):
                fout.write(line)


def gff_to_bed(src: Path, dst: Path, region: tuple[str, int, int]) -> None:
    with src.open() as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _src, feature, start, end, _score, strand, _phase, attrs = parts
            start_i = int(start)
            end_i = int(end)
            if not overlaps(chrom, start_i, end_i, region):
                continue
            name = attrs.split(";")[0].replace("ID=", "")
            fout.write(f"{chrom}\t{start_i - 1}\t{end_i}\t{name or feature}\t0\t{strand}\n")


def mapped_primary_reads(bam: Path) -> int:
    cmd = [str(SAMTOOLS), "view", "-c", "-F", "2308", str(bam)]
    return int(subprocess.check_output(cmd, text=True).strip())


def bedgraph_from_depth(region: tuple[str, int, int], bams: list[Path], dst: Path, mapped_counts: list[int]) -> None:
    chrom, start, end = region
    region_str = f"{chrom}:{start}-{end}"
    cmd = [str(SAMTOOLS), "depth", "-aa", "-r", region_str] + [str(b) for b in bams]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    if proc.stdout is None:
        raise RuntimeError("Failed to capture samtools depth output")

    rows: list[tuple[int, float]] = []
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        pos = int(parts[1])
        depths = [int(x) for x in parts[2:]]
        scaled = [(d * 1_000_000.0 / n) if n > 0 else 0.0 for d, n in zip(depths, mapped_counts)]
        mean_cov = sum(scaled) / len(scaled)
        rows.append((pos, mean_cov))
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"samtools depth failed with exit code {ret}")

    with dst.open("w") as fout:
        current_start = None
        current_end = None
        current_val = None
        for pos, cov in rows:
            if cov == 0:
                if current_start is not None:
                    fout.write(f"{chrom}\t{current_start - 1}\t{current_end}\t{current_val:.6f}\n")
                    current_start = current_end = current_val = None
                continue
            cov = round(cov, 6)
            if current_start is None:
                current_start = current_end = pos
                current_val = cov
                continue
            if pos == current_end + 1 and cov == current_val:
                current_end = pos
            else:
                fout.write(f"{chrom}\t{current_start - 1}\t{current_end}\t{current_val:.6f}\n")
                current_start = current_end = pos
                current_val = cov
        if current_start is not None:
            fout.write(f"{chrom}\t{current_start - 1}\t{current_end}\t{current_val:.6f}\n")


def write_ini(path: Path, region_key: str) -> None:
    if region_key == "H1":
        files = {
            "coverage": OUTDIR / "INV1936_H1_sunflowerL_meanRPM.bedgraph",
            "gtf": OUTDIR / "INV1936_H1_region.braker.gtf",
            "repeat": OUTDIR / "INV1936_H1_region.repeats.bed",
            "gene": H1_GENE_BED,
            "breakpoint": H1_BREAKPOINT_BED,
        }
        title = "INV1936 H1"
    else:
        files = {
            "gtf": OUTDIR / "INV1936_H2_region.braker.gtf",
            "repeat": OUTDIR / "INV1936_H2_region.repeats.bed",
            "gene": H2_GENE_BED,
            "breakpoint": H2_BREAKPOINT_BED,
        }
        title = "INV1936 H2"

    lines = []
    if region_key == "H1":
        lines.extend(
            [
                "[coverage]",
                f"file = {files['coverage']}",
                "title = Mean RPM-normalized RNA coverage",
                "height = 3.0",
                "file_type = bedgraph",
                "color = #2c7fb8",
                "alpha = 0.8",
                "min_value = 0",
                "nans_to_zeros = true",
                "show_data_range = true",
                "number_of_bins = 500",
                "type = fill",
                "",
            ]
        )

    lines.extend(
        [
            "[repeats]",
            f"file = {files['repeat']}",
            "title = Repeats",
            "height = 1.2",
            "file_type = bed",
            "display = collapsed",
            "labels = false",
            "color = #ef8a62",
            "border_color = none",
            "",
            "[gene_models]",
            f"file = {files['gtf']}",
            f"title = {title} gene model",
            "height = 2.4",
            "file_type = gtf",
            "prefered_name = gene_id",
            "merge_transcripts = false",
            "style = UCSC",
            "labels = true",
            "fontsize = 10",
            "color = #2166ac",
            "border_color = black",
            "",
            "[highlight_gene]",
            f"file = {files['gene']}",
            "title = ",
            "height = 0.7",
            "file_type = bed",
            "display = collapsed",
            "labels = false",
            "color = #fdb863",
            "border_color = #b35806",
            "",
            "[breakpoint]",
            f"file = {files['breakpoint']}",
            "type = vlines",
            "line_style = dashed",
            "color = #d7301f",
            "line_width = 1.8",
            "",
            "[x-axis]",
            "",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def write_manifest(path: Path, mapped_counts: list[int]) -> None:
    with path.open("w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["bam", "mapped_primary_reads", "normalization"])
        for bam, count in zip(RNA_BAMS, mapped_counts):
            writer.writerow([str(bam), count, "RPM-scaled depth"])


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    subset_text_track(H1_GTF, OUTDIR / "INV1936_H1_region.braker.gtf", H1_REGION, {"gene", "transcript", "exon", "CDS", "intron", "start_codon", "stop_codon"})
    subset_text_track(H2_GTF, OUTDIR / "INV1936_H2_region.braker.gtf", H2_REGION, {"gene", "transcript", "exon", "CDS", "intron", "start_codon", "stop_codon"})
    gff_to_bed(H1_REPEAT, OUTDIR / "INV1936_H1_region.repeats.bed", H1_REGION)
    gff_to_bed(H2_REPEAT, OUTDIR / "INV1936_H2_region.repeats.bed", H2_REGION)

    mapped_counts = [mapped_primary_reads(bam) for bam in RNA_BAMS]
    bedgraph_from_depth(H1_REGION, RNA_BAMS, OUTDIR / "INV1936_H1_sunflowerL_meanRPM.bedgraph", mapped_counts)
    write_manifest(OUTDIR / "INV1936_H1_sunflowerL_meanRPM.manifest.tsv", mapped_counts)

    write_ini(OUTDIR / "INV1936_H1_tracks.ini", "H1")
    write_ini(OUTDIR / "INV1936_H2_tracks.ini", "H2")


if __name__ == "__main__":
    main()
