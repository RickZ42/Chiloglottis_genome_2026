#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import subprocess
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import pandas as pd


ROOT = Path("/g/data/xf3/zz3507")
SAMTOOLS = Path("/g/data/xf3/miniconda/envs/common-tools/bin/samtools")
MINIMAP2 = Path("/g/data/xf3/miniconda/bin/minimap2")

INV_TABLE = ROOT / "Output/20260127Genome/syri/syri_asm10/H1_vs_H2syri.highconfINV.tsv"
H1_REF = ROOT / "Output/20260127Genome/H1/H1_20260127.FINAL.top20.fa"
H2_REF = ROOT / "Output/20260127Genome/H2/H2_20260127.FINAL.top20.ordered.renamed.fa"
H1_GTF = ROOT / "Output/20260127Genome/H1/breaker/braker.gtf"
H2_GTF = ROOT / "Output/20260127Genome/H2/breaker_H2_2026_0211/braker.gtf"
H1_GFF = ROOT / "Output/20260127Genome/H1/RepeatMasker/H1_20260127.FINAL.top20.fa.out.gff"
H2_GFF = ROOT / "Output/20260127Genome/H2/RepeatMasker/H2_20260127.FINAL.top20.ordered.renamed.fa.out.gff"
H1_T1_BED = ROOT / "compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H1t1.bed"
H2_T1_BED = ROOT / "compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H2t1.bed"

DEFAULT_H1_BAMS = [
    ROOT / "Output/20260127Genome/H1/inversion_rna_minimap2/bam/141L_S13.sorted.bam",
    ROOT / "Output/20260127Genome/H1/inversion_rna_minimap2/bam/142L_S14.sorted.bam",
    ROOT / "Output/20260127Genome/H1/inversion_rna_minimap2/bam/143L_S15.sorted.bam",
]


@dataclass(frozen=True)
class GeneBedRecord:
    chrom: str
    start1: int
    end1: int
    tx_id: str
    strand: str

    @property
    def gene_id(self) -> str:
        return self.tx_id.rsplit(".", 1)[0]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inv-id", required=True)
    p.add_argument(
        "--outdir",
        type=Path,
        default=ROOT / "Output/20260127Genome/other_syri_inversion",
    )
    p.add_argument("--flank-bp", type=int, default=10000)
    p.add_argument("--padding-bp", type=int, default=5000)
    p.add_argument("--zoom-bp", type=int, default=2500)
    return p.parse_args()


def load_scaffold_lengths(fai_path: Path) -> Dict[str, int]:
    d: Dict[str, int] = {}
    with fai_path.open() as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            d[fields[0]] = int(fields[1])
    return d


def load_bed_records(path: Path) -> Dict[str, List[GeneBedRecord]]:
    by_chr: Dict[str, List[GeneBedRecord]] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start0, end1, tx_id, *_rest = line.rstrip("\n").split("\t")
            strand = _rest[1] if len(_rest) >= 2 else "+"
            rec = GeneBedRecord(chrom=chrom, start1=int(start0) + 1, end1=int(end1), tx_id=tx_id, strand=strand)
            by_chr.setdefault(chrom, []).append(rec)
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda x: (x.start1, x.end1, x.tx_id))
    return by_chr


def gene_distance(rec: GeneBedRecord, pos: int) -> int:
    if rec.start1 <= pos <= rec.end1:
        return 0
    if pos < rec.start1:
        return rec.start1 - pos
    return pos - rec.end1


def nearby_genes(records: Sequence[GeneBedRecord], left_bp: int, right_bp: int, flank_bp: int) -> List[GeneBedRecord]:
    return [r for r in records if min(gene_distance(r, left_bp), gene_distance(r, right_bp)) <= flank_bp]


def expand_region(chrom_len: int, starts: Sequence[int], ends: Sequence[int], padding_bp: int) -> Tuple[int, int]:
    return max(1, min(starts) - padding_bp), min(chrom_len, max(ends) + padding_bp)


def write_single_record_fai(fasta_path: Path) -> None:
    with fasta_path.open("rb") as fh:
        header = fh.readline()
        if not header.startswith(b">"):
            raise ValueError(f"{fasta_path} is not FASTA")
        name = header[1:].strip().split()[0].decode()
        seq_offset = fh.tell()
        first_seq = fh.readline()
        if not first_seq:
            raise ValueError(f"{fasta_path} has no sequence")
        linebases = len(first_seq.rstrip(b"\r\n"))
        linewidth = len(first_seq)
        total_len = linebases
        for line in fh:
            total_len += len(line.rstrip(b"\r\n"))
    with fasta_path.with_suffix(fasta_path.suffix + ".fai").open("w") as out:
        out.write(f"{name}\t{total_len}\t{seq_offset}\t{linebases}\t{linewidth}\n")


def extract_scaffold(ref_fa: Path, out_fa: Path, scaffold: str) -> None:
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
    write_single_record_fai(out_fa)


def _iter_lines(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            yield line


def overlap(start0: int, end0: int, qstart0: int, qend0: int) -> bool:
    return start0 < qend0 and end0 > qstart0


def subset_gtf(in_path: Path, out_path: Path, chrom: str, start1: int, end1: int) -> None:
    with out_path.open("w") as out:
        for line in _iter_lines(in_path):
            if not line.strip() or line.startswith("#"):
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != chrom:
                continue
            if overlap(int(parts[3]) - 1, int(parts[4]), start1 - 1, end1):
                out.write(line)


def subset_gff(in_path: Path, out_path: Path, chrom: str, start1: int, end1: int) -> None:
    with out_path.open("w") as out:
        for line in _iter_lines(in_path):
            if not line.strip() or line.startswith("#"):
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != chrom:
                continue
            if overlap(int(parts[3]) - 1, int(parts[4]), start1 - 1, end1):
                out.write(line)


def subset_bed_records(records: Sequence[GeneBedRecord], out_path: Path, chrom: str, start1: int, end1: int) -> List[GeneBedRecord]:
    kept: List[GeneBedRecord] = []
    with out_path.open("w") as out:
        for rec in records:
            if rec.chrom != chrom:
                continue
            if overlap(rec.start1 - 1, rec.end1, start1 - 1, end1):
                out.write(f"{rec.chrom}\t{rec.start1 - 1}\t{rec.end1}\t{rec.tx_id}\t0\t{rec.strand}\n")
                kept.append(rec)
    return kept


def write_interval_beds(outdir: Path, prefix: str, chrom: str, start1: int, end1: int) -> None:
    with (outdir / f"{prefix}_interval.bed").open("w") as fh:
        fh.write(f"{chrom}\t{start1 - 1}\t{end1}\t{prefix}_interval\t0\t+\n")
    with (outdir / f"{prefix}_boundaries.bed").open("w") as fh:
        fh.write(f"{chrom}\t{start1 - 1}\t{start1}\t{prefix}_left_boundary\t0\t+\n")
        fh.write(f"{chrom}\t{end1 - 1}\t{end1}\t{prefix}_right_boundary\t0\t+\n")


def write_gene_bed(outdir: Path, prefix: str, records: Sequence[GeneBedRecord]) -> None:
    with (outdir / f"{prefix}_genes.bed").open("w") as fh:
        for rec in records:
            fh.write(f"{rec.chrom}\t{rec.start1 - 1}\t{rec.end1}\t{rec.gene_id}\t0\t{rec.strand}\n")


def slice_bam(in_bam: Path, out_bam: Path, region: str) -> None:
    subprocess.run([str(SAMTOOLS), "view", "-b", "-o", str(out_bam), str(in_bam), region], check=True)
    subprocess.run([str(SAMTOOLS), "index", str(out_bam)], check=True)


def bam_to_fastq(in_bam: Path, out_fastq_gz: Path) -> None:
    with gzip.open(out_fastq_gz, "wb") as fh:
        proc = subprocess.run(
            [str(SAMTOOLS), "fastq", "-0", "/dev/stdout", "-n", str(in_bam)],
            check=True,
            stdout=subprocess.PIPE,
        )
        fh.write(proc.stdout)


def remap_to_h2(h2_ref: Path, fastq_gz: Path, out_bam: Path) -> None:
    cmd = (
        f"{MINIMAP2} -ax splice --secondary=no {h2_ref} {fastq_gz} "
        f"| {SAMTOOLS} view -b - "
        f"| {SAMTOOLS} sort -o {out_bam}"
    )
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    subprocess.run([str(SAMTOOLS), "index", str(out_bam)], check=True)


def count_reads(bam: Path) -> int:
    res = subprocess.run([str(SAMTOOLS), "view", "-c", str(bam)], check=True, capture_output=True, text=True)
    return int(res.stdout.strip())


def write_readme(
    out_path: Path,
    inv: pd.Series,
    h1_genes: Sequence[GeneBedRecord],
    h2_genes: Sequence[GeneBedRecord],
    h1_region: Tuple[int, int],
    h2_region: Tuple[int, int],
    zoom_bp: int,
    read_counts: Dict[str, Dict[str, int]],
) -> None:
    text = f"""IGV-ready H1/H2 comparison bundle for {inv['ID']}

Samples
- 141L_S13
- 142L_S14
- 143L_S15

Inversion coordinates
- H1: {inv['RefChr']}:{int(inv['RefStart']):,}-{int(inv['RefEnd']):,}
- H2: {inv['QryChr']}:{int(inv['QryStart']):,}-{int(inv['QryEnd']):,}

Suggested IGV windows
- H1 broad: {inv['RefChr']}:{h1_region[0]:,}-{h1_region[1]:,}
- H1 left-breakpoint zoom: {inv['RefChr']}:{max(1, int(inv['RefStart']) - zoom_bp):,}-{int(inv['RefStart']) + zoom_bp:,}
- H1 right-breakpoint zoom: {inv['RefChr']}:{max(1, int(inv['RefEnd']) - zoom_bp):,}-{int(inv['RefEnd']) + zoom_bp:,}
- H2 broad: {inv['QryChr']}:{h2_region[0]:,}-{h2_region[1]:,}
- H2 left-breakpoint zoom: {inv['QryChr']}:{max(1, int(inv['QryStart']) - zoom_bp):,}-{int(inv['QryStart']) + zoom_bp:,}
- H2 right-breakpoint zoom: {inv['QryChr']}:{max(1, int(inv['QryEnd']) - zoom_bp):,}-{int(inv['QryEnd']) + zoom_bp:,}

Local H1 genes
{chr(10).join('- ' + rec.gene_id + f' ({rec.start1:,}-{rec.end1:,}, {rec.strand})' for rec in h1_genes) if h1_genes else '- none'}

Local H2 genes
{chr(10).join('- ' + rec.gene_id + f' ({rec.start1:,}-{rec.end1:,}, {rec.strand})' for rec in h2_genes) if h2_genes else '- none'}

Read counts
{chr(10).join(f"- {sample}: H1={vals['H1']}, H2={vals['H2']}" for sample, vals in read_counts.items())}

Suggested load order
1. Load the scaffold FASTA for one haplotype.
2. Load the local BRAKER GTF.
3. Load the local RepeatMasker GFF.
4. Load breakpoint / boundaries / interval BED.
5. Load local t1 BED and genes BED.
6. Load the three BAMs.
"""
    out_path.write_text(text)


def main() -> None:
    args = parse_args()
    inv_df = pd.read_csv(INV_TABLE, sep="\t", dtype={"ID": str})
    row = inv_df[inv_df["ID"].astype(str) == args.inv_id]
    if row.empty:
        raise SystemExit(f"Could not find {args.inv_id}")
    inv = row.iloc[0]

    bundle_root = args.outdir / f"{args.inv_id}_H1_H2_IGV_comparison_bundle"
    h1_dir = bundle_root / "H1"
    h2_dir = bundle_root / "H2"
    work_dir = bundle_root / "_work"
    for d in (h1_dir, h2_dir, work_dir):
        d.mkdir(parents=True, exist_ok=True)

    h1_lengths = load_scaffold_lengths(H1_REF.with_suffix(H1_REF.suffix + ".fai"))
    h2_lengths = load_scaffold_lengths(H2_REF.with_suffix(H2_REF.suffix + ".fai"))
    h1_bed = load_bed_records(H1_T1_BED)
    h2_bed = load_bed_records(H2_T1_BED)

    h1_chr = str(inv["RefChr"])
    h2_chr = str(inv["QryChr"])
    h1_start = int(inv["RefStart"])
    h1_end = int(inv["RefEnd"])
    h2_start = int(inv["QryStart"])
    h2_end = int(inv["QryEnd"])

    h1_local = nearby_genes(h1_bed.get(h1_chr, []), h1_start, h1_end, args.flank_bp)
    h2_local = nearby_genes(h2_bed.get(h2_chr, []), h2_start, h2_end, args.flank_bp)

    h1_region = expand_region(
        h1_lengths[h1_chr],
        [h1_start, h1_end] + [g.start1 for g in h1_local],
        [h1_start, h1_end] + [g.end1 for g in h1_local],
        args.padding_bp,
    )
    h2_region = expand_region(
        h2_lengths[h2_chr],
        [h2_start, h2_end] + [g.start1 for g in h2_local],
        [h2_start, h2_end] + [g.end1 for g in h2_local],
        args.padding_bp,
    )

    extract_scaffold(H1_REF, h1_dir / f"H1_20260127.FINAL.{h1_chr}.fa", h1_chr)
    extract_scaffold(H2_REF, h2_dir / f"H2_20260127.FINAL.{h2_chr}.fa", h2_chr)
    subset_gtf(H1_GTF, h1_dir / f"braker.{args.inv_id}_H1_broad_region.gtf", h1_chr, *h1_region)
    subset_gtf(H2_GTF, h2_dir / f"braker.{args.inv_id}_H2_broad_region.gtf", h2_chr, *h2_region)
    subset_gff(H1_GFF, h1_dir / f"H1_20260127.FINAL.top20.fa.out.{args.inv_id}_H1_broad_region.gff", h1_chr, *h1_region)
    subset_gff(H2_GFF, h2_dir / f"H2_20260127.FINAL.top20.ordered.renamed.fa.out.{args.inv_id}_H2_broad_region.gff", h2_chr, *h2_region)
    h1_kept = subset_bed_records(h1_bed.get(h1_chr, []), h1_dir / f"H1_{args.inv_id}_local.t1.bed", h1_chr, *h1_region)
    h2_kept = subset_bed_records(h2_bed.get(h2_chr, []), h2_dir / f"H2_{args.inv_id}_local.t1.bed", h2_chr, *h2_region)
    write_interval_beds(h1_dir, f"{args.inv_id}_H1", h1_chr, h1_start, h1_end)
    write_interval_beds(h2_dir, f"{args.inv_id}_H2", h2_chr, h2_start, h2_end)
    write_gene_bed(h1_dir, f"{args.inv_id}_H1", h1_kept)
    write_gene_bed(h2_dir, f"{args.inv_id}_H2", h2_kept)

    read_counts: Dict[str, Dict[str, int]] = {}
    h1_region_str = f"{h1_chr}:{h1_region[0]}-{h1_region[1]}"

    for bam in DEFAULT_H1_BAMS:
        sample = bam.name.replace(".sorted.bam", "")
        h1_slice = h1_dir / f"{sample}.{args.inv_id}.broad_region.bam"
        h2_remap = h2_dir / f"{sample}.{args.inv_id}.H2_targeted.bam"
        fastq_gz = work_dir / f"{sample}.{args.inv_id}.fastq.gz"
        slice_bam(bam, h1_slice, h1_region_str)
        bam_to_fastq(h1_slice, fastq_gz)
        remap_to_h2(h2_dir / f"H2_20260127.FINAL.{h2_chr}.fa", fastq_gz, h2_remap)
        read_counts[sample] = {"H1": count_reads(h1_slice), "H2": count_reads(h2_remap)}

    write_readme(bundle_root / "README.txt", inv, h1_kept, h2_kept, h1_region, h2_region, args.zoom_bp, read_counts)

    tar_path = bundle_root.with_suffix(".tar.gz")
    with tarfile.open(tar_path, "w:gz") as tar:
        for path in [h1_dir, h2_dir, bundle_root / "README.txt"]:
            tar.add(path, arcname=f"{bundle_root.name}/{path.relative_to(bundle_root)}")

    print(f"Wrote bundle directory: {bundle_root}")
    print(f"Wrote tarball: {tar_path}")
    print(f"H1 broad region: {h1_chr}:{h1_region[0]}-{h1_region[1]}")
    print(f"H2 broad region: {h2_chr}:{h2_region[0]}-{h2_region[1]}")
    for sample, vals in read_counts.items():
        print(f"{sample}: H1={vals['H1']} H2={vals['H2']}")


if __name__ == "__main__":
    main()
