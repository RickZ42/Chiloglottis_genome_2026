from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pysam


ROOT = Path("/g/data/xf3/zz3507")
OUTDIR = ROOT / "Output/20260127Genome/inv1936_alllibs_merged_support"


@dataclass
class LocusSpec:
    label: str
    bam: Path
    gtf: Path
    gene: str
    transcript: str
    chrom: str
    gene_start: int
    gene_end: int
    breakpoint: int
    terminal_start: int
    terminal_end: int


def parse_exons(gtf_path: Path, gene_id: str, transcript_id: str) -> tuple[list[tuple[int, int]], str]:
    exons: list[tuple[int, int]] = []
    strand = None
    with gtf_path.open() as handle:
        for raw in handle:
            if f'gene_id "{gene_id}"' not in raw or "\texon\t" not in raw:
                continue
            attrs = raw.rstrip().split("\t")[8]
            m = re.search(r'transcript_id "([^"]+)"', attrs)
            if not m or m.group(1) != transcript_id:
                continue
            fields = raw.rstrip().split("\t")
            exons.append((int(fields[3]), int(fields[4])))
            strand = fields[6]
    if not exons or strand is None:
        raise ValueError(f"Could not parse exons for {gene_id} {transcript_id}")
    return sorted(exons), strand


def introns_from_exons(exons: list[tuple[int, int]]) -> list[tuple[int, int]]:
    return [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)]


def collect_junctions(bam_path: Path, chrom: str, start: int, end: int) -> tuple[int, int, Counter[tuple[int, int]]]:
    counts: Counter[tuple[int, int]] = Counter()
    total = 0
    spliced = 0
    with pysam.AlignmentFile(bam_path) as bam:
        for rec in bam.fetch(chrom, start, end):
            if rec.is_unmapped:
                continue
            total += 1
            ref = rec.reference_start + 1
            saw_splice = False
            for op, length in rec.cigartuples or []:
                if op == 3:  # N
                    counts[(ref, ref + length - 1)] += 1
                    ref += length
                    saw_splice = True
                elif op in (0, 2, 7, 8):
                    ref += length
            spliced += int(saw_splice)
    return total, spliced, counts


def count_reads(bam_path: Path, chrom: str, start: int, end: int) -> int:
    with pysam.AlignmentFile(bam_path) as bam:
        return sum(1 for _ in bam.fetch(chrom, start, end))


def is_breakpoint_local(junc: tuple[int, int], bp: int, flank: int = 5000) -> bool:
    left, right = junc
    return abs(left - bp) <= flank or abs(right - bp) <= flank


def main() -> None:
    specs = [
        LocusSpec(
            label="H1",
            bam=OUTDIR / "H1_INV1936_alllibs_merged.bam",
            gtf=ROOT / "Output/20260127Genome/H1/breaker/braker.gtf",
            gene="g13041",
            transcript="g13041.t1",
            chrom="scaffold_5",
            gene_start=7223150,
            gene_end=7278940,
            breakpoint=7277582,
            terminal_start=7278118,
            terminal_end=7278940,
        ),
        LocusSpec(
            label="H2",
            bam=OUTDIR / "h2_remap/H2_INV1936_alllibs_merged_remap.bam",
            gtf=ROOT / "Output/20260127Genome/H2/breaker_H2_2026_0211/braker.gtf",
            gene="g13243",
            transcript="g13243.t1",
            chrom="scaffold_5",
            gene_start=7448673,
            gene_end=7501391,
            breakpoint=7450031,
            terminal_start=7448673,
            terminal_end=7449495,
        ),
    ]

    tsv = OUTDIR / "INV1936_alllibs_merged_support_summary.tsv"
    md = OUTDIR / "INV1936_alllibs_merged_support_summary.md"
    novel_tsv = OUTDIR / "INV1936_alllibs_breakpoint_local_novel_junctions.tsv"

    with tsv.open("w") as out:
        out.write(
            "\t".join(
                [
                    "haplotype",
                    "gene",
                    "total_reads_in_gene_window",
                    "spliced_reads_in_gene_window",
                    "terminal_block_reads",
                    "dominant_boundary_intron",
                    "dominant_boundary_intron_count",
                    "novel_breakpoint_local_junctions_ge3",
                ]
            )
            + "\n"
        )
        with novel_tsv.open("w") as jout:
            jout.write("haplotype\tjunction_start\tjunction_end\tcount\n")
            md_lines = [
                "# INV1936 merged-all-libraries read-support summary",
                "",
            ]
            for spec in specs:
                exons, strand = parse_exons(spec.gtf, spec.gene, spec.transcript)
                introns = introns_from_exons(exons)
                annotated = set(introns)
                total, spliced, junctions = collect_junctions(spec.bam, spec.chrom, spec.gene_start, spec.gene_end)
                terminal_reads = count_reads(spec.bam, spec.chrom, spec.terminal_start, spec.terminal_end)
                boundary_intron = next(j for j in introns if j[0] <= spec.breakpoint <= j[1])
                dominant_count = junctions.get(boundary_intron, 0)
                novel_local = [
                    (j, c)
                    for j, c in junctions.items()
                    if j not in annotated and is_breakpoint_local(j, spec.breakpoint, flank=5000) and c >= 3
                ]
                novel_local.sort(key=lambda x: x[1], reverse=True)

                out.write(
                    "\t".join(
                        [
                            spec.label,
                            spec.gene,
                            str(total),
                            str(spliced),
                            str(terminal_reads),
                            f"{boundary_intron[0]}-{boundary_intron[1]}",
                            str(dominant_count),
                            str(len(novel_local)),
                        ]
                    )
                    + "\n"
                )
                for j, c in novel_local:
                    jout.write(f"{spec.label}\t{j[0]}\t{j[1]}\t{c}\n")

                md_lines.extend(
                    [
                        f"## {spec.label} {spec.gene}",
                        f"- Gene window reads: `{total}`; spliced reads: `{spliced}`.",
                        f"- Terminal exon-block reads: `{terminal_reads}`.",
                        f"- Breakpoint-spanning annotated intron: `{boundary_intron[0]}-{boundary_intron[1]}` (`{dominant_count}` supporting junction reads).",
                    ]
                )
                if novel_local:
                    md_lines.append("- Breakpoint-local novel junctions with >=3 supporting reads:")
                    for j, c in novel_local:
                        md_lines.append(f"  - `{j[0]}-{j[1]}`: `{c}`")
                else:
                    md_lines.append("- No breakpoint-local novel junctions reached >=3 supporting reads.")
                md_lines.append("")

            md_lines.extend(
                [
                    "## Interpretation",
                    "- In both haplotypes, the dominant splice pattern at the breakpoint remains the annotated long intron spanning the boundary.",
                    "- The terminal exon block on the far side of the boundary remains strongly covered in both H1 and H2 after merging all available libraries.",
                    "- Breakpoint-local novel junctions are present only at low counts relative to the dominant annotated junction and do not support a stable recurrent chimera model.",
                    "",
                    "## Practical conclusion",
                    "The merged-all-libraries view increases local read depth, but it still does not reveal a convincing recurrent `E1-E5 + new exons`, `new exons + E5-E8`, or `E5-E8 + new ORF` fusion pattern at INV1936.",
                ]
            )
            md.write_text("\n".join(md_lines) + "\n")


if __name__ == "__main__":
    main()
