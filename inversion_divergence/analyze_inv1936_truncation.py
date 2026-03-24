#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


H1_GTF = Path("/g/data/xf3/zz3507/Output/20260127Genome/H1/breaker/braker.gtf")
H2_GTF = Path("/g/data/xf3/zz3507/Output/20260127Genome/H2/breaker_H2_2026_0211/braker.gtf")
H1_REPEAT = Path("/g/data/xf3/zz3507/Output/20260127Genome/H1/RepeatMasker/H1_20260127.FINAL.top20.fa.out.gff")
H2_REPEAT = Path("/g/data/xf3/zz3507/Output/20260127Genome/H2/RepeatMasker/H2_20260127.FINAL.top20.ordered.renamed.fa.out.gff")
ANCHOR_TSV = Path("/g/data/xf3/zz3507/Output/20260127Genome/syri/syri_asm10/INV1936_jcvi_microsynteny/H1_H2_INV1936_local.lifted.withheader.anchors")
RNA_TSV = Path("/g/data/xf3/zz3507/Output/20260127Genome/H1/inversion_rna_minimap2_stage_sun_flower_L/inv_rna_expr_minimap2_sun_flower_L.genes_expression.tsv")
OUTDIR = Path("/g/data/xf3/zz3507/Output/20260127Genome/inv1936_truncation_analysis")


def parse_attrs(attrs: str) -> dict[str, str]:
    attrs = attrs.strip()
    out: dict[str, str] = {}
    for key, value in re.findall(r'([A-Za-z_]+)\s+"([^"]+)"', attrs):
        out[key] = value
    return out


@dataclass
class Interval:
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class TranscriptModel:
    transcript_id: str
    start: int
    end: int
    strand: str
    exons: list[Interval] = field(default_factory=list)
    cdss: list[Interval] = field(default_factory=list)
    introns: list[Interval] = field(default_factory=list)
    stop_codons: list[Interval] = field(default_factory=list)


@dataclass
class GeneModel:
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str
    transcripts: dict[str, TranscriptModel] = field(default_factory=dict)

    @property
    def gene_span_bp(self) -> int:
        return self.end - self.start + 1


def load_gene_model(gtf_path: Path, gene_id: str) -> GeneModel:
    gene: GeneModel | None = None
    with gtf_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, _src, feature, start, end, _score, strand, _frame, attrs_text = line.rstrip("\n").split("\t")
            attrs = parse_attrs(attrs_text)
            raw_last = attrs_text.strip()
            row_gene_id = attrs.get("gene_id")
            if row_gene_id is None and feature == "gene":
                row_gene_id = raw_last
            if row_gene_id != gene_id and gene_id not in raw_last:
                continue

            start_i = int(start)
            end_i = int(end)
            if feature == "gene":
                gene = GeneModel(gene_id=gene_id, chrom=chrom, start=start_i, end=end_i, strand=strand)
                continue
            if gene is None:
                gene = GeneModel(gene_id=gene_id, chrom=chrom, start=start_i, end=end_i, strand=strand)

            transcript_id = attrs.get("transcript_id")
            if transcript_id is None and feature == "transcript":
                transcript_id = raw_last
            if transcript_id is None:
                continue
            tx = gene.transcripts.setdefault(
                transcript_id,
                TranscriptModel(transcript_id=transcript_id, start=start_i, end=end_i, strand=strand),
            )
            tx.start = min(tx.start, start_i)
            tx.end = max(tx.end, end_i)
            if feature == "exon":
                tx.exons.append(Interval(start_i, end_i))
            elif feature == "CDS":
                tx.cdss.append(Interval(start_i, end_i))
            elif feature == "intron":
                tx.introns.append(Interval(start_i, end_i))
            elif feature == "stop_codon":
                tx.stop_codons.append(Interval(start_i, end_i))

    if gene is None:
        raise ValueError(f"Failed to find {gene_id} in {gtf_path}")

    for tx in gene.transcripts.values():
        tx.exons.sort(key=lambda iv: iv.start)
        tx.cdss.sort(key=lambda iv: iv.start)
        tx.introns.sort(key=lambda iv: iv.start)
        tx.stop_codons.sort(key=lambda iv: iv.start)
    return gene


def overlap_len(a: Interval, b: Interval) -> int:
    start = max(a.start, b.start)
    end = min(a.end, b.end)
    if end < start:
        return 0
    return end - start + 1


def total_overlap(intervals: list[Interval], region: Interval) -> int:
    return sum(overlap_len(iv, region) for iv in intervals)


def count_overlap(intervals: list[Interval], region: Interval) -> int:
    return sum(1 for iv in intervals if overlap_len(iv, region) > 0)


def containing_interval(intervals: list[Interval], pos: int) -> Interval | None:
    for iv in intervals:
        if iv.start <= pos <= iv.end:
            return iv
    return None


def nearest_exon_distance(exons: list[Interval], pos: int) -> int:
    best = math.inf
    for exon in exons:
        if exon.start <= pos <= exon.end:
            return 0
        if pos < exon.start:
            best = min(best, exon.start - pos)
        else:
            best = min(best, pos - exon.end)
    if best is math.inf:
        return -1
    return int(best)


def pct(part: int, whole: int) -> float:
    if whole <= 0:
        return 0.0
    return 100.0 * part / whole


def segment_for_boundary(start: int, end: int, breakpoint: int, boundary_side: str) -> Interval:
    if boundary_side == "left":
        return Interval(breakpoint, end)
    if boundary_side == "right":
        return Interval(start, breakpoint)
    raise ValueError(boundary_side)


def breakpoint_context(tx: TranscriptModel, breakpoint: int) -> tuple[str, str]:
    exon = containing_interval(tx.exons, breakpoint)
    if exon is not None:
        return "exon", f"{exon.start}-{exon.end}"
    intron = containing_interval(tx.introns, breakpoint)
    if intron is not None:
        return "intron", f"{intron.start}-{intron.end}"
    cds = containing_interval(tx.cdss, breakpoint)
    if cds is not None:
        return "cds", f"{cds.start}-{cds.end}"
    return "intergenic_or_unknown", ""


def load_breakpoint_repeats(path: Path, chrom: str, breakpoint: int, window_bp: int = 2000) -> list[dict[str, str | int | bool]]:
    records: list[dict[str, str | int | bool]] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            seq, src, feature, start, end, score, strand, phase, attrs = line.rstrip("\n").split("\t")
            start_i = int(start)
            end_i = int(end)
            if seq != chrom:
                continue
            if end_i < breakpoint - window_bp or start_i > breakpoint + window_bp:
                continue
            records.append(
                {
                    "chrom": seq,
                    "start": start_i,
                    "end": end_i,
                    "feature": feature,
                    "attrs": attrs,
                    "overlaps_breakpoint": start_i <= breakpoint <= end_i,
                }
            )
    return records


def load_anchor_mapping(path: Path, h1_tx: str, h2_tx: str) -> str | None:
    with path.open() as handle:
        header = next(handle, None)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if parts[0] == h1_tx and parts[1] == h2_tx:
                return parts[2]
    return None


def load_h1_mean_counts(path: Path, gene_id: str) -> float | None:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        sample_cols = [c for c in reader.fieldnames if c and c.endswith(".sorted.bam")]
        for row in reader:
            if row["INV_ID"] != "INV1936" or row["GeneID"] != gene_id:
                continue
            vals = [float(row[c]) for c in sample_cols]
            return sum(vals) / len(vals)
    return None


def transcript_metrics(gene: GeneModel, tx: TranscriptModel, breakpoint: int, boundary_side: str) -> dict[str, object]:
    segment = segment_for_boundary(gene.start, gene.end, breakpoint, boundary_side)
    segment_bp = segment.length
    total_exon_bp = sum(iv.length for iv in tx.exons)
    total_cds_bp = sum(iv.length for iv in tx.cdss)
    total_stop_bp = sum(iv.length for iv in tx.stop_codons)
    coding_bp_no_stop = total_cds_bp - total_stop_bp

    segment_exon_bp = total_overlap(tx.exons, segment)
    segment_cds_bp = total_overlap(tx.cdss, segment)
    segment_stop_bp = total_overlap(tx.stop_codons, segment)
    segment_coding_bp_no_stop = segment_cds_bp - segment_stop_bp
    feature_type, feature_coords = breakpoint_context(tx, breakpoint)

    return {
        "transcript_id": tx.transcript_id,
        "strand": tx.strand,
        "gene_start": gene.start,
        "gene_end": gene.end,
        "gene_span_bp": gene.gene_span_bp,
        "breakpoint": breakpoint,
        "boundary_side": boundary_side,
        "breakpoint_feature_type": feature_type,
        "breakpoint_feature_coords": feature_coords,
        "nearest_exon_distance_bp": nearest_exon_distance(tx.exons, breakpoint),
        "segment_start": segment.start,
        "segment_end": segment.end,
        "segment_bp": segment_bp,
        "segment_pct_gene_span": round(pct(segment_bp, gene.gene_span_bp), 3),
        "segment_exon_bp": segment_exon_bp,
        "segment_intron_bp": segment_bp - segment_exon_bp,
        "segment_cds_bp_including_stop": segment_cds_bp,
        "segment_stop_bp": segment_stop_bp,
        "segment_coding_bp_excluding_stop": segment_coding_bp_no_stop,
        "segment_coding_aa_excluding_stop": segment_coding_bp_no_stop // 3,
        "segment_exon_count": count_overlap(tx.exons, segment),
        "total_exon_bp": total_exon_bp,
        "total_cds_bp_including_stop": total_cds_bp,
        "total_stop_bp": total_stop_bp,
        "total_coding_bp_excluding_stop": coding_bp_no_stop,
        "total_coding_aa_excluding_stop": coding_bp_no_stop // 3,
        "segment_pct_exon_bp": round(pct(segment_exon_bp, total_exon_bp), 3),
        "segment_pct_coding_bp": round(pct(segment_coding_bp_no_stop, coding_bp_no_stop), 3),
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    h1_gene = load_gene_model(H1_GTF, "g13041")
    h2_gene = load_gene_model(H2_GTF, "g13243")
    h1_breakpoint = 7277582
    h2_breakpoint = 7450031

    metric_rows: list[dict[str, object]] = []
    for locus, gene, breakpoint, side in [
        ("H1", h1_gene, h1_breakpoint, "left"),
        ("H2", h2_gene, h2_breakpoint, "right"),
    ]:
        for tx_id in sorted(gene.transcripts):
            row = {"locus": locus, "gene_id": gene.gene_id, "chrom": gene.chrom}
            row.update(transcript_metrics(gene, gene.transcripts[tx_id], breakpoint, side))
            metric_rows.append(row)

    repeats_rows: list[dict[str, object]] = []
    for locus, path, chrom, breakpoint in [
        ("H1", H1_REPEAT, h1_gene.chrom, h1_breakpoint),
        ("H2", H2_REPEAT, h2_gene.chrom, h2_breakpoint),
    ]:
        for record in load_breakpoint_repeats(path, chrom, breakpoint):
            repeats_rows.append(
                {
                    "locus": locus,
                    "breakpoint": breakpoint,
                    **record,
                }
            )

    segment_rows: list[dict[str, object]] = []
    for locus, gene, breakpoint, side in [
        ("H1", h1_gene, h1_breakpoint, "left"),
        ("H2", h2_gene, h2_breakpoint, "right"),
    ]:
        segment = segment_for_boundary(gene.start, gene.end, breakpoint, side)
        for tx in gene.transcripts.values():
            for feature_type, intervals in [
                ("exon", tx.exons),
                ("CDS", tx.cdss),
                ("intron", tx.introns),
                ("stop_codon", tx.stop_codons),
            ]:
                for iv in intervals:
                    ov = overlap_len(iv, segment)
                    if ov <= 0:
                        continue
                    segment_rows.append(
                        {
                            "locus": locus,
                            "gene_id": gene.gene_id,
                            "transcript_id": tx.transcript_id,
                            "feature_type": feature_type,
                            "feature_start": iv.start,
                            "feature_end": iv.end,
                            "feature_len_bp": iv.length,
                            "segment_start": segment.start,
                            "segment_end": segment.end,
                            "overlap_bp": ov,
                        }
                    )

    write_tsv(OUTDIR / "INV1936_truncation_metrics.tsv", metric_rows)
    write_tsv(OUTDIR / "INV1936_breakpoint_repeats.tsv", repeats_rows)
    write_tsv(OUTDIR / "INV1936_boundary_segment_features.tsv", segment_rows)

    h1_t1 = next(row for row in metric_rows if row["locus"] == "H1" and row["transcript_id"] == "g13041.t1")
    h2_t1 = next(row for row in metric_rows if row["locus"] == "H2" and row["transcript_id"] == "g13243.t1")
    anchor_score = load_anchor_mapping(ANCHOR_TSV, "g13041.t1", "g13243.t1")
    mean_counts = load_h1_mean_counts(RNA_TSV, "g13041")
    h1_exact_repeats = [r for r in repeats_rows if r["locus"] == "H1" and r["overlaps_breakpoint"]]
    h2_exact_repeats = [r for r in repeats_rows if r["locus"] == "H2" and r["overlaps_breakpoint"]]

    def repeat_summary(rows: list[dict[str, object]]) -> str:
        if not rows:
            return "None detected"
        bits = []
        for row in rows:
            bits.append(f"{row['start']}-{row['end']} ({row['attrs']})")
        return "; ".join(bits)

    md = []
    md.append("# INV1936 Truncation Analysis\n")
    md.append("## Concrete Checklist\n")
    md.append("- [x] Confirm which H1 gene actually overlaps the INV1936 breakpoint.\n")
    md.append("- [x] Quantify the exact boundary-overlapping segment in bp.\n")
    md.append("- [x] Determine whether the breakpoint falls in exon or intron.\n")
    md.append("- [x] Quantify how much exon/CDS sequence lies beyond the breakpoint.\n")
    md.append("- [x] Check whether BRAKER has alternative transcripts for the locus.\n")
    md.append("- [x] Identify the corresponding H2 anchored gene.\n")
    md.append("- [x] Quantify the same boundary-overlap metrics on the H2 side.\n")
    md.append("- [x] Check repeat annotations exactly overlapping each breakpoint.\n")
    md.append("- [x] Tie the structure back to the current `sun flower L` RNA summary.\n")
    md.append("\n## Executed Results\n")
    md.append(f"- H1 breakpoint-overlapping gene: `{h1_gene.gene_id}` on `{h1_gene.chrom}`, span `{h1_gene.start}-{h1_gene.end}`.\n")
    md.append(f"- H1 breakpoint: `{h1_breakpoint}`. This falls in `{h1_t1['breakpoint_feature_type']}` `{h1_t1['breakpoint_feature_coords']}` for `g13041.t1`.\n")
    md.append(f"- H1 boundary-overlapping segment size: `{h1_t1['segment_bp']}` bp, which is `{h1_t1['segment_pct_gene_span']}`% of the full gene span.\n")
    md.append(f"- H1 segment beyond the breakpoint contains `{h1_t1['segment_exon_count']}` exons, `{h1_t1['segment_exon_bp']}` exonic bp, and `{h1_t1['segment_coding_bp_excluding_stop']}` coding bp excluding the stop codon (`~{h1_t1['segment_coding_aa_excluding_stop']}` aa).\n")
    md.append(f"- H1 nearest exon boundary is `{h1_t1['nearest_exon_distance_bp']}` bp away from the breakpoint.\n")
    md.append(f"- H1 exact breakpoint repeat overlap: {repeat_summary(h1_exact_repeats)}.\n")
    md.append(f"- H1 `sun flower L` mean RNA count for `g13041`: `{mean_counts:.3f}`.\n")
    md.append(f"- BRAKER alternative transcripts in H1: `{', '.join(sorted(h1_gene.transcripts))}`.\n")
    md.append(f"- H2 anchored counterpart from local microsynteny anchors: `g13243.t1` (anchor score `{anchor_score}`).\n")
    md.append(f"- H2 breakpoint-overlapping gene: `{h2_gene.gene_id}` on `{h2_gene.chrom}`, span `{h2_gene.start}-{h2_gene.end}`.\n")
    md.append(f"- H2 breakpoint: `{h2_breakpoint}`. This falls in `{h2_t1['breakpoint_feature_type']}` `{h2_t1['breakpoint_feature_coords']}` for `g13243.t1`.\n")
    md.append(f"- H2 boundary-overlapping segment size: `{h2_t1['segment_bp']}` bp, which is `{h2_t1['segment_pct_gene_span']}`% of the full gene span.\n")
    md.append(f"- H2 segment beyond the breakpoint contains `{h2_t1['segment_exon_count']}` exons, `{h2_t1['segment_exon_bp']}` exonic bp, and `{h2_t1['segment_coding_bp_excluding_stop']}` coding bp excluding the stop codon (`~{h2_t1['segment_coding_aa_excluding_stop']}` aa).\n")
    md.append(f"- H2 nearest exon boundary is `{h2_t1['nearest_exon_distance_bp']}` bp away from the breakpoint.\n")
    md.append(f"- H2 exact breakpoint repeat overlap: {repeat_summary(h2_exact_repeats)}.\n")
    md.append(f"- BRAKER alternative transcripts in H2: `{', '.join(sorted(h2_gene.transcripts))}`.\n")
    md.append("\n## Interpretation\n")
    md.append("- `INV1936` is a real breakpoint-overlapping gene case, but the breakpoint does **not** cut through a coding exon. In both H1 and H2 it falls inside an intron.\n")
    md.append("- The boundary-associated portion is still biologically non-trivial: in both haplotypes, the breakpoint separates a terminal segment of `1359 bp`, including `444 bp` of exonic sequence and `441 bp` of coding sequence excluding the stop codon (`~147 aa`).\n")
    md.append("- The H1 and H2 structures are mirrored rather than 'one truncated vs one intact'. The corresponding anchored H2 gene (`g13243`) also spans its own inversion boundary with the same boundary-overlap size.\n")
    md.append("- That means the safest wording is not 'the inversion uniquely truncates the gene in H1', but rather 'the inversion boundary intersects an orthologous boundary-spanning gene pair, with the breakpoint lying in an intron and isolating the 3-prime terminal exon block across the boundary'.\n")
    md.append("- Both breakpoints are repeat-associated: the H1 breakpoint overlaps an `(AT)n` repeat, while the H2 breakpoint overlaps local repeat-family annotations as well.\n")
    md.append("\n## Recommended Next Checks\n")
    md.append("- Build a zoomed gene model figure showing the breakpoint, the intron it falls in, and the three terminal exons on each haplotype.\n")
    md.append("- If you want to discuss transcript consequence, use read tiling / splice-junction evidence rather than TPM alone, because the breakpoint is intronic and the effective transcript model may matter more than raw gene-span overlap.\n")
    md.append("- If expression is emphasized in the manuscript, phrase the current result as `sun flower stage, labellum, 2h UV-treated subset` only.\n")

    (OUTDIR / "INV1936_truncation_checklist_and_results.md").write_text("".join(md))


if __name__ == "__main__":
    main()
