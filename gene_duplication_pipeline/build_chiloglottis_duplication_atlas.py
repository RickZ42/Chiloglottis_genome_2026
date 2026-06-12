#!/usr/bin/env python3
"""Build a first-pass Chiloglottis duplicated gene atlas.

This script combines the current H1/H2 gene models, MapMan/Mercator
functional annotations, available OrthoFinder results, flower/labellum
expression counts for H1, and high-quality H1/H2 inversions.

It intentionally labels genomic duplication classes as "first_pass":
they are coordinate heuristics within orthogroups, not a replacement for
MCScanX or DupGen_finder mode classification.
"""

import argparse
import csv
import json
import math
import os
import re
from collections import defaultdict


ROOT = "/g/data/xf3/zz3507"

DEFAULTS = {
    "h1_gtf": ROOT + "/compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H1t1.gtf",
    "h2_gtf": ROOT + "/compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H2t1.gtf",
    "h1_func": ROOT + "/Output/20260127Genome/functional_annotation/H1.Functional_Annotation_results.txt",
    "h2_func": ROOT + "/Output/20260127Genome/functional_annotation/H2.Functional_Annotation_results.txt",
    "h1_hog_n0": ROOT + "/database/GeneFamilyExpentionContractionDatasetH12026/OrthoFinder/Results_Feb16/Phylogenetic_Hierarchical_Orthogroups/N0.tsv",
    "h2_orthogroups": ROOT + "/database/GeneFamilyExpentionContractionDatasetH22026/OrthoFinder/Results_Feb15/Orthogroups/Orthogroups.tsv",
    "h2_gene_counts": ROOT + "/database/GeneFamilyExpentionContractionDatasetH22026/OrthoFinder/Results_Feb15/Orthogroups/Orthogroups.GeneCount.tsv",
    "h1_expression": ROOT + "/Output/20260127Genome/H1/inversion_rna_minimap2_stage_sun_flower_L/rna_gene_counts.sun_flower_L.tsv",
    "strict_inversions": ROOT + "/Output/20260127Genome/other_syri_inversion/hic_high_quality_inversions.strict.tsv",
    "relaxed_inversions": ROOT + "/Output/20260127Genome/other_syri_inversion/hic_high_quality_inversions.relaxed.tsv",
    "outdir": ROOT + "/Output/20260127Genome/duplication_atlas",
}

SELECTED_ORCHID_COLUMNS = [
    "Anoectochilus_roxburghii",
    "Apostasia_shenzhenica",
    "Chiloglottis_H2",
    "Cymbidium_ensifolium",
    "Dendrobium_catenatum",
    "Dendrobium_hrysotoxum",
    "Dendrobium_nobile",
    "Gastrodia_elata",
    "Ophrys_sphegodes",
    "Phalaenopsis_equestris",
    "Platanthera_guangdongensis",
    "Platanthera_zijinensis",
    "Vanilla_planifolia",
]

KEYWORDS = {
    "scent_or_specialized_metabolism": [
        "cytochrome p450", "cyp", "methyltransferase", "o-methyltransferase",
        "acyltransferase", "bahd", "acyl-coa", "ketoacyl", "fatty acid",
        "desaturase", "reductase", "dehydrogenase", "oxidase", "monooxygenase",
        "polyketide", "chalcone", "phenylpropanoid", "benzenoid", "terpene",
        "sulfur", "cysteine", "methionine", "lipid", "cuticular lipid",
    ],
    "floral_development_or_regulation": [
        "flower", "floral", "pollen", "petal", "labellum", "mads", "tcp",
        "myb", "bhlh", "yabby", "agamous", "apetala", "sepallata",
        "jasmonate", "auxin", "ethylene", "gibberellin", "transcription factor",
    ],
    "fungal_or_immune_secondary": [
        "chitinase", "lysm", "lectin", "disease resistance", "nbs-lrr",
        "nlr", "receptor-like kinase", "lrr", "immune", "pathogen",
        "fungal", "mycorrh", "symbiosis", "sweet", "phosphate transporter",
        "ammonium transporter", "amino acid transporter", "gras",
    ],
}


def clean_cell(value):
    if value is None:
        return ""
    value = value.strip()
    if len(value) >= 2 and value[0] == "'" and value[-1] == "'":
        value = value[1:-1]
    return value


def base_gene_id(transcript_id):
    return re.sub(r"\.t\d+$", "", transcript_id)


def parse_gene_id_attr(attr):
    attr = attr.strip()
    m = re.search(r'gene_id "([^"]+)"', attr)
    if m:
        return m.group(1)
    if attr and ";" not in attr and " " not in attr:
        return attr
    return attr.split(";")[0].strip()


def parse_gtf_genes(path, haplotype):
    genes = {}
    with open(path, "r", newline="") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            gene_id = parse_gene_id_attr(parts[8])
            start = int(parts[3])
            end = int(parts[4])
            if end < start:
                start, end = end, start
            genes[gene_id] = {
                "haplotype": haplotype,
                "gene_id": gene_id,
                "transcript_id": gene_id + ".t1",
                "prefixed_gene_id": haplotype + "|" + gene_id,
                "prefixed_transcript_id": haplotype + "|" + gene_id + ".t1",
                "chr": parts[0],
                "start": start,
                "end": end,
                "strand": parts[6],
                "length_bp": end - start + 1,
            }
    return genes


def parse_function_annotation(path):
    by_gene = defaultdict(lambda: {
        "bincodes": set(),
        "names": set(),
        "descriptions": set(),
        "types": set(),
    })
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ident = clean_cell(row.get("IDENTIFIER", ""))
            if not ident:
                continue
            gene_id = base_gene_id(ident)
            entry = by_gene[gene_id]
            for out_key, in_key in [
                ("bincodes", "BINCODE"),
                ("names", "NAME"),
                ("descriptions", "DESCRIPTION"),
                ("types", "TYPE"),
            ]:
                val = clean_cell(row.get(in_key, ""))
                if val:
                    entry[out_key].add(val)

    out = {}
    for gene_id, entry in by_gene.items():
        out[gene_id] = {
            "function_bincodes": "; ".join(sorted(entry["bincodes"])),
            "function_names": "; ".join(sorted(entry["names"])),
            "function_descriptions": "; ".join(sorted(entry["descriptions"])),
            "function_types": "; ".join(sorted(entry["types"])),
        }
    return out


def split_gene_list(cell):
    cell = clean_cell(cell)
    if not cell:
        return []
    return [x.strip() for x in cell.split(",") if x.strip()]


def parse_h1_hog_n0(path):
    gene_to_og = {}
    og_to_genes = defaultdict(set)
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "Chiloglottis_H1" not in reader.fieldnames:
            return gene_to_og, {}
        for row in reader:
            og = row.get("OG", "")
            for transcript_id in split_gene_list(row.get("Chiloglottis_H1", "")):
                if re.match(r"^g\d+\.t\d+$", transcript_id):
                    gene_id = base_gene_id(transcript_id)
                    gene_to_og[gene_id] = og
                    og_to_genes[og].add(gene_id)
    counts = {og: len(genes) for og, genes in og_to_genes.items()}
    return gene_to_og, counts


def parse_h2_orthogroups(path):
    gene_to_og = {}
    og_to_genes = defaultdict(set)
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "Chiloglottis_H2" not in reader.fieldnames:
            return gene_to_og, {}
        for row in reader:
            og = row.get("Orthogroup", "")
            for transcript_id in split_gene_list(row.get("Chiloglottis_H2", "")):
                if re.match(r"^g\d+\.t\d+$", transcript_id):
                    gene_id = base_gene_id(transcript_id)
                    gene_to_og[gene_id] = og
                    og_to_genes[og].add(gene_id)
    counts = {og: len(genes) for og, genes in og_to_genes.items()}
    return gene_to_og, counts


def parse_gene_counts(path):
    counts = {}
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            og = row.get("Orthogroup", "")
            if not og:
                continue
            clean = {}
            for key, val in row.items():
                if key in ("Orthogroup", "Total") or val in (None, ""):
                    continue
                try:
                    clean[key] = int(float(val))
                except ValueError:
                    pass
            counts[og] = clean
    return counts


def median(values):
    values = sorted([v for v in values if v is not None])
    if not values:
        return ""
    mid = len(values) // 2
    if len(values) % 2:
        return float(values[mid])
    return (values[mid - 1] + values[mid]) / 2.0


def parse_featurecounts_expression(path):
    rows = []
    sample_names = []
    with open(path, "r", newline="") as handle:
        reader = csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t")
        header = next(reader)
        sample_names = [os.path.basename(x).replace(".sorted.bam", "") for x in header[6:]]
        for row in reader:
            if not row:
                continue
            gene_id = row[0]
            try:
                length = float(row[5])
                counts = [float(x) for x in row[6:]]
            except ValueError:
                continue
            rows.append((gene_id, length, counts))

    rpk_sums = [0.0 for _ in sample_names]
    for _, length, counts in rows:
        if length <= 0:
            continue
        length_kb = length / 1000.0
        for i, count in enumerate(counts):
            rpk_sums[i] += count / length_kb

    expr = {}
    for gene_id, length, counts in rows:
        if length <= 0:
            continue
        length_kb = length / 1000.0
        tpms = []
        for i, count in enumerate(counts):
            rpk = count / length_kb
            tpms.append((rpk / rpk_sums[i] * 1e6) if rpk_sums[i] else 0.0)
        mean_count = sum(counts) / len(counts) if counts else 0.0
        mean_tpm = sum(tpms) / len(tpms) if tpms else 0.0
        expr[gene_id] = {
            "expr_source": "featureCounts_sun_flower_L_H1",
            "expr_n_samples": str(len(sample_names)),
            "mean_raw_count": mean_count,
            "mean_tpm": mean_tpm,
            "max_tpm": max(tpms) if tpms else 0.0,
            "sample_names": ",".join(sample_names),
        }
    return expr


def parse_inversions(path, strictness_label):
    invs = []
    if not os.path.exists(path):
        return invs
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                invs.append({
                    "id": row["ID"],
                    "strictness": strictness_label,
                    "ref_chr": row["RefChr"],
                    "ref_start": int(float(row["RefStart"])),
                    "ref_end": int(float(row["RefEnd"])),
                    "qry_chr": row["QryChr"],
                    "qry_start": int(float(row["QryStart"])),
                    "qry_end": int(float(row["QryEnd"])),
                })
            except (KeyError, ValueError):
                continue
    return invs


def index_inversions(strict_path, relaxed_path):
    merged = {}
    for inv in parse_inversions(relaxed_path, "relaxed"):
        merged[inv["id"]] = inv
    for inv in parse_inversions(strict_path, "strict"):
        inv["strictness"] = "strict"
        merged[inv["id"]] = inv

    by_hap_chr = {"H1": defaultdict(list), "H2": defaultdict(list)}
    for inv in merged.values():
        rstart, rend = sorted([inv["ref_start"], inv["ref_end"]])
        qstart, qend = sorted([inv["qry_start"], inv["qry_end"]])
        by_hap_chr["H1"][inv["ref_chr"]].append((rstart, rend, inv["id"], inv["strictness"]))
        by_hap_chr["H2"][inv["qry_chr"]].append((qstart, qend, inv["id"], inv["strictness"]))
    return by_hap_chr


def find_overlapping_inversions(gene, inv_index):
    overlaps = []
    for start, end, inv_id, strictness in inv_index.get(gene["chr"], []):
        if gene["end"] >= start and gene["start"] <= end:
            overlaps.append((inv_id, strictness))
    if not overlaps:
        return "", ""
    inv_ids = sorted(set(x[0] for x in overlaps))
    strictness = "strict" if any(x[1] == "strict" for x in overlaps) else "relaxed"
    return ";".join(inv_ids), strictness


def keyword_categories(function_text):
    text = function_text.lower()
    hits = []
    for category, kws in KEYWORDS.items():
        if any(kw in text for kw in kws):
            hits.append(category)
    return hits


def classify_duplication(genes_by_key):
    classes = {}
    for key, genes in genes_by_key.items():
        if len(genes) <= 1:
            for gene in genes:
                classes[(gene["haplotype"], gene["gene_id"])] = "singleton_in_current_orthogroup"
            continue
        by_chr = defaultdict(list)
        for gene in genes:
            by_chr[gene["chr"]].append(gene)
        for chrom in by_chr:
            by_chr[chrom].sort(key=lambda x: (x["start"], x["end"]))
        for gene in genes:
            nearest = None
            for other in by_chr[gene["chr"]]:
                if other["gene_id"] == gene["gene_id"] and other["haplotype"] == gene["haplotype"]:
                    continue
                dist = max(0, max(other["start"], gene["start"]) - min(other["end"], gene["end"]))
                nearest = dist if nearest is None else min(nearest, dist)
            if nearest is not None and nearest <= 100000:
                label = "tandem_candidate_first_pass"
            elif nearest is not None and nearest <= 1000000:
                label = "proximal_candidate_first_pass"
            else:
                label = "dispersed_or_segmental_candidate_first_pass"
            classes[(gene["haplotype"], gene["gene_id"])] = label
    return classes


def fmt_num(value, digits=3):
    if value == "" or value is None:
        return ""
    if isinstance(value, int):
        return str(value)
    try:
        value = float(value)
    except (TypeError, ValueError):
        return str(value)
    if math.isnan(value):
        return ""
    return ("{0:." + str(digits) + "f}").format(value)


def build_atlas(args):
    h1_genes = parse_gtf_genes(args.h1_gtf, "H1")
    h2_genes = parse_gtf_genes(args.h2_gtf, "H2")
    h1_func = parse_function_annotation(args.h1_func)
    h2_func = parse_function_annotation(args.h2_func)
    h1_gene_to_og, h1_og_counts = parse_h1_hog_n0(args.h1_hog_n0)
    h2_gene_to_og, h2_og_counts = parse_h2_orthogroups(args.h2_orthogroups)
    h2_counts_all = parse_gene_counts(args.h2_gene_counts)
    h1_expr = parse_featurecounts_expression(args.h1_expression)
    inv_index = index_inversions(args.strict_inversions, args.relaxed_inversions)

    all_rows = []
    genes_by_key = defaultdict(list)

    for hap, genes, funcs, gene_to_og, og_counts, of_source in [
        ("H1", h1_genes, h1_func, h1_gene_to_og, h1_og_counts, "H1_OrthoFinder_N0_HOG"),
        ("H2", h2_genes, h2_func, h2_gene_to_og, h2_og_counts, "H2_OrthoFinder_Orthogroups"),
    ]:
        for gene_id, gene in genes.items():
            og = gene_to_og.get(gene_id, "")
            func = funcs.get(gene_id, {})
            expr = h1_expr.get(gene_id, {}) if hap == "H1" else {}
            inv_ids, inv_strictness = find_overlapping_inversions(gene, inv_index[hap])

            h2_ref_stats = {}
            if hap == "H2" and og in h2_counts_all:
                counts = h2_counts_all[og]
                orchid_values = [
                    counts.get(col) for col in SELECTED_ORCHID_COLUMNS
                    if col in counts and col != "Chiloglottis_H2"
                ]
                orchid_median = median(orchid_values)
                chilo_count = counts.get("Chiloglottis_H2", "")
                if orchid_median != "" and chilo_count != "":
                    expansion_delta = float(chilo_count) - orchid_median
                    expansion_ratio = (float(chilo_count) / orchid_median) if orchid_median else ""
                else:
                    expansion_delta = ""
                    expansion_ratio = ""
                h2_ref_stats = {
                    "selected_orchid_median_copy_count": orchid_median,
                    "h2_vs_orchid_median_delta": expansion_delta,
                    "h2_vs_orchid_median_ratio": expansion_ratio,
                    "ophrys_copy_count": counts.get("Ophrys_sphegodes", ""),
                    "apostasia_copy_count": counts.get("Apostasia_shenzhenica", ""),
                    "gastrodia_copy_count": counts.get("Gastrodia_elata", ""),
                }

            function_text = " ".join([
                func.get("function_names", ""),
                func.get("function_descriptions", ""),
            ])
            categories = keyword_categories(function_text)

            row = dict(gene)
            row.update({
                "orthogroup": og,
                "orthogroup_source": of_source if og else "",
                "copy_count_in_haplotype_orthogroup": og_counts.get(og, "") if og else "",
                "inversion_overlap_ids": inv_ids,
                "inversion_support": inv_strictness,
                "candidate_categories": ";".join(categories),
                "mean_raw_count_sunflower_L": expr.get("mean_raw_count", ""),
                "mean_tpm_sunflower_L": expr.get("mean_tpm", ""),
                "max_tpm_sunflower_L": expr.get("max_tpm", ""),
                "expr_n_samples": expr.get("expr_n_samples", ""),
                "function_bincodes": func.get("function_bincodes", ""),
                "function_names": func.get("function_names", ""),
                "function_descriptions": func.get("function_descriptions", ""),
                "function_types": func.get("function_types", ""),
                "selected_orchid_median_copy_count": h2_ref_stats.get("selected_orchid_median_copy_count", ""),
                "h2_vs_orchid_median_delta": h2_ref_stats.get("h2_vs_orchid_median_delta", ""),
                "h2_vs_orchid_median_ratio": h2_ref_stats.get("h2_vs_orchid_median_ratio", ""),
                "ophrys_copy_count": h2_ref_stats.get("ophrys_copy_count", ""),
                "apostasia_copy_count": h2_ref_stats.get("apostasia_copy_count", ""),
                "gastrodia_copy_count": h2_ref_stats.get("gastrodia_copy_count", ""),
            })
            all_rows.append(row)
            if og:
                genes_by_key[(hap, og)].append(row)

    dup_classes = classify_duplication(genes_by_key)
    for row in all_rows:
        row["duplication_class_first_pass"] = dup_classes.get(
            (row["haplotype"], row["gene_id"]),
            "no_orthogroup_assignment",
        )
        row["candidate_score"] = score_row(row)

    family_summary = summarize_families(all_rows)
    top_candidates = [
        row for row in all_rows
        if row["candidate_score"] >= 3 and (
            row["candidate_categories"] or row["inversion_overlap_ids"]
        )
    ]
    top_candidates.sort(
        key=lambda r: (
            int(r["candidate_score"]),
            float(r["mean_tpm_sunflower_L"] or 0),
            int(r["copy_count_in_haplotype_orthogroup"] or 0),
        ),
        reverse=True,
    )
    h1_expressed_scent_floral = [
        row for row in all_rows
        if row["haplotype"] == "H1"
        and float(row.get("mean_tpm_sunflower_L") or 0) >= args.min_tpm
        and int(row.get("copy_count_in_haplotype_orthogroup") or 0) >= args.min_copy
        and (
            "scent_or_specialized_metabolism" in row.get("candidate_categories", "")
            or "floral_development_or_regulation" in row.get("candidate_categories", "")
        )
    ]
    h1_expressed_scent_floral.sort(
        key=lambda r: (
            int(r["candidate_score"]),
            float(r["mean_tpm_sunflower_L"] or 0),
            int(r["copy_count_in_haplotype_orthogroup"] or 0),
        ),
        reverse=True,
    )
    inversion_candidates = [
        row for row in all_rows
        if row.get("inversion_overlap_ids") and row.get("candidate_categories")
    ]
    inversion_candidates.sort(
        key=lambda r: (
            int(r["candidate_score"]),
            r.get("inversion_support") == "strict",
            int(r.get("copy_count_in_haplotype_orthogroup") or 0),
        ),
        reverse=True,
    )
    candidate_families = [
        row for row in family_summary
        if int(row.get("max_candidate_score") or 0) >= 6 and row.get("candidate_categories")
    ]

    os.makedirs(args.outdir, exist_ok=True)
    write_tsv(os.path.join(args.outdir, "chiloglottis_duplication_gene_atlas.tsv"), all_rows, GENE_COLUMNS)
    write_tsv(os.path.join(args.outdir, "chiloglottis_duplication_family_summary.tsv"), family_summary, FAMILY_COLUMNS)
    write_tsv(os.path.join(args.outdir, "chiloglottis_top_candidate_genes.tsv"), top_candidates[:500], GENE_COLUMNS)
    write_tsv(
        os.path.join(args.outdir, "chiloglottis_H1_expressed_scent_floral_candidates.tsv"),
        h1_expressed_scent_floral,
        GENE_COLUMNS,
    )
    write_tsv(
        os.path.join(args.outdir, "chiloglottis_inversion_linked_candidate_genes.tsv"),
        inversion_candidates,
        GENE_COLUMNS,
    )
    write_tsv(
        os.path.join(args.outdir, "chiloglottis_candidate_families.tsv"),
        candidate_families,
        FAMILY_COLUMNS,
    )
    write_metadata(
        args,
        all_rows,
        family_summary,
        top_candidates,
        h1_expressed_scent_floral,
        inversion_candidates,
        candidate_families,
    )


def score_row(row):
    score = 0
    try:
        copy_count = int(row.get("copy_count_in_haplotype_orthogroup") or 0)
    except ValueError:
        copy_count = 0
    if copy_count >= 5:
        score += 2
    elif copy_count >= 2:
        score += 1

    dup_class = row.get("duplication_class_first_pass", "")
    if dup_class.startswith("tandem") or dup_class.startswith("proximal"):
        score += 1

    if row.get("inversion_support") == "strict":
        score += 2
    elif row.get("inversion_support") == "relaxed":
        score += 1

    try:
        mean_tpm = float(row.get("mean_tpm_sunflower_L") or 0)
    except ValueError:
        mean_tpm = 0
    if mean_tpm >= 50:
        score += 2
    elif mean_tpm >= 10:
        score += 1

    categories = set((row.get("candidate_categories") or "").split(";"))
    if "scent_or_specialized_metabolism" in categories:
        score += 3
    if "floral_development_or_regulation" in categories:
        score += 2
    if "fungal_or_immune_secondary" in categories:
        score += 1

    try:
        delta = float(row.get("h2_vs_orchid_median_delta") or 0)
        ratio = float(row.get("h2_vs_orchid_median_ratio") or 0)
    except ValueError:
        delta = 0
        ratio = 0
    if delta >= 3:
        score += 2
    elif ratio >= 2 and copy_count >= 2:
        score += 1

    return score


def summarize_families(rows):
    grouped = defaultdict(list)
    for row in rows:
        og = row.get("orthogroup", "")
        if og:
            grouped[(row["haplotype"], og)].append(row)

    summaries = []
    for (hap, og), members in grouped.items():
        categories = sorted(set(
            cat
            for row in members
            for cat in (row.get("candidate_categories") or "").split(";")
            if cat
        ))
        inv_ids = sorted(set(
            inv
            for row in members
            for inv in (row.get("inversion_overlap_ids") or "").split(";")
            if inv
        ))
        mean_tpms = [float(row.get("mean_tpm_sunflower_L") or 0) for row in members]
        scores = [int(row.get("candidate_score") or 0) for row in members]
        funcs = sorted(set(row.get("function_names", "") for row in members if row.get("function_names", "")))
        summaries.append({
            "haplotype": hap,
            "orthogroup": og,
            "n_genes": len(members),
            "n_inversion_genes": sum(1 for row in members if row.get("inversion_overlap_ids")),
            "inversion_ids": ";".join(inv_ids),
            "candidate_categories": ";".join(categories),
            "duplication_classes": ";".join(sorted(set(row.get("duplication_class_first_pass", "") for row in members))),
            "max_candidate_score": max(scores) if scores else 0,
            "mean_candidate_score": sum(scores) / len(scores) if scores else 0,
            "mean_tpm_sunflower_L_family": sum(mean_tpms) / len(mean_tpms) if mean_tpms else "",
            "max_tpm_sunflower_L_family": max(mean_tpms) if mean_tpms else "",
            "representative_functions": " | ".join(funcs[:5]),
            "gene_ids": ",".join(sorted(row["gene_id"] for row in members)),
        })
    summaries.sort(key=lambda r: (int(r["max_candidate_score"]), int(r["n_genes"])), reverse=True)
    return summaries


GENE_COLUMNS = [
    "haplotype",
    "gene_id",
    "transcript_id",
    "prefixed_gene_id",
    "chr",
    "start",
    "end",
    "strand",
    "length_bp",
    "orthogroup",
    "orthogroup_source",
    "copy_count_in_haplotype_orthogroup",
    "duplication_class_first_pass",
    "selected_orchid_median_copy_count",
    "h2_vs_orchid_median_delta",
    "h2_vs_orchid_median_ratio",
    "ophrys_copy_count",
    "apostasia_copy_count",
    "gastrodia_copy_count",
    "inversion_overlap_ids",
    "inversion_support",
    "mean_raw_count_sunflower_L",
    "mean_tpm_sunflower_L",
    "max_tpm_sunflower_L",
    "expr_n_samples",
    "candidate_categories",
    "candidate_score",
    "function_bincodes",
    "function_names",
    "function_descriptions",
    "function_types",
]

FAMILY_COLUMNS = [
    "haplotype",
    "orthogroup",
    "n_genes",
    "n_inversion_genes",
    "inversion_ids",
    "candidate_categories",
    "duplication_classes",
    "max_candidate_score",
    "mean_candidate_score",
    "mean_tpm_sunflower_L_family",
    "max_tpm_sunflower_L_family",
    "representative_functions",
    "gene_ids",
]


def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            clean = {}
            for col in columns:
                value = row.get(col, "")
                if isinstance(value, float):
                    clean[col] = fmt_num(value)
                else:
                    clean[col] = value
            writer.writerow(clean)


def write_metadata(
    args,
    atlas_rows,
    family_rows,
    top_rows,
    h1_expressed_scent_floral,
    inversion_candidates,
    candidate_families,
):
    metadata = {
        "script": os.path.abspath(__file__),
        "inputs": {key: getattr(args, key) for key in vars(args) if key != "outdir"},
        "outputs": {
            "gene_atlas": os.path.join(args.outdir, "chiloglottis_duplication_gene_atlas.tsv"),
            "family_summary": os.path.join(args.outdir, "chiloglottis_duplication_family_summary.tsv"),
            "top_candidate_genes": os.path.join(args.outdir, "chiloglottis_top_candidate_genes.tsv"),
            "h1_expressed_scent_floral_candidates": os.path.join(args.outdir, "chiloglottis_H1_expressed_scent_floral_candidates.tsv"),
            "inversion_linked_candidate_genes": os.path.join(args.outdir, "chiloglottis_inversion_linked_candidate_genes.tsv"),
            "candidate_families": os.path.join(args.outdir, "chiloglottis_candidate_families.tsv"),
        },
        "row_counts": {
            "gene_atlas": len(atlas_rows),
            "family_summary": len(family_rows),
            "top_candidate_genes_written_max_500": min(500, len(top_rows)),
            "top_candidate_genes_total_before_cap": len(top_rows),
            "h1_expressed_scent_floral_candidates": len(h1_expressed_scent_floral),
            "inversion_linked_candidate_genes": len(inversion_candidates),
            "candidate_families": len(candidate_families),
        },
        "limitations": [
            "H1 and H2 orthogroups come from separate OrthoFinder runs; OG IDs are not treated as comparable across haplotypes.",
            "duplication_class_first_pass is a coordinate heuristic within each haplotype orthogroup; use MCScanX or DupGen_finder for publication-grade duplication mode calls.",
            "Expression is currently H1-only, derived from the sun_flower_L featureCounts subset and converted to TPM-like values from gene-level counts.",
            "H2 expansion statistics are based on the H2 OrthoFinder GeneCount table; equivalent H1 reference-copy statistics require a standard H1 GeneCount table or a unified OrthoFinder run.",
        ],
    }
    with open(os.path.join(args.outdir, "chiloglottis_duplication_atlas.metadata.json"), "w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    for key, value in DEFAULTS.items():
        parser.add_argument("--" + key.replace("_", "-"), default=value)
    parser.add_argument("--min-tpm", type=float, default=10.0)
    parser.add_argument("--min-copy", type=int, default=2)
    return parser.parse_args()


if __name__ == "__main__":
    build_atlas(parse_args())
