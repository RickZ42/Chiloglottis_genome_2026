#!/usr/bin/env python3
"""
Map RNA expression values to genes overlapping selected inversion candidates.

Inputs:
1) Inversion-gene overlap table (e.g. H1_INV_gene_overlaps.tsv)
2) RNA expression matrix (featureCounts output or generic gene-by-sample table)
3) Inversion IDs (e.g. 1855,1963,1915,1926,1903,1820,1936)

Outputs:
- per-gene expression table for selected inversions
- per-inversion per-sample summary stats
- matrix of mean expression per inversion
- missing-gene report
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import statistics
import sys
from collections import defaultdict
from typing import Dict, Iterable, List, Sequence, Set, Tuple


DEFAULT_INV_IDS = "1855,1963,1915,1926,1903,1820,1936"


def eprint(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, **kwargs)


def infer_delimiter(header_line: str) -> str:
    if "\t" in header_line:
        return "\t"
    if "," in header_line:
        return ","
    return "\t"


def normalize_inv_id(token: str) -> str:
    t = token.strip()
    if not t:
        return ""
    t_up = t.upper()
    if t_up.startswith("INV"):
        suffix = t_up[3:].strip()
        if suffix.isdigit():
            return f"INV{int(suffix)}"
        return t_up
    if t_up.isdigit():
        return f"INV{int(t_up)}"
    m = re.match(r"^INV[_-]?(\d+)$", t_up)
    if m:
        return f"INV{int(m.group(1))}"
    return t_up


def parse_inv_ids(inv_ids_arg: str) -> List[str]:
    tokens: List[str] = []
    if os.path.isfile(inv_ids_arg):
        with open(inv_ids_arg, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                tokens.extend(re.split(r"[\s,;]+", line))
    else:
        tokens = re.split(r"[\s,;]+", inv_ids_arg.strip())

    norm = [normalize_inv_id(x) for x in tokens if x.strip()]
    # Keep order, remove duplicates
    seen: Set[str] = set()
    ordered = []
    for x in norm:
        if x and x not in seen:
            ordered.append(x)
            seen.add(x)
    return ordered


def load_overlap_table(path: str, inv_ids: Set[str]) -> Tuple[List[dict], Set[str], List[str]]:
    with open(path, "r", encoding="utf-8") as f:
        first = f.readline()
        if not first:
            raise ValueError(f"Empty overlap table: {path}")
        delim = infer_delimiter(first)
        f.seek(0)
        reader = csv.DictReader(f, delimiter=delim)
        if reader.fieldnames is None:
            raise ValueError(f"Failed to parse header in overlap table: {path}")

        required = {"INV_ID", "GeneID"}
        missing_cols = required - set(reader.fieldnames)
        if missing_cols:
            raise ValueError(
                f"Missing required columns in overlap table {path}: {sorted(missing_cols)}"
            )

        selected = []
        seen_inv = set()
        for row in reader:
            inv = normalize_inv_id(row.get("INV_ID", ""))
            if not inv:
                continue
            if inv_ids and inv not in inv_ids:
                continue
            row["INV_ID"] = inv
            selected.append(row)
            seen_inv.add(inv)

    return selected, seen_inv, reader.fieldnames


def parse_float(x: str) -> float:
    s = str(x).strip()
    if s == "" or s == "NA" or s == "NaN":
        return math.nan
    try:
        return float(s)
    except ValueError:
        return math.nan


def load_expression_table(path: str, gene_col_arg: str | None) -> Tuple[List[str], Dict[str, List[float]]]:
    # featureCounts files often start with commented metadata lines (starting with '#')
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    data_lines = [ln for ln in lines if not ln.startswith("#") and ln.strip()]
    if not data_lines:
        raise ValueError(f"No data lines found in expression file: {path}")

    delim = infer_delimiter(data_lines[0])
    reader = csv.DictReader(data_lines, delimiter=delim)
    if reader.fieldnames is None:
        raise ValueError(f"Failed to parse header in expression table: {path}")
    header = [h.strip() for h in reader.fieldnames]

    if gene_col_arg:
        if gene_col_arg not in header:
            raise ValueError(f"--gene-col '{gene_col_arg}' not found in expression header.")
        gene_col = gene_col_arg
    else:
        candidates = ["Geneid", "GeneID", "gene_id", "gene", "Gene", "ID"]
        found = [c for c in candidates if c in header]
        gene_col = found[0] if found else header[0]

    # featureCounts metadata columns
    meta_cols = {gene_col, "Chr", "Start", "End", "Strand", "Length"}
    sample_cols = [c for c in header if c not in meta_cols]
    if not sample_cols:
        sample_cols = [c for c in header if c != gene_col]
    if not sample_cols:
        raise ValueError("No sample columns detected in expression table.")

    expr: Dict[str, List[float]] = {}
    duplicates = 0
    for row in reader:
        g = row.get(gene_col, "").strip()
        if not g:
            continue
        vals = [parse_float(row.get(c, "")) for c in sample_cols]
        if g in expr:
            duplicates += 1
            # Keep first occurrence by default.
            continue
        expr[g] = vals

    if duplicates:
        eprint(f"[WARN] {duplicates} duplicated gene IDs in expression table; kept first occurrence.")

    return sample_cols, expr


def nanmean(vals: Sequence[float]) -> float:
    v = [x for x in vals if not math.isnan(x)]
    return float(sum(v) / len(v)) if v else math.nan


def nanmedian(vals: Sequence[float]) -> float:
    v = [x for x in vals if not math.isnan(x)]
    return float(statistics.median(v)) if v else math.nan


def nansum(vals: Sequence[float]) -> float:
    v = [x for x in vals if not math.isnan(x)]
    return float(sum(v)) if v else math.nan


def nanmin(vals: Sequence[float]) -> float:
    v = [x for x in vals if not math.isnan(x)]
    return float(min(v)) if v else math.nan


def nanmax(vals: Sequence[float]) -> float:
    v = [x for x in vals if not math.isnan(x)]
    return float(max(v)) if v else math.nan


def safe_fmt(x: float) -> str:
    return "NA" if (isinstance(x, float) and math.isnan(x)) else f"{x:.6f}"


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Map RNA expression values to genes overlapping inversion candidates."
    )
    ap.add_argument(
        "--inv-overlaps",
        required=True,
        help="Inversion-gene overlap table (e.g. H1_INV_gene_overlaps.tsv).",
    )
    ap.add_argument(
        "--expr",
        required=True,
        help="RNA expression table (featureCounts output or generic gene x sample table).",
    )
    ap.add_argument(
        "--inv-ids",
        default=DEFAULT_INV_IDS,
        help=(
            "Comma/space-separated inversion IDs (supports '1855' or 'INV1855'), "
            f"default: {DEFAULT_INV_IDS}, or provide a file path with IDs."
        ),
    )
    ap.add_argument(
        "--gene-col",
        default=None,
        help="Gene ID column in expression table (auto-detect by default).",
    )
    ap.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix path, e.g. output/inv_rna_map",
    )
    args = ap.parse_args()

    inv_id_list = parse_inv_ids(args.inv_ids)
    if not inv_id_list:
        raise ValueError("No valid inversion IDs provided.")
    requested_set = set(inv_id_list)
    eprint(f"[INFO] Requested inversions ({len(inv_id_list)}): {', '.join(inv_id_list)}")

    overlap_rows, seen_inv, overlap_cols = load_overlap_table(args.inv_overlaps, requested_set)
    if not overlap_rows:
        raise ValueError("No overlap rows found for requested inversion IDs.")

    missing_inv = [x for x in inv_id_list if x not in seen_inv]
    if missing_inv:
        eprint(f"[WARN] Not found in overlap table: {', '.join(missing_inv)}")

    sample_cols, expr = load_expression_table(args.expr, args.gene_col)
    eprint(f"[INFO] Loaded expression for {len(expr)} genes across {len(sample_cols)} samples.")

    # Merge overlaps with expression
    genes_expr_rows: List[dict] = []
    missing_gene_rows: List[dict] = []
    seen_pair: Set[Tuple[str, str]] = set()
    per_inv_gene_values: Dict[str, Dict[str, List[float]]] = defaultdict(dict)

    for row in overlap_rows:
        inv = row["INV_ID"]
        gene = row["GeneID"]
        pair = (inv, gene)
        if pair in seen_pair:
            continue
        seen_pair.add(pair)

        values = expr.get(gene)
        if values is None:
            missing_gene_rows.append(row)
            continue

        out = dict(row)
        for s, v in zip(sample_cols, values):
            out[s] = v
        genes_expr_rows.append(out)
        per_inv_gene_values[inv][gene] = values

    # Prepare output paths
    out_prefix = args.out_prefix
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    out_gene_expr = f"{out_prefix}.genes_expression.tsv"
    out_summary = f"{out_prefix}.per_inversion_sample_summary.tsv"
    out_matrix = f"{out_prefix}.per_inversion_mean_matrix.tsv"
    out_missing_genes = f"{out_prefix}.missing_genes_in_expression.tsv"
    out_inv_meta = f"{out_prefix}.selected_inversions.tsv"

    # Write per-gene expression table
    gene_expr_fields = [c for c in overlap_cols if c in genes_expr_rows[0]] + sample_cols
    with open(out_gene_expr, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=gene_expr_fields, delimiter="\t")
        w.writeheader()
        for row in genes_expr_rows:
            rec = dict(row)
            for s in sample_cols:
                rec[s] = safe_fmt(rec[s])
            w.writerow(rec)

    # Write missing genes table
    with open(out_missing_genes, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=overlap_cols, delimiter="\t")
        w.writeheader()
        for row in missing_gene_rows:
            w.writerow(row)

    # Summaries and matrix
    summary_rows: List[dict] = []
    matrix_rows: List[dict] = []

    for inv in inv_id_list:
        gene_map = per_inv_gene_values.get(inv, {})
        gene_values = list(gene_map.values())

        matrix_row = {"INV_ID": inv, "GenesWithExpression": str(len(gene_values))}
        for si, sample in enumerate(sample_cols):
            vals = [gv[si] for gv in gene_values]
            matrix_row[sample] = safe_fmt(nanmean(vals))

            summary_rows.append(
                {
                    "INV_ID": inv,
                    "Sample": sample,
                    "GenesWithExpression": str(len([x for x in vals if not math.isnan(x)])),
                    "Mean": safe_fmt(nanmean(vals)),
                    "Median": safe_fmt(nanmedian(vals)),
                    "Sum": safe_fmt(nansum(vals)),
                    "Min": safe_fmt(nanmin(vals)),
                    "Max": safe_fmt(nanmax(vals)),
                }
            )
        matrix_rows.append(matrix_row)

    with open(out_summary, "w", encoding="utf-8", newline="") as f:
        cols = ["INV_ID", "Sample", "GenesWithExpression", "Mean", "Median", "Sum", "Min", "Max"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in summary_rows:
            w.writerow(row)

    with open(out_matrix, "w", encoding="utf-8", newline="") as f:
        cols = ["INV_ID", "GenesWithExpression"] + sample_cols
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in matrix_rows:
            w.writerow(row)

    # Selected inversion metadata table
    seen_meta: Set[str] = set()
    inv_meta_rows: List[dict] = []
    meta_cols = [c for c in ["INV_ID", "RefChr", "RefStart", "RefEnd", "QryChr", "QryStart", "QryEnd"] if c in overlap_cols]
    for row in overlap_rows:
        inv = row["INV_ID"]
        if inv in seen_meta:
            continue
        seen_meta.add(inv)
        inv_meta_rows.append({k: row.get(k, "") for k in meta_cols})
    with open(out_inv_meta, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=meta_cols, delimiter="\t")
        w.writeheader()
        for row in inv_meta_rows:
            w.writerow(row)

    eprint("[INFO] Done.")
    eprint(f"[INFO] Gene-level output: {out_gene_expr}")
    eprint(f"[INFO] Per-sample summary: {out_summary}")
    eprint(f"[INFO] Mean matrix: {out_matrix}")
    eprint(f"[INFO] Missing genes: {out_missing_genes}")
    eprint(f"[INFO] Inversion metadata: {out_inv_meta}")


if __name__ == "__main__":
    main()
