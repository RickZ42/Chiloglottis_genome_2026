#!/usr/bin/env python3
"""
Combine inversion Hi-C summary tables and RNA inversion expression summaries into one XLSX workbook.

Default behavior:
- Combine Hi-C summaries from joint + other folders (171 inversions total if present)
- Add strict HQ flag from the strict_HQ Hi-C subset
- Merge RNA overview stats (52 inversions with gene overlaps / expression)
- Optionally add overlap-gene count summary
- Write multiple sheets for downstream manual review
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd


def shorten_sample_name(x: str) -> str:
    if not isinstance(x, str):
        return x
    base = os.path.basename(x)
    if base.endswith(".sorted.bam"):
        return base[:-11]
    if base.endswith(".bam"):
        return base[:-4]
    return base


def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)


def normalize_inv_id_col(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "INV_ID" in out.columns:
        out["INV_ID"] = out["INV_ID"].astype(str)
        return out
    if "ID" in out.columns:
        out = out.rename(columns={"ID": "INV_ID"})
        out["INV_ID"] = out["INV_ID"].astype(str)
        return out
    raise ValueError("No INV_ID/ID column found.")


def clean_rna_mean_matrix(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    rename_map = {}
    for c in out.columns:
        if c in {"INV_ID", "GenesWithExpression"}:
            continue
        rename_map[c] = shorten_sample_name(c)
    return out.rename(columns=rename_map)


def clean_rna_sample_summary(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "Sample" in out.columns:
        out["Sample"] = out["Sample"].map(shorten_sample_name)
    return out


def build_hic_all(joint_path: str, other_path: str, strict_hq_path: str | None) -> pd.DataFrame:
    dfs = []
    for src, p in [("joint", joint_path), ("other", other_path)]:
        if not p:
            continue
        if not Path(p).exists():
            raise FileNotFoundError(p)
        df = normalize_inv_id_col(read_tsv(p))
        df["HiC_source"] = src
        dfs.append(df)

    if not dfs:
        raise ValueError("No Hi-C summary tables provided.")

    hic = pd.concat(dfs, ignore_index=True)
    hic = hic.drop_duplicates(subset=["INV_ID"], keep="first").copy()

    strict_ids = set()
    if strict_hq_path and Path(strict_hq_path).exists():
        strict_df = normalize_inv_id_col(read_tsv(strict_hq_path))
        strict_ids = set(strict_df["INV_ID"].astype(str))
    hic["HiC_strict_HQ"] = hic["INV_ID"].astype(str).isin(strict_ids)

    return hic


def build_readme_rows(hic: pd.DataFrame, rna_overview: pd.DataFrame, merged: pd.DataFrame) -> pd.DataFrame:
    rna_id_set = set(rna_overview["INV_ID"].astype(str))
    merged_rna_hits = int(merged["RNA_has_expression_summary"].fillna(False).astype(bool).sum())
    rows = [
        ["HiC inversions (combined)", len(hic)],
        ["RNA inversions with gene-expression summary", len(rna_overview)],
        ["Merged rows", len(merged)],
        ["Merged rows with RNA summary", merged_rna_hits],
        ["HiC-only rows (no RNA summary)", len(merged) - merged_rna_hits],
        ["RNA IDs missing from Hi-C combined table", int((~rna_overview["INV_ID"].isin(hic["INV_ID"])).sum())],
        ["Note", "RNA summary only exists for inversions with overlapping annotated genes"],
    ]
    return pd.DataFrame(rows, columns=["Field", "Value"])


def main() -> None:
    ap = argparse.ArgumentParser(description="Combine inversion Hi-C and RNA summaries into one XLSX file.")
    ap.add_argument("--hic-joint", required=True, help="Hi-C summary TSV for joint set (e.g. 6 inversions).")
    ap.add_argument("--hic-other", required=True, help="Hi-C summary TSV for other inversions.")
    ap.add_argument("--hic-strict", default=None, help="Hi-C summary TSV for strict HQ subset (optional).")
    ap.add_argument("--rna-overview", required=True, help="RNA all-inversion overview TSV.")
    ap.add_argument("--rna-mean-matrix", required=True, help="RNA per-inversion mean matrix TSV.")
    ap.add_argument("--rna-sample-summary", required=True, help="RNA per-inversion sample summary TSV.")
    ap.add_argument(
        "--inv-gene-overlap-summary",
        default=None,
        help="Optional H1_INV_gene_overlap_summary.tsv to add OverlappingGeneCount.",
    )
    ap.add_argument("--out-xlsx", required=True, help="Output XLSX path.")
    args = ap.parse_args()

    hic = build_hic_all(args.hic_joint, args.hic_other, args.hic_strict)

    rna_overview = normalize_inv_id_col(read_tsv(args.rna_overview))
    rna_mean = normalize_inv_id_col(read_tsv(args.rna_mean_matrix))
    rna_sample = normalize_inv_id_col(read_tsv(args.rna_sample_summary))

    # Prefix RNA overview columns before merge to avoid clashes with shared coordinate columns.
    rna_overview_pref = rna_overview.copy()
    rna_overview_pref = rna_overview_pref.rename(
        columns={c: f"RNA_{c}" for c in rna_overview_pref.columns if c != "INV_ID"}
    )
    rna_overview_pref["RNA_has_expression_summary"] = True

    merged = hic.merge(rna_overview_pref, on="INV_ID", how="left")
    if "RNA_has_expression_summary" not in merged.columns:
        merged["RNA_has_expression_summary"] = False
    merged["RNA_has_expression_summary"] = merged["RNA_has_expression_summary"].fillna(False)

    # Optional overlap gene count summary table (171 inversions, including zero-gene overlaps)
    overlap_summary = None
    if args.inv_gene_overlap_summary and Path(args.inv_gene_overlap_summary).exists():
        overlap_summary = normalize_inv_id_col(read_tsv(args.inv_gene_overlap_summary))
        # Avoid coordinate duplication from overlap summary; keep count and maybe coords if missing.
        keep_cols = ["INV_ID"]
        for c in ["OverlappingGeneCount"]:
            if c in overlap_summary.columns:
                keep_cols.append(c)
        overlap_summary = overlap_summary[keep_cols]
        merged = merged.merge(overlap_summary, on="INV_ID", how="left")

    # Reorder merged columns to keep key fields up front.
    front = ["INV_ID", "HiC_source", "HiC_strict_HQ"]
    for c in ["RefChr", "RefStart", "RefEnd", "QryChr", "QryStart", "QryEnd", "MinLen"]:
        if c in merged.columns:
            front.append(c)
    if "OverlappingGeneCount" in merged.columns:
        front.append("OverlappingGeneCount")
    if "RNA_has_expression_summary" in merged.columns:
        front.append("RNA_has_expression_summary")
    rest = [c for c in merged.columns if c not in front]
    merged = merged[front + rest]

    rna_mean_clean = clean_rna_mean_matrix(rna_mean)
    rna_sample_clean = clean_rna_sample_summary(rna_sample)
    readme = build_readme_rows(hic, rna_overview, merged)

    out_path = Path(args.out_xlsx)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(out_path, engine="xlsxwriter") as xw:
        merged.to_excel(xw, sheet_name="Merged_HiC_RNA", index=False)
        hic.to_excel(xw, sheet_name="HiC_All", index=False)
        if overlap_summary is not None:
            overlap_summary.to_excel(xw, sheet_name="INV_GeneOverlap", index=False)
        rna_overview.to_excel(xw, sheet_name="RNA_Overview", index=False)
        rna_mean_clean.to_excel(xw, sheet_name="RNA_MeanMatrix", index=False)
        rna_sample_clean.to_excel(xw, sheet_name="RNA_SampleSummary", index=False)
        readme.to_excel(xw, sheet_name="README", index=False)

        # Light formatting for readability.
        wb = xw.book
        header_fmt = wb.add_format({"bold": True, "bg_color": "#DCE6F1", "border": 1})
        for sheet_name, df in [
            ("Merged_HiC_RNA", merged),
            ("HiC_All", hic),
            ("INV_GeneOverlap", overlap_summary if overlap_summary is not None else None),
            ("RNA_Overview", rna_overview),
            ("RNA_MeanMatrix", rna_mean_clean),
            ("RNA_SampleSummary", rna_sample_clean),
            ("README", readme),
        ]:
            if df is None:
                continue
            ws = xw.sheets[sheet_name]
            ws.freeze_panes(1, 0)
            for col_idx, col_name in enumerate(df.columns):
                ws.write(0, col_idx, col_name, header_fmt)
                width = min(max(len(str(col_name)) + 2, 12), 28)
                if sheet_name == "RNA_MeanMatrix":
                    width = min(max(len(str(col_name)) + 2, 12), 20)
                ws.set_column(col_idx, col_idx, width)

    print(f"[INFO] Wrote workbook: {out_path}")
    print(f"[INFO] HiC rows: {len(hic)}")
    print(f"[INFO] RNA overview rows: {len(rna_overview)}")
    print(f"[INFO] Merged rows: {len(merged)}")
    print(
        "[INFO] Merged rows with RNA summary: "
        f"{int(merged['RNA_has_expression_summary'].fillna(False).astype(bool).sum())}"
    )


if __name__ == "__main__":
    main()
