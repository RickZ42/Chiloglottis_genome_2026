#!/usr/bin/env python3
"""
Subset a featureCounts-style matrix to selected sample columns.

The script preserves the standard featureCounts metadata columns:
Geneid, Chr, Start, End, Strand, Length

Selected samples can be provided as a comma-separated list or a text file
containing one sample per line. Sample names are matched flexibly against:
1) the original header cell
2) the basename of the header cell
3) the basename with trailing .bam / .sorted.bam removed
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Dict, List, Sequence


META_COLS = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def infer_delimiter(line: str) -> str:
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return "\t"


def load_requested_samples(value: str) -> List[str]:
    if os.path.isfile(value):
        out: List[str] = []
        with open(value, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                out.append(line)
        return out
    return [x.strip() for x in value.split(",") if x.strip()]


def sample_aliases(header_name: str) -> List[str]:
    aliases = [header_name]
    base = os.path.basename(header_name)
    aliases.append(base)
    for suffix in (".sorted.bam", ".bam"):
        if base.endswith(suffix):
            aliases.append(base[: -len(suffix)])
    return aliases


def build_sample_map(header: Sequence[str]) -> Dict[str, str]:
    sample_map: Dict[str, str] = {}
    for col in header:
        if col in META_COLS:
            continue
        for alias in sample_aliases(col):
            sample_map.setdefault(alias, col)
    return sample_map


def main() -> None:
    parser = argparse.ArgumentParser(description="Subset featureCounts matrix by sample.")
    parser.add_argument("--input", required=True, help="Input featureCounts-style TSV.")
    parser.add_argument(
        "--samples",
        required=True,
        help="Comma-separated samples or a file with one sample per line.",
    )
    parser.add_argument("--output", required=True, help="Output TSV path.")
    args = parser.parse_args()

    requested = load_requested_samples(args.samples)
    if not requested:
        raise SystemExit("No samples requested.")

    with open(args.input, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    comment_lines = [line for line in lines if line.startswith("#")]
    data_lines = [line for line in lines if line.strip() and not line.startswith("#")]
    if not data_lines:
        raise SystemExit(f"No table data found in {args.input}")

    delim = infer_delimiter(data_lines[0])
    reader = csv.DictReader(data_lines, delimiter=delim)
    if reader.fieldnames is None:
        raise SystemExit(f"Failed to parse header from {args.input}")

    header = list(reader.fieldnames)
    missing_meta = [c for c in META_COLS if c not in header]
    if missing_meta:
        raise SystemExit(f"Missing featureCounts metadata columns: {missing_meta}")

    sample_map = build_sample_map(header)
    selected_actual: List[str] = []
    missing_requested: List[str] = []
    for sample in requested:
        actual = sample_map.get(sample)
        if actual is None:
            missing_requested.append(sample)
            continue
        if actual not in selected_actual:
            selected_actual.append(actual)

    if missing_requested:
        raise SystemExit(f"Requested samples not found: {', '.join(missing_requested)}")
    if not selected_actual:
        raise SystemExit("No sample columns matched the request.")

    out_header = META_COLS + selected_actual
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        for line in comment_lines:
            handle.write(line)
        writer = csv.DictWriter(handle, fieldnames=out_header, delimiter="\t")
        writer.writeheader()
        for row in reader:
            writer.writerow({col: row.get(col, "") for col in out_header})

    eprint(f"[INFO] Wrote subset matrix: {args.output}")
    eprint(f"[INFO] Selected samples ({len(selected_actual)}):")
    for col in selected_actual:
        eprint(f"  - {col}")


if __name__ == "__main__":
    main()
