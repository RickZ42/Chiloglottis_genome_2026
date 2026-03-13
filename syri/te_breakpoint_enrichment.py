#!/usr/bin/env python3
"""
TE/repeat enrichment analysis at inversion breakpoints.

Compares repeat class/family composition and density in windows around
inversion breakpoints versus matched random genomic windows (same
chromosome, same size).

Inputs:
  - RepeatMasker .out file (H1 genome)
  - SyRI inversion coordinates TSV
  - H1 genome .fai (chromosome sizes)

Outputs:
  - Breakpoint vs. background TE density comparison (TSV + plot)
  - TE class composition at breakpoints vs. background (TSV + plot)
  - Statistical test results (Fisher's exact, Mann-Whitney U)
"""

from __future__ import annotations

import argparse
import csv
import os
import random
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from scipy import stats


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class RepeatEntry:
    chrom: str
    start: int
    end: int
    repeat_class: str
    repeat_family: str


@dataclass
class Inversion:
    inv_id: str
    ref_chr: str
    ref_start: int
    ref_end: int
    qry_chr: str
    qry_start: int
    qry_end: int


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_repeatmasker_out(path: str) -> List[RepeatEntry]:
    """Parse RepeatMasker .out file into RepeatEntry list."""
    entries: List[RepeatEntry] = []
    with open(path) as f:
        for i, line in enumerate(f):
            # Skip header lines (first 3 lines)
            if i < 3:
                continue
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 11:
                continue
            try:
                chrom = parts[4]
                start = int(parts[5])
                end = int(parts[6])
                class_family = parts[10]
            except (ValueError, IndexError):
                continue
            if "/" in class_family:
                rclass, rfamily = class_family.split("/", 1)
            else:
                rclass = class_family
                rfamily = class_family
            entries.append(RepeatEntry(
                chrom=chrom, start=start, end=end,
                repeat_class=rclass, repeat_family=rfamily,
            ))
    return entries


def parse_fai(path: str) -> Dict[str, int]:
    """Parse .fai to get chromosome sizes."""
    sizes: Dict[str, int] = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            sizes[parts[0]] = int(parts[1])
    return sizes


def parse_inversions(path: str) -> List[Inversion]:
    """Parse inversion coordinates TSV."""
    inversions: List[Inversion] = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            inversions.append(Inversion(
                inv_id=row.get("ID", row.get("id", "")),
                ref_chr=row["RefChr"],
                ref_start=int(row["RefStart"]),
                ref_end=int(row["RefEnd"]),
                qry_chr=row.get("QryChr", ""),
                qry_start=int(row.get("QryStart", 0)),
                qry_end=int(row.get("QryEnd", 0)),
            ))
    return inversions


# ---------------------------------------------------------------------------
# Index repeats by chromosome for fast interval queries
# ---------------------------------------------------------------------------

def index_repeats_by_chrom(
    repeats: List[RepeatEntry],
) -> Dict[str, List[RepeatEntry]]:
    idx: Dict[str, List[RepeatEntry]] = defaultdict(list)
    for r in repeats:
        idx[r.chrom].append(r)
    # Sort by start position for binary search
    for chrom in idx:
        idx[chrom].sort(key=lambda x: x.start)
    return idx


def repeats_in_window(
    repeat_idx: Dict[str, List[RepeatEntry]],
    chrom: str,
    win_start: int,
    win_end: int,
) -> List[RepeatEntry]:
    """Return repeats overlapping a genomic window."""
    entries = repeat_idx.get(chrom, [])
    result = []
    for r in entries:
        if r.end < win_start:
            continue
        if r.start > win_end:
            break
        result.append(r)
    return result


def compute_repeat_coverage(
    overlapping: List[RepeatEntry],
    win_start: int,
    win_end: int,
) -> Tuple[int, float]:
    """Compute total repeat bp covered in window (merging overlaps)."""
    win_len = win_end - win_start
    if win_len <= 0:
        return 0, 0.0
    intervals = []
    for r in overlapping:
        s = max(r.start, win_start)
        e = min(r.end, win_end)
        if s < e:
            intervals.append((s, e))
    if not intervals:
        return 0, 0.0
    # Merge overlapping intervals
    intervals.sort()
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    covered = sum(e - s for s, e in merged)
    return covered, covered / win_len


def compute_class_composition(
    overlapping: List[RepeatEntry],
    win_start: int,
    win_end: int,
) -> Dict[str, int]:
    """Compute bp covered by each repeat class in window."""
    class_bp: Dict[str, int] = Counter()
    for r in overlapping:
        s = max(r.start, win_start)
        e = min(r.end, win_end)
        if s < e:
            class_bp[r.repeat_class] += e - s
    return dict(class_bp)


# ---------------------------------------------------------------------------
# Random control windows
# ---------------------------------------------------------------------------

def generate_control_windows(
    chrom_sizes: Dict[str, int],
    n_per_chrom: Dict[str, int],
    window_size: int,
    exclude_regions: List[Tuple[str, int, int]],
    n_controls_per_bp: int = 100,
    seed: int = 42,
) -> List[Tuple[str, int, int]]:
    """Generate random control windows, matched by chromosome."""
    rng = random.Random(seed)
    exclude_by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for chrom, s, e in exclude_regions:
        exclude_by_chrom[chrom].append((s, e))

    controls: List[Tuple[str, int, int]] = []
    for chrom, n_bp in n_per_chrom.items():
        chrom_len = chrom_sizes.get(chrom, 0)
        if chrom_len < window_size:
            continue
        n_needed = n_bp * n_controls_per_bp
        excl = sorted(exclude_by_chrom.get(chrom, []))
        max_start = chrom_len - window_size
        attempts = 0
        found = 0
        while found < n_needed and attempts < n_needed * 20:
            attempts += 1
            s = rng.randint(1, max_start)
            e = s + window_size
            # Check not overlapping excluded regions
            overlap = False
            for es, ee in excl:
                if s < ee and e > es:
                    overlap = True
                    break
            if not overlap:
                controls.append((chrom, s, e))
                found += 1
    return controls


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="TE enrichment analysis at inversion breakpoints."
    )
    parser.add_argument("--rm-out", required=True,
                        help="RepeatMasker .out file")
    parser.add_argument("--inv-tsv", required=True,
                        help="Inversion coordinates TSV (with RefChr, RefStart, RefEnd, ID)")
    parser.add_argument("--fai", required=True,
                        help="Genome .fai file for chromosome sizes")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    parser.add_argument("--window-sizes", default="5000,10000,20000",
                        help="Comma-separated breakpoint window sizes (bp) (default: 5000,10000,20000)")
    parser.add_argument("--n-controls", type=int, default=1000,
                        help="Number of random control windows per breakpoint (default: 1000)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--min-inv-len", type=int, default=0,
                        help="Minimum inversion length to include (default: 0, all)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    window_sizes = [int(x) for x in args.window_sizes.split(",")]

    print("[INFO] Parsing RepeatMasker output...", file=sys.stderr)
    repeats = parse_repeatmasker_out(args.rm_out)
    print(f"[INFO] Loaded {len(repeats):,} repeat entries", file=sys.stderr)
    repeat_idx = index_repeats_by_chrom(repeats)

    print("[INFO] Parsing inversions...", file=sys.stderr)
    inversions = parse_inversions(args.inv_tsv)
    if args.min_inv_len > 0:
        inversions = [inv for inv in inversions
                      if (inv.ref_end - inv.ref_start) >= args.min_inv_len]
    print(f"[INFO] {len(inversions)} inversions after filtering", file=sys.stderr)

    chrom_sizes = parse_fai(args.fai)

    # Build breakpoint list (two breakpoints per inversion on H1 ref)
    breakpoints: List[Tuple[str, int, str]] = []  # (chrom, pos, label)
    for inv in inversions:
        breakpoints.append((inv.ref_chr, inv.ref_start, f"{inv.inv_id}_left"))
        breakpoints.append((inv.ref_chr, inv.ref_end, f"{inv.inv_id}_right"))

    # Exclude inversion interior regions from control sampling
    exclude_regions = [(inv.ref_chr, inv.ref_start, inv.ref_end) for inv in inversions]

    # For each window size, compute breakpoint vs. control TE stats
    for ws in window_sizes:
        half = ws // 2
        print(f"\n[INFO] === Window size: ±{half//1000}kb ({ws} bp) ===", file=sys.stderr)

        # Breakpoint windows
        bp_records = []
        bp_class_totals: Dict[str, int] = Counter()
        bp_total_bp = 0

        for chrom, pos, label in breakpoints:
            chrom_len = chrom_sizes.get(chrom, 0)
            win_start = max(1, pos - half)
            win_end = min(chrom_len, pos + half)
            overlapping = repeats_in_window(repeat_idx, chrom, win_start, win_end)
            covered_bp, frac = compute_repeat_coverage(overlapping, win_start, win_end)
            class_comp = compute_class_composition(overlapping, win_start, win_end)

            bp_records.append({
                "label": label,
                "chrom": chrom,
                "pos": pos,
                "win_start": win_start,
                "win_end": win_end,
                "repeat_bp": covered_bp,
                "repeat_frac": frac,
                "class_comp": class_comp,
            })
            for cls, bp_count in class_comp.items():
                bp_class_totals[cls] += bp_count
            bp_total_bp += (win_end - win_start)

        bp_fracs = [r["repeat_frac"] for r in bp_records]

        # Control windows
        bp_per_chrom: Dict[str, int] = Counter()
        for chrom, pos, label in breakpoints:
            bp_per_chrom[chrom] += 1

        print(f"[INFO] Generating {args.n_controls} control windows per breakpoint...",
              file=sys.stderr)
        total_controls_needed = len(breakpoints) * args.n_controls
        controls = generate_control_windows(
            chrom_sizes, bp_per_chrom, ws, exclude_regions,
            n_controls_per_bp=args.n_controls, seed=args.seed,
        )
        print(f"[INFO] Generated {len(controls)} control windows", file=sys.stderr)

        ctrl_fracs = []
        ctrl_class_totals: Dict[str, int] = Counter()
        ctrl_total_bp = 0

        for chrom, win_start, win_end in controls:
            overlapping = repeats_in_window(repeat_idx, chrom, win_start, win_end)
            covered_bp, frac = compute_repeat_coverage(overlapping, win_start, win_end)
            class_comp = compute_class_composition(overlapping, win_start, win_end)
            ctrl_fracs.append(frac)
            for cls, bp_count in class_comp.items():
                ctrl_class_totals[cls] += bp_count
            ctrl_total_bp += (win_end - win_start)

        # Statistical tests
        stat_mw, pval_mw = stats.mannwhitneyu(
            bp_fracs, ctrl_fracs, alternative="two-sided"
        )

        bp_mean = np.mean(bp_fracs)
        ctrl_mean = np.mean(ctrl_fracs)
        bp_median = np.median(bp_fracs)
        ctrl_median = np.median(ctrl_fracs)

        print(f"[INFO] Breakpoint repeat fraction: mean={bp_mean:.4f}, median={bp_median:.4f}",
              file=sys.stderr)
        print(f"[INFO] Control repeat fraction:    mean={ctrl_mean:.4f}, median={ctrl_median:.4f}",
              file=sys.stderr)
        print(f"[INFO] Mann-Whitney U: stat={stat_mw:.1f}, p={pval_mw:.2e}", file=sys.stderr)

        # Write per-breakpoint detail
        detail_path = os.path.join(args.outdir, f"breakpoint_repeat_detail.ws{ws}.tsv")
        with open(detail_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["label", "chrom", "pos", "win_start", "win_end",
                         "repeat_bp", "repeat_frac"])
            for r in bp_records:
                w.writerow([r["label"], r["chrom"], r["pos"], r["win_start"],
                            r["win_end"], r["repeat_bp"], f"{r['repeat_frac']:.6f}"])

        # Write summary stats
        summary_path = os.path.join(args.outdir, f"breakpoint_vs_control_summary.ws{ws}.tsv")
        with open(summary_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["metric", "breakpoint_windows", "control_windows"])
            w.writerow(["n_windows", len(bp_records), len(controls)])
            w.writerow(["mean_repeat_frac", f"{bp_mean:.6f}", f"{ctrl_mean:.6f}"])
            w.writerow(["median_repeat_frac", f"{bp_median:.6f}", f"{ctrl_median:.6f}"])
            w.writerow(["MannWhitneyU_stat", f"{stat_mw:.1f}", ""])
            w.writerow(["MannWhitneyU_pval", f"{pval_mw:.2e}", ""])

        # TE class composition comparison
        all_classes = sorted(set(list(bp_class_totals.keys()) + list(ctrl_class_totals.keys())))
        class_path = os.path.join(args.outdir, f"te_class_composition.ws{ws}.tsv")
        with open(class_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["TE_class", "breakpoint_bp", "breakpoint_frac",
                         "control_bp", "control_frac", "enrichment_ratio",
                         "fisher_pval"])
            for cls in all_classes:
                bp_cls_bp = bp_class_totals.get(cls, 0)
                ctrl_cls_bp = ctrl_class_totals.get(cls, 0)
                bp_frac_cls = bp_cls_bp / bp_total_bp if bp_total_bp > 0 else 0
                ctrl_frac_cls = ctrl_cls_bp / ctrl_total_bp if ctrl_total_bp > 0 else 0
                enrichment = (bp_frac_cls / ctrl_frac_cls) if ctrl_frac_cls > 0 else float("inf")

                # Fisher's exact test: TE-class bp vs non-TE bp in breakpoint vs control
                a = bp_cls_bp
                b = bp_total_bp - bp_cls_bp
                c = ctrl_cls_bp
                d = ctrl_total_bp - ctrl_cls_bp
                # Scale down large numbers for Fisher's exact
                scale = max(1, (a + b + c + d) // 1_000_000)
                _, fisher_p = stats.fisher_exact([
                    [max(1, a // scale), max(1, b // scale)],
                    [max(1, c // scale), max(1, d // scale)],
                ], alternative="two-sided")

                w.writerow([cls, bp_cls_bp, f"{bp_frac_cls:.6f}",
                            ctrl_cls_bp, f"{ctrl_frac_cls:.6f}",
                            f"{enrichment:.4f}", f"{fisher_p:.2e}"])

        # --------------- Plots ---------------

        # 1. Repeat fraction distribution: breakpoints vs. controls
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(ctrl_fracs, bins=50, alpha=0.6, label=f"Control (n={len(controls):,})",
                color="#888888", density=True)
        ax.hist(bp_fracs, bins=20, alpha=0.75, label=f"Breakpoints (n={len(bp_records)})",
                color="#E74C3C", density=True)
        ax.axvline(bp_median, color="#E74C3C", linestyle="--", linewidth=1.5,
                   label=f"BP median={bp_median:.3f}")
        ax.axvline(ctrl_median, color="#555555", linestyle="--", linewidth=1.5,
                   label=f"Ctrl median={ctrl_median:.3f}")
        ax.set_xlabel("Repeat fraction in window", fontsize=12, fontweight="bold")
        ax.set_ylabel("Density", fontsize=12, fontweight="bold")
        ax.set_title(f"Repeat density at inversion breakpoints vs. random control\n"
                     f"(window ±{half//1000} kb, Mann-Whitney p = {pval_mw:.2e})",
                     fontsize=13, fontweight="bold")
        ax.legend(fontsize=10)
        ax.set_xlim(0, 1.05)
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, f"repeat_density_distribution.ws{ws}.png"),
                    dpi=300)
        fig.savefig(os.path.join(args.outdir, f"repeat_density_distribution.ws{ws}.pdf"))
        plt.close(fig)

        # 2. TE class composition barplot
        # Top classes by total bp
        top_classes = sorted(all_classes,
                             key=lambda c: bp_class_totals.get(c, 0) + ctrl_class_totals.get(c, 0),
                             reverse=True)[:10]
        bp_vals = [bp_class_totals.get(c, 0) / bp_total_bp * 100 if bp_total_bp > 0 else 0
                   for c in top_classes]
        ctrl_vals = [ctrl_class_totals.get(c, 0) / ctrl_total_bp * 100 if ctrl_total_bp > 0 else 0
                     for c in top_classes]

        x = np.arange(len(top_classes))
        width = 0.35
        fig, ax = plt.subplots(figsize=(10, 5.5))
        ax.bar(x - width / 2, bp_vals, width, label="Breakpoints", color="#E74C3C", alpha=0.85)
        ax.bar(x + width / 2, ctrl_vals, width, label="Control", color="#888888", alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(top_classes, rotation=35, ha="right", fontsize=10)
        ax.set_ylabel("% of window bp", fontsize=12, fontweight="bold")
        ax.set_title(f"TE class composition: breakpoints vs. control (±{half//1000} kb)",
                     fontsize=13, fontweight="bold")
        ax.legend(fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, f"te_class_composition.ws{ws}.png"), dpi=300)
        fig.savefig(os.path.join(args.outdir, f"te_class_composition.ws{ws}.pdf"))
        plt.close(fig)

        # 3. TE family-level composition (top 15 families)
        bp_family_totals: Dict[str, int] = Counter()
        ctrl_family_totals: Dict[str, int] = Counter()

        for r in bp_records:
            overlapping = repeats_in_window(repeat_idx, r["chrom"], r["win_start"], r["win_end"])
            for rep in overlapping:
                s = max(rep.start, r["win_start"])
                e = min(rep.end, r["win_end"])
                if s < e:
                    full_label = f"{rep.repeat_class}/{rep.repeat_family}" if rep.repeat_class != rep.repeat_family else rep.repeat_class
                    bp_family_totals[full_label] += e - s

        for chrom, win_start, win_end in controls:
            overlapping = repeats_in_window(repeat_idx, chrom, win_start, win_end)
            for rep in overlapping:
                s = max(rep.start, win_start)
                e = min(rep.end, win_end)
                if s < e:
                    full_label = f"{rep.repeat_class}/{rep.repeat_family}" if rep.repeat_class != rep.repeat_family else rep.repeat_class
                    ctrl_family_totals[full_label] += e - s

        all_families = sorted(
            set(list(bp_family_totals.keys()) + list(ctrl_family_totals.keys())),
            key=lambda f: bp_family_totals.get(f, 0) + ctrl_family_totals.get(f, 0),
            reverse=True,
        )[:15]

        family_path = os.path.join(args.outdir, f"te_family_composition.ws{ws}.tsv")
        with open(family_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["TE_family", "breakpoint_bp", "breakpoint_pct",
                         "control_bp", "control_pct", "enrichment_ratio"])
            for fam in all_families:
                bp_v = bp_family_totals.get(fam, 0)
                ct_v = ctrl_family_totals.get(fam, 0)
                bp_pct = bp_v / bp_total_bp * 100 if bp_total_bp > 0 else 0
                ct_pct = ct_v / ctrl_total_bp * 100 if ctrl_total_bp > 0 else 0
                enr = (bp_pct / ct_pct) if ct_pct > 0 else float("inf")
                w.writerow([fam, bp_v, f"{bp_pct:.4f}", ct_v, f"{ct_pct:.4f}", f"{enr:.4f}"])

        bp_fam_vals = [bp_family_totals.get(f, 0) / bp_total_bp * 100 if bp_total_bp > 0 else 0
                       for f in all_families]
        ctrl_fam_vals = [ctrl_family_totals.get(f, 0) / ctrl_total_bp * 100 if ctrl_total_bp > 0 else 0
                         for f in all_families]

        x = np.arange(len(all_families))
        fig, ax = plt.subplots(figsize=(12, 5.5))
        ax.bar(x - width / 2, bp_fam_vals, width, label="Breakpoints", color="#E74C3C", alpha=0.85)
        ax.bar(x + width / 2, ctrl_fam_vals, width, label="Control", color="#888888", alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(all_families, rotation=40, ha="right", fontsize=9)
        ax.set_ylabel("% of window bp", fontsize=12, fontweight="bold")
        ax.set_title(f"TE family composition: breakpoints vs. control (±{half//1000} kb)",
                     fontsize=13, fontweight="bold")
        ax.legend(fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, f"te_family_composition.ws{ws}.png"), dpi=300)
        fig.savefig(os.path.join(args.outdir, f"te_family_composition.ws{ws}.pdf"))
        plt.close(fig)

    print(f"\n[INFO] All outputs written to {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
