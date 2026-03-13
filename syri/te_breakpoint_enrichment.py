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
import bisect
import csv
import os
import random
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


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


def parse_repeatmasker_out(path: str) -> List[RepeatEntry]:
    entries: List[RepeatEntry] = []
    with open(path) as f:
        for i, line in enumerate(f):
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
    sizes: Dict[str, int] = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            sizes[parts[0]] = int(parts[1])
    return sizes


def parse_inversions(path: str) -> List[Inversion]:
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


class RepeatIndex:
    """Index repeats by chromosome with binary search on start positions."""

    def __init__(self, repeats: List[RepeatEntry]):
        self._data: Dict[str, List[RepeatEntry]] = defaultdict(list)
        self._starts: Dict[str, List[int]] = {}
        for r in repeats:
            self._data[r.chrom].append(r)
        for chrom in self._data:
            self._data[chrom].sort(key=lambda x: x.start)
            self._starts[chrom] = [r.start for r in self._data[chrom]]

    def query(self, chrom: str, win_start: int, win_end: int) -> List[RepeatEntry]:
        entries = self._data.get(chrom, [])
        if not entries:
            return []
        starts = self._starts.get(chrom, [])
        # Find first entry that could overlap (start <= win_end)
        right = bisect.bisect_right(starts, win_end)
        result = []
        for i in range(right):
            r = entries[i]
            if r.end >= win_start:
                result.append(r)
        return result


def analyze_window(
    ridx: RepeatIndex, chrom: str, win_start: int, win_end: int,
) -> Tuple[float, Dict[str, int], Dict[str, int]]:
    """Return (repeat_frac, class_bp, family_bp) for a window."""
    win_len = win_end - win_start
    if win_len <= 0:
        return 0.0, {}, {}
    overlapping = ridx.query(chrom, win_start, win_end)
    intervals = []
    class_bp: Dict[str, int] = Counter()
    family_bp: Dict[str, int] = Counter()
    for r in overlapping:
        s = max(r.start, win_start)
        e = min(r.end, win_end)
        if s < e:
            intervals.append((s, e))
            bp = e - s
            class_bp[r.repeat_class] += bp
            full_label = (f"{r.repeat_class}/{r.repeat_family}"
                          if r.repeat_class != r.repeat_family
                          else r.repeat_class)
            family_bp[full_label] += bp
    if not intervals:
        return 0.0, {}, {}
    intervals.sort()
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    covered = sum(e - s for s, e in merged)
    return covered / win_len, dict(class_bp), dict(family_bp)


def generate_control_windows(
    chrom_sizes: Dict[str, int],
    n_per_chrom: Dict[str, int],
    window_size: int,
    exclude_regions: List[Tuple[str, int, int]],
    n_controls_per_bp: int = 100,
    seed: int = 42,
) -> List[Tuple[str, int, int]]:
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
            overlap = False
            for es, ee in excl:
                if s < ee and e > es:
                    overlap = True
                    break
            if not overlap:
                controls.append((chrom, s, e))
                found += 1
    return controls


def main():
    parser = argparse.ArgumentParser(
        description="TE enrichment analysis at inversion breakpoints."
    )
    parser.add_argument("--rm-out", required=True)
    parser.add_argument("--inv-tsv", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--window-sizes", default="5000,10000,20000")
    parser.add_argument("--n-controls", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-inv-len", type=int, default=0)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    window_sizes = [int(x) for x in args.window_sizes.split(",")]

    print("[INFO] Parsing RepeatMasker output...", file=sys.stderr)
    repeats = parse_repeatmasker_out(args.rm_out)
    print(f"[INFO] Loaded {len(repeats):,} repeat entries", file=sys.stderr)
    ridx = RepeatIndex(repeats)

    print("[INFO] Parsing inversions...", file=sys.stderr)
    inversions = parse_inversions(args.inv_tsv)
    if args.min_inv_len > 0:
        inversions = [inv for inv in inversions
                      if (inv.ref_end - inv.ref_start) >= args.min_inv_len]
    print(f"[INFO] {len(inversions)} inversions after filtering", file=sys.stderr)

    chrom_sizes = parse_fai(args.fai)

    breakpoints: List[Tuple[str, int, str]] = []
    for inv in inversions:
        breakpoints.append((inv.ref_chr, inv.ref_start, f"{inv.inv_id}_left"))
        breakpoints.append((inv.ref_chr, inv.ref_end, f"{inv.inv_id}_right"))

    exclude_regions = [(inv.ref_chr, inv.ref_start, inv.ref_end) for inv in inversions]

    for ws in window_sizes:
        half = ws // 2
        print(f"\n[INFO] === Window size: +/-{half//1000}kb ({ws} bp) ===", file=sys.stderr)

        # ---- Breakpoint windows ----
        bp_fracs = []
        bp_class_totals: Dict[str, int] = Counter()
        bp_family_totals: Dict[str, int] = Counter()
        bp_total_bp = 0
        bp_records = []

        for chrom, pos, label in breakpoints:
            chrom_len = chrom_sizes.get(chrom, 0)
            win_start = max(1, pos - half)
            win_end = min(chrom_len, pos + half)
            frac, cls_bp, fam_bp = analyze_window(ridx, chrom, win_start, win_end)
            bp_fracs.append(frac)
            for k, v in cls_bp.items():
                bp_class_totals[k] += v
            for k, v in fam_bp.items():
                bp_family_totals[k] += v
            bp_total_bp += (win_end - win_start)
            bp_records.append({
                "label": label, "chrom": chrom, "pos": pos,
                "win_start": win_start, "win_end": win_end,
                "repeat_frac": frac,
            })

        # ---- Control windows ----
        bp_per_chrom: Dict[str, int] = Counter()
        for chrom, pos, label in breakpoints:
            bp_per_chrom[chrom] += 1

        print(f"[INFO] Generating {args.n_controls} control windows per breakpoint...",
              file=sys.stderr)
        controls = generate_control_windows(
            chrom_sizes, bp_per_chrom, ws, exclude_regions,
            n_controls_per_bp=args.n_controls, seed=args.seed,
        )
        print(f"[INFO] Generated {len(controls)} control windows", file=sys.stderr)

        ctrl_fracs = []
        ctrl_class_totals: Dict[str, int] = Counter()
        ctrl_family_totals: Dict[str, int] = Counter()
        ctrl_total_bp = 0

        for i, (chrom, win_start, win_end) in enumerate(controls):
            frac, cls_bp, fam_bp = analyze_window(ridx, chrom, win_start, win_end)
            ctrl_fracs.append(frac)
            for k, v in cls_bp.items():
                ctrl_class_totals[k] += v
            for k, v in fam_bp.items():
                ctrl_family_totals[k] += v
            ctrl_total_bp += (win_end - win_start)
            if (i + 1) % 50000 == 0:
                print(f"[INFO]   processed {i+1}/{len(controls)} control windows...",
                      file=sys.stderr)

        # ---- Stats ----
        stat_mw, pval_mw = stats.mannwhitneyu(bp_fracs, ctrl_fracs, alternative="two-sided")
        bp_mean = np.mean(bp_fracs)
        ctrl_mean = np.mean(ctrl_fracs)
        bp_median = np.median(bp_fracs)
        ctrl_median = np.median(ctrl_fracs)

        print(f"[INFO] Breakpoint repeat fraction: mean={bp_mean:.4f}, median={bp_median:.4f}",
              file=sys.stderr)
        print(f"[INFO] Control repeat fraction:    mean={ctrl_mean:.4f}, median={ctrl_median:.4f}",
              file=sys.stderr)
        print(f"[INFO] Mann-Whitney U: stat={stat_mw:.1f}, p={pval_mw:.2e}", file=sys.stderr)

        # ---- Write per-breakpoint detail ----
        detail_path = os.path.join(args.outdir, f"breakpoint_repeat_detail.ws{ws}.tsv")
        with open(detail_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["label", "chrom", "pos", "win_start", "win_end", "repeat_frac"])
            for r in bp_records:
                w.writerow([r["label"], r["chrom"], r["pos"],
                            r["win_start"], r["win_end"], f"{r['repeat_frac']:.6f}"])

        # ---- Write summary ----
        summary_path = os.path.join(args.outdir, f"breakpoint_vs_control_summary.ws{ws}.tsv")
        with open(summary_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["metric", "breakpoint_windows", "control_windows"])
            w.writerow(["n_windows", len(bp_records), len(controls)])
            w.writerow(["mean_repeat_frac", f"{bp_mean:.6f}", f"{ctrl_mean:.6f}"])
            w.writerow(["median_repeat_frac", f"{bp_median:.6f}", f"{ctrl_median:.6f}"])
            w.writerow(["MannWhitneyU_stat", f"{stat_mw:.1f}", ""])
            w.writerow(["MannWhitneyU_pval", f"{pval_mw:.2e}", ""])

        # ---- TE class composition ----
        all_classes = sorted(set(list(bp_class_totals.keys()) + list(ctrl_class_totals.keys())))
        class_path = os.path.join(args.outdir, f"te_class_composition.ws{ws}.tsv")
        with open(class_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["TE_class", "breakpoint_bp", "breakpoint_frac",
                         "control_bp", "control_frac", "enrichment_ratio"])
            for cls in all_classes:
                bp_v = bp_class_totals.get(cls, 0)
                ct_v = ctrl_class_totals.get(cls, 0)
                bp_f = bp_v / bp_total_bp if bp_total_bp > 0 else 0
                ct_f = ct_v / ctrl_total_bp if ctrl_total_bp > 0 else 0
                enr = (bp_f / ct_f) if ct_f > 0 else float("inf")
                w.writerow([cls, bp_v, f"{bp_f:.6f}", ct_v, f"{ct_f:.6f}", f"{enr:.4f}"])

        # ---- TE family composition ----
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

        # ---- Plots ----

        # 1. Repeat fraction distribution
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
                     f"(window +/-{half//1000} kb, Mann-Whitney p = {pval_mw:.2e})",
                     fontsize=13, fontweight="bold")
        ax.legend(fontsize=10)
        ax.set_xlim(0, 1.05)
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, f"repeat_density_distribution.ws{ws}.png"), dpi=300)
        fig.savefig(os.path.join(args.outdir, f"repeat_density_distribution.ws{ws}.pdf"))
        plt.close(fig)

        # 2. TE class composition barplot
        top_classes = sorted(all_classes,
                             key=lambda c: bp_class_totals.get(c, 0) + ctrl_class_totals.get(c, 0),
                             reverse=True)[:10]
        bp_vals = [bp_class_totals.get(c, 0) / bp_total_bp * 100 for c in top_classes]
        ctrl_vals = [ctrl_class_totals.get(c, 0) / ctrl_total_bp * 100 for c in top_classes]

        x = np.arange(len(top_classes))
        width = 0.35
        fig, ax = plt.subplots(figsize=(10, 5.5))
        ax.bar(x - width / 2, bp_vals, width, label="Breakpoints", color="#E74C3C", alpha=0.85)
        ax.bar(x + width / 2, ctrl_vals, width, label="Control", color="#888888", alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(top_classes, rotation=35, ha="right", fontsize=10)
        ax.set_ylabel("% of window bp", fontsize=12, fontweight="bold")
        ax.set_title(f"TE class composition: breakpoints vs. control (+/-{half//1000} kb)",
                     fontsize=13, fontweight="bold")
        ax.legend(fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, f"te_class_composition.ws{ws}.png"), dpi=300)
        fig.savefig(os.path.join(args.outdir, f"te_class_composition.ws{ws}.pdf"))
        plt.close(fig)

        # 3. TE family composition barplot (top 15)
        bp_fam_vals = [bp_family_totals.get(f, 0) / bp_total_bp * 100 for f in all_families]
        ctrl_fam_vals = [ctrl_family_totals.get(f, 0) / ctrl_total_bp * 100 for f in all_families]

        x = np.arange(len(all_families))
        fig, ax = plt.subplots(figsize=(12, 5.5))
        ax.bar(x - width / 2, bp_fam_vals, width, label="Breakpoints", color="#E74C3C", alpha=0.85)
        ax.bar(x + width / 2, ctrl_fam_vals, width, label="Control", color="#888888", alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(all_families, rotation=40, ha="right", fontsize=9)
        ax.set_ylabel("% of window bp", fontsize=12, fontweight="bold")
        ax.set_title(f"TE family composition: breakpoints vs. control (+/-{half//1000} kb)",
                     fontsize=13, fontweight="bold")
        ax.legend(fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, f"te_family_composition.ws{ws}.png"), dpi=300)
        fig.savefig(os.path.join(args.outdir, f"te_family_composition.ws{ws}.pdf"))
        plt.close(fig)

    print(f"\n[INFO] All outputs written to {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
