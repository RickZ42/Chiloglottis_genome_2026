#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


LABEL_MAP = {
    "Chiloglottis_H1": "Chiloglottis H1",
    "Chiloglottis_H2": "Chiloglottis H2",
    "Anoectochilus_roxburghii_HapA": "A. roxburghii (HapA)",
    "Ophrys_sphegodes": "O. sphegodes",
    "Platanthera_zijinensis": "P. zijinensis",
    "Dactylorhiza_incarnata": "D. incarnata",
}


def read_summary(path: Path):
    rows = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            if r.get("status") != "ok":
                continue
            rows.append(
                {
                    "label": r["label"],
                    "display": LABEL_MAP.get(r["label"], r["label"].replace("_", " ")),
                    "lineage": r["lineage"],
                    "C": float(r["C_pct"]),
                    "S": float(r["S_pct"]),
                    "D": float(r["D_pct"]),
                    "F": float(r["F_pct"]),
                    "M": float(r["M_pct"]),
                }
            )
    return rows


def plot(rows, out_png: Path, out_pdf: Path, title: str | None = None):
    if not rows:
        raise SystemExit("No valid rows (status=ok) found in summary file.")

    # Keep manifest order but show first item at top
    rows_plot = list(reversed(rows))
    y = list(range(len(rows_plot)))

    s_vals = [r["S"] for r in rows_plot]
    d_vals = [r["D"] for r in rows_plot]
    f_vals = [r["F"] for r in rows_plot]
    m_vals = [r["M"] for r in rows_plot]
    labels = [r["display"] for r in rows_plot]
    lineage = rows[0]["lineage"]

    colors = {
        "S": "#5DADE2",
        "D": "#2E86C1",
        "F": "#F4D03F",
        "M": "#EC7063",
    }

    fig_h = max(4.5, 1.0 + 0.78 * len(rows_plot))
    fig, ax = plt.subplots(figsize=(12, fig_h))

    bar_h = 0.82
    ax.barh(y, s_vals, color=colors["S"], edgecolor="white", height=bar_h)
    ax.barh(y, d_vals, left=s_vals, color=colors["D"], edgecolor="white", height=bar_h)
    left_f = [s + d for s, d in zip(s_vals, d_vals)]
    ax.barh(y, f_vals, left=left_f, color=colors["F"], edgecolor="white", height=bar_h)
    left_m = [s + d + f for s, d, f in zip(s_vals, d_vals, f_vals)]
    ax.barh(y, m_vals, left=left_m, color=colors["M"], edgecolor="white", height=bar_h)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlim(0, 100)
    ax.set_xlabel("%BUSCOs", fontsize=11)
    ax.set_title(title or "BUSCO-style Completeness Assessment on Genomes", fontsize=15, weight="bold", pad=18)

    ax.xaxis.grid(True, linestyle="--", linewidth=0.7, alpha=0.35)
    ax.set_axisbelow(True)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(1.8)
    ax.spines["bottom"].set_linewidth(1.8)

    legend_handles = [
        Patch(color=colors["S"], label="Complete single-copy (S)"),
        Patch(color=colors["D"], label="Complete duplicated (D)"),
        Patch(color=colors["F"], label="Fragmented (F)"),
        Patch(color=colors["M"], label="Missing (M)"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.08),
        ncol=2,
        frameon=False,
        fontsize=10,
        handlelength=1.1,
        columnspacing=1.2,
    )

    # Database note on left
    fig.text(0.08, 0.90, f"Database:\n{lineage}", ha="left", va="top", fontsize=10)

    # Inline annotation similar to the example figure
    for yi, r in zip(y, rows_plot):
        txt = f"C:{r['C']:.1f}% [S:{r['S']:.1f}% D:{r['D']:.1f}%] F:{r['F']:.1f}% M:{r['M']:.1f}%"
        ax.text(2.0, yi, txt, va="center", ha="left", fontsize=9.4, color="black")

    plt.tight_layout(rect=(0.12, 0.06, 0.98, 0.94))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Plot BUSCO-style stacked bar chart from compleasm summary TSV.")
    ap.add_argument(
        "-i",
        "--input",
        default="/g/data/xf3/zz3507/Output/20260127Genome/compleasm_orchidoideae_odb12/compleasm_orchidoideae_odb12_summary.tsv",
        help="Input summary TSV from run_orchidoideae_compleasm_odb12.sh",
    )
    ap.add_argument(
        "--out-prefix",
        default="/g/data/xf3/zz3507/Output/20260127Genome/compleasm_orchidoideae_odb12/compleasm_orchidoideae_odb12.busco_style",
        help="Output prefix (without extension)",
    )
    ap.add_argument("--title", default=None, help="Custom plot title")
    args = ap.parse_args()

    rows = read_summary(Path(args.input))
    out_prefix = Path(args.out_prefix)
    out_png = Path(str(out_prefix) + ".png")
    out_pdf = Path(str(out_prefix) + ".pdf")
    plot(rows, out_png, out_pdf, args.title)
    print(out_png)
    print(out_pdf)


if __name__ == "__main__":
    main()
