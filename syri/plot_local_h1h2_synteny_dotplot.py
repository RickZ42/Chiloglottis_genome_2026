#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_paf(paf_path, min_aln_bp=10000, min_mapq=0):
    recs = []
    qlen = tlen = None
    with open(paf_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qname = f[0]
            ql = int(f[1])
            qs = int(f[2])
            qe = int(f[3])
            strand = f[4]
            tname = f[5]
            tl = int(f[6])
            ts = int(f[7])
            te = int(f[8])
            alnlen = int(f[10])
            mapq = int(f[11])
            if alnlen < min_aln_bp or mapq < min_mapq:
                continue
            if qlen is None:
                qlen = ql
            if tlen is None:
                tlen = tl
            recs.append(
                dict(
                    qname=qname,
                    tname=tname,
                    qlen=ql,
                    tlen=tl,
                    qs=qs,
                    qe=qe,
                    ts=ts,
                    te=te,
                    strand=strand,
                    alnlen=alnlen,
                    mapq=mapq,
                )
            )
    return recs, qlen, tlen


def main():
    ap = argparse.ArgumentParser(description="Plot local H1/H2 synteny dotplot from PAF")
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--title", default="INV1936 local synteny dotplot (H1 vs H2)")
    ap.add_argument("--min-aln-bp", type=int, default=10000)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--x-label", default="H1 local region (bp)")
    ap.add_argument("--y-label", default="H2 local region (bp)")
    ap.add_argument("--h1-inv-start", type=int, required=True, help="H1 inversion start in local coordinates (0-based approx ok)")
    ap.add_argument("--h1-inv-end", type=int, required=True, help="H1 inversion end in local coordinates")
    ap.add_argument("--h2-inv-start", type=int, required=True)
    ap.add_argument("--h2-inv-end", type=int, required=True)
    args = ap.parse_args()

    recs, qlen, tlen = parse_paf(args.paf, min_aln_bp=args.min_aln_bp, min_mapq=args.min_mapq)
    if not recs:
        raise SystemExit("No alignments left after filtering")

    # x = target (H1), y = query (H2)
    fig, ax = plt.subplots(figsize=(8.5, 8.0))

    # light highlight for inversion windows
    ax.axvspan(args.h1_inv_start, args.h1_inv_end, color="#d6eaf8", alpha=0.35, zorder=0)
    ax.axhspan(args.h2_inv_start, args.h2_inv_end, color="#fdebd0", alpha=0.25, zorder=0)

    # boundary lines
    for x in (args.h1_inv_start, args.h1_inv_end):
        ax.axvline(x, color="#1f77b4", linestyle="--", linewidth=1.0, alpha=0.9)
    for y in (args.h2_inv_start, args.h2_inv_end):
        ax.axhline(y, color="#d35400", linestyle="--", linewidth=1.0, alpha=0.9)

    # plot alignments, sorted by length so long blocks are visible
    recs = sorted(recs, key=lambda r: r["alnlen"])
    for r in recs:
        x1, x2 = r["ts"], r["te"]
        if r["strand"] == "+":
            y1, y2 = r["qs"], r["qe"]
            c = "#2E86C1"  # syntenic (+)
        else:
            y1, y2 = r["qe"], r["qs"]  # descending line for inversion (-)
            c = "#C0392B"  # inverted (-)
        lw = 0.5 if r["alnlen"] < 100000 else 1.1
        alpha = 0.35 if r["alnlen"] < 50000 else 0.8
        ax.plot([x1, x2], [y1, y2], color=c, linewidth=lw, alpha=alpha)

    ax.set_xlim(0, tlen)
    ax.set_ylim(0, qlen)
    ax.set_xlabel(args.x_label)
    ax.set_ylabel(args.y_label)
    ax.set_title(args.title, fontsize=12)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.4)

    # Annotation text
    ax.text(
        0.01,
        0.99,
        "Blue: same orientation (+)\nRed: inverted orientation (-)",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.6", alpha=0.85),
    )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    png = Path(str(out_prefix) + ".png")
    pdf = Path(str(out_prefix) + ".pdf")
    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()

