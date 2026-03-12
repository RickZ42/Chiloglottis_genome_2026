#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Create publication-grade inter-haplotype Hi-C context maps using Juicer-compatible processing.

This script:
1) Builds a joint H1+H2 reference with prefixed chromosome IDs.
2) Maps Hi-C reads with BWA-MEM (`-SP5M`) as in Juicer single-CPU workflow.
3) Runs Juicer core text processing (`chimeric_blacklist.awk`, fragment assignment, `dups.awk`).
4) Builds normalized `.hic` with `juicer_tools pre`.
5) Extracts local H1-vs-H2 windows with `juicer_tools dump observed` and plots standardized panels.

Usage:
  bash run_joint_hic_inversion_context_map.sh [options]

Options:
  --h1-fa PATH              H1 FASTA
  --h2-fa PATH              H2 FASTA
  --inv-tsv PATH            inversion TSV (required: ID,RefChr,RefStart,RefEnd,QryChr,QryStart,QryEnd; optional: H1WindowStart,H1WindowEnd,H2WindowStart,H2WindowEnd)
  --r1 PATH                 Hi-C read1 FASTQ(.gz)
  --r2 PATH                 Hi-C read2 FASTQ(.gz)
  --outdir PATH             output directory
  --threads INT             threads for BWA/sort (default: 32)
  --mapq INT                MAPQ cutoff used by juicer_tools pre (default: 0; equivalent to retaining all signals)
  --bin-size INT            matrix bin size in bp (default: 100000)
  --flank-bp INT            flank size around inversion in bp (default: 3000000)
  --min-inv-len INT         min inversion length to keep (default: 500000)
  --max-inversions INT      max number of inversions to plot (default: 8)
  --priority-scaffolds CSV  prioritize these scaffolds (default: scaffold_5,scaffold_10,scaffold_19,scaffold_20)
  --inv-ids CSV             explicit inversion IDs to plot (default: auto)
  --norm STR                matrix normalization for dump (default: KR)
  --fallback-norm STR       fallback normalization when --norm unavailable (default: VC_SQRT)
  --skip-juicer-core        skip mapping/dedup/pre and only re-dump/replot from existing .hic
  --force-rebuild           rebuild Juicer core files even if existing
  --help                    show this help

Environment overrides:
  BWA, SAMTOOLS, PYTHON, JAVA, JUICER_COMMON_DIR, JUICER_JAR, JAVA_XMX

Outputs:
  juicer_core/inter_q<MAPQ>.hic
  hic_inv_context.summary.tsv
  hic_inv_context_panels.png
  hic_inv_context_panels.pdf
  per_inversion_png/*.png
  per_inversion_pdf/*.pdf
  matrices/*.npz
USAGE
}

# defaults for current project
H1_FA="/g/data/xf3/zz3507/Output/20260127Genome/H1/H1_20260127.FINAL.top20.fa"
H2_FA="/g/data/xf3/zz3507/Output/20260127Genome/H2/H2_20260127.FINAL.top20.ordered.renamed.fa"
INV_TSV="/g/data/xf3/zz3507/Output/20260127Genome/syri/syri_asm10/H1_vs_H2syri.highconfINV10kb_coordinates.tsv"
R1="/g/data/xf3/zz3507/RawData/ctrapeziformis_hic/trimmed_reads/40769_R1_001_val_1.fq.gz"
R2="/g/data/xf3/zz3507/RawData/ctrapeziformis_hic/trimmed_reads/40769_R2_001_val_2.fq.gz"
OUTDIR="/g/data/xf3/zz3507/Output/20260127Genome/hic_inv_context_joint"

THREADS=32
MAPQ=0
BIN_SIZE=100000
FLANK_BP=3000000
MIN_INV_LEN=500000
MAX_INV=8
PRIORITY_SCAFFOLDS="scaffold_5,scaffold_10,scaffold_19,scaffold_20"
INV_IDS=""
NORM="KR"
FALLBACK_NORM="VC_SQRT"
SKIP_JUICER_CORE=0
FORCE_REBUILD=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --h1-fa) H1_FA="$2"; shift 2 ;;
    --h2-fa) H2_FA="$2"; shift 2 ;;
    --inv-tsv) INV_TSV="$2"; shift 2 ;;
    --r1) R1="$2"; shift 2 ;;
    --r2) R2="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --mapq) MAPQ="$2"; shift 2 ;;
    --bin-size) BIN_SIZE="$2"; shift 2 ;;
    --flank-bp) FLANK_BP="$2"; shift 2 ;;
    --min-inv-len) MIN_INV_LEN="$2"; shift 2 ;;
    --max-inversions) MAX_INV="$2"; shift 2 ;;
    --priority-scaffolds) PRIORITY_SCAFFOLDS="$2"; shift 2 ;;
    --inv-ids) INV_IDS="$2"; shift 2 ;;
    --norm) NORM="$2"; shift 2 ;;
    --fallback-norm) FALLBACK_NORM="$2"; shift 2 ;;
    --skip-juicer-core) SKIP_JUICER_CORE=1; shift ;;
    --force-rebuild) FORCE_REBUILD=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

BWA="${BWA:-/apps/bwa/0.7.17/bin/bwa}"
SAMTOOLS="${SAMTOOLS:-/g/data/xf3/miniconda/envs/common-tools/bin/samtools}"
PYTHON="${PYTHON:-/g/data/xf3/miniconda/envs/jcvi_env/bin/python}"
JAVA="${JAVA:-/bin/java}"
JUICER_COMMON_DIR="${JUICER_COMMON_DIR:-/g/data/xf3/zz3507/Output/Assembly/assemblyQ720KOntHifiHiCVeryImportantNotDelete/assemblyQ720KOntHifiHiCBiggerThan50KContigVeryImportantNotDelete/hap2/JuicerYanVeryImportantNotDelete/scripts/common}"
JUICER_JAR="${JUICER_JAR:-$JUICER_COMMON_DIR/juicer_tools.jar}"
JAVA_XMX="${JAVA_XMX:-64g}"

CHIMERIC_AWK="$JUICER_COMMON_DIR/chimeric_blacklist.awk"
DUPS_AWK="$JUICER_COMMON_DIR/dups.awk"

for f in "$H1_FA" "$H2_FA" "$INV_TSV" "$R1" "$R2"; do
  [[ -s "$f" ]] || { echo "ERROR: Missing input file: $f" >&2; exit 1; }
done
for b in "$BWA" "$SAMTOOLS" "$PYTHON" "$JAVA"; do
  [[ -x "$b" ]] || { echo "ERROR: Missing executable: $b" >&2; exit 1; }
done
[[ -s "$JUICER_JAR" ]] || { echo "ERROR: Missing Juicer jar: $JUICER_JAR" >&2; exit 1; }
for f in "$CHIMERIC_AWK" "$DUPS_AWK"; do
  [[ -s "$f" ]] || { echo "ERROR: Missing Juicer helper file: $f" >&2; exit 1; }
done

mkdir -p "$OUTDIR"/{log,reference,juicer_core,tmp,dumps,matrices,per_inversion_png,per_inversion_pdf}
LOG="$OUTDIR/log/run_context_juicer.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "[INFO] Start: $(date)"
echo "[INFO] H1=$H1_FA"
echo "[INFO] H2=$H2_FA"
echo "[INFO] INV=$INV_TSV"
echo "[INFO] R1 =$R1"
echo "[INFO] R2 =$R2"
echo "[INFO] OUT=$OUTDIR"
echo "[INFO] THREADS=$THREADS MAPQ=$MAPQ BIN_SIZE=$BIN_SIZE FLANK_BP=$FLANK_BP MIN_INV_LEN=$MIN_INV_LEN MAX_INV=$MAX_INV"
echo "[INFO] PRIORITY_SCAFFOLDS=$PRIORITY_SCAFFOLDS"
echo "[INFO] INV_IDS=${INV_IDS:-AUTO}"
echo "[INFO] NORM=$NORM FALLBACK_NORM=$FALLBACK_NORM"

echo "[INFO] Juicer helpers: $JUICER_COMMON_DIR"
echo "[INFO] Juicer jar: $JUICER_JAR (JAVA_XMX=$JAVA_XMX)"

JOINT_FA="$OUTDIR/reference/H1H2_joint.fa"
JOINT_FAI="$JOINT_FA.fai"
CHROM_SIZES="$OUTDIR/reference/H1H2_joint.chrom.sizes"

if [[ ! -s "$JOINT_FA" || "$FORCE_REBUILD" -eq 1 ]]; then
  echo "[INFO] Build joint FASTA: $JOINT_FA"
  awk '/^>/{h=substr($0,2); split(h,a,/[	 ]+/); print ">H1|"a[1]; next} {print}' "$H1_FA" > "$JOINT_FA"
  awk '/^>/{h=substr($0,2); split(h,a,/[	 ]+/); print ">H2|"a[1]; next} {print}' "$H2_FA" >> "$JOINT_FA"
else
  echo "[INFO] Reuse joint FASTA: $JOINT_FA"
fi

if [[ ! -s "$JOINT_FAI" || "$FORCE_REBUILD" -eq 1 ]]; then
  echo "[INFO] samtools faidx"
  "$SAMTOOLS" faidx "$JOINT_FA"
fi
awk 'BEGIN{OFS="\t"}{print $1,$2}' "$JOINT_FAI" > "$CHROM_SIZES"

if [[ "$SKIP_JUICER_CORE" -eq 0 ]]; then
  if [[ ! -s "${JOINT_FA}.bwt" || ! -s "${JOINT_FA}.sa" || "$FORCE_REBUILD" -eq 1 ]]; then
    echo "[INFO] bwa index"
    "$BWA" index "$JOINT_FA"
  fi
else
  echo "[INFO] --skip-juicer-core: skip bwa index"
fi

JUICER_CORE_DIR="$OUTDIR/juicer_core"
NORM_TXT="$JUICER_CORE_DIR/readpairs.norm.txt"
ABNORM_SAM="$JUICER_CORE_DIR/readpairs.abnormal.sam"
UNMAPPED_SAM="$JUICER_CORE_DIR/readpairs.unmapped.sam"
FRAG_TXT="$JUICER_CORE_DIR/readpairs.frag.txt"
MERGED_SORT="$JUICER_CORE_DIR/merged_sort.txt"
MERGED_NODUPS="$JUICER_CORE_DIR/merged_nodups.txt"
HIC_FILE="$JUICER_CORE_DIR/inter_q${MAPQ}.hic"

if [[ "$SKIP_JUICER_CORE" -eq 1 ]]; then
  echo "[INFO] --skip-juicer-core enabled; expecting existing $HIC_FILE"
else
  if [[ "$FORCE_REBUILD" -eq 1 ]]; then
    echo "[INFO] Force rebuild Juicer core outputs"
    rm -f "$NORM_TXT" "$ABNORM_SAM" "$UNMAPPED_SAM" "$FRAG_TXT" "$MERGED_SORT" "$MERGED_NODUPS" "$HIC_FILE" \
      "$JUICER_CORE_DIR/dups.txt" "$JUICER_CORE_DIR/optdups.txt" "$JUICER_CORE_DIR/opt_dups.txt"
  fi

  if [[ ! -s "$MERGED_NODUPS" ]]; then
    echo "[INFO] BWA-MEM + Juicer chimeric filtering"
    "$BWA" mem -SP5M -t "$THREADS" "$JOINT_FA" "$R1" "$R2" \
      | awk -v "fname1=$NORM_TXT" -v "fname2=$ABNORM_SAM" -v "fname3=$UNMAPPED_SAM" -f "$CHIMERIC_AWK"

    [[ -s "$NORM_TXT" ]] || { echo "ERROR: no normal readpairs produced: $NORM_TXT" >&2; exit 1; }

    echo "[INFO] Juicer fragment assignment (site=none)"
    awk '{printf("%s %s %s %d %s %s %s %d", $1, $2, $3, 0, $4, $5, $6, 1); for (i=7; i<=NF; i++) {printf(" %s",$i);} printf("\n");}' \
      "$NORM_TXT" > "$FRAG_TXT"

    echo "[INFO] Sort readpairs"
    LC_ALL=C sort -T "$OUTDIR/tmp" -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n "$FRAG_TXT" > "$MERGED_SORT"

    echo "[INFO] Deduplicate readpairs with Juicer dups.awk"
    awk -f "$DUPS_AWK" -v "name=${JUICER_CORE_DIR}/" "$MERGED_SORT"
    [[ -s "$MERGED_NODUPS" ]] || { echo "ERROR: dedup did not produce $MERGED_NODUPS" >&2; exit 1; }
    if [[ -s "$JUICER_CORE_DIR/optdups.txt" && ! -e "$JUICER_CORE_DIR/opt_dups.txt" ]]; then
      mv "$JUICER_CORE_DIR/optdups.txt" "$JUICER_CORE_DIR/opt_dups.txt"
    fi
  else
    echo "[INFO] Reuse existing $MERGED_NODUPS"
  fi

  if [[ ! -s "$HIC_FILE" ]]; then
    echo "[INFO] Build normalized .hic with juicer_tools pre (MAPQ >= $MAPQ)"
    "$JAVA" -Djava.awt.headless=true -Xmx"$JAVA_XMX" -jar "$JUICER_JAR" pre \
      -q "$MAPQ" -t "$OUTDIR/tmp" "$MERGED_NODUPS" "$HIC_FILE" "$CHROM_SIZES"
  else
    echo "[INFO] Reuse existing .hic: $HIC_FILE"
  fi
fi

[[ -s "$HIC_FILE" ]] || { echo "ERROR: Missing .hic file: $HIC_FILE" >&2; exit 1; }

SUMMARY_TSV="$OUTDIR/hic_inv_context.summary.tsv"
PANEL_PNG="$OUTDIR/hic_inv_context_panels.png"
PANEL_PDF="$OUTDIR/hic_inv_context_panels.pdf"

echo "[INFO] Dump normalized local matrices and plot panels..."
"$PYTHON" - "$JAVA" "$JAVA_XMX" "$JUICER_JAR" "$HIC_FILE" "$INV_TSV" "$JOINT_FAI" "$OUTDIR/dumps" "$OUTDIR/matrices" "$OUTDIR/per_inversion_png" "$OUTDIR/per_inversion_pdf" "$SUMMARY_TSV" "$PANEL_PNG" "$PANEL_PDF" "$BIN_SIZE" "$FLANK_BP" "$MIN_INV_LEN" "$MAX_INV" "$PRIORITY_SCAFFOLDS" "$INV_IDS" "$NORM" "$FALLBACK_NORM" <<'PY'
import csv
import math
import os
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if len(sys.argv) != 22:
    raise SystemExit(f"Expected 21 args, got {len(sys.argv)-1}")

java_bin = sys.argv[1]
java_xmx = sys.argv[2]
juicer_jar = sys.argv[3]
hic_file = sys.argv[4]
inv_tsv = sys.argv[5]
joint_fai = sys.argv[6]
dump_dir = sys.argv[7]
mat_dir = sys.argv[8]
png_dir = sys.argv[9]
pdf_dir = sys.argv[10]
summary_tsv = sys.argv[11]
panel_png = sys.argv[12]
panel_pdf = sys.argv[13]
bin_size = int(sys.argv[14])
flank_bp = int(sys.argv[15])
min_inv_len = int(sys.argv[16])
max_inv = int(sys.argv[17])
priority_scaffolds = [x.strip() for x in sys.argv[18].split(",") if x.strip()]
inv_ids_csv = sys.argv[19].strip()
primary_norm = sys.argv[20].strip() or "KR"
fallback_norm = sys.argv[21].strip()

explicit_ids = set(x.strip() for x in inv_ids_csv.split(",") if x.strip()) if inv_ids_csv else set()

def load_fai(path):
    d = {}
    with open(path) as f:
        for ln in f:
            if not ln.strip():
                continue
            c = ln.rstrip("\n").split("\t")
            d[c[0]] = int(c[1])
    return d

def parse_inv(path):
    out = []
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        need = {"ID", "RefChr", "RefStart", "RefEnd", "QryChr", "QryStart", "QryEnd"}
        miss = need - set(rd.fieldnames or [])
        if miss:
            raise SystemExit(f"Missing inversion columns: {sorted(miss)}")
        for r in rd:
            rs, re = int(r["RefStart"]), int(r["RefEnd"])
            qs, qe = int(r["QryStart"]), int(r["QryEnd"])
            rl = abs(re - rs) + 1
            ql = abs(qe - qs) + 1
            out.append({
                "ID": r["ID"],
                "RefChr": r["RefChr"], "RefStart": rs, "RefEnd": re, "RefLen": rl,
                "QryChr": r["QryChr"], "QryStart": qs, "QryEnd": qe, "QryLen": ql,
                "MinLen": min(rl, ql),
                "H1WindowStart": parse_opt_int(r.get("H1WindowStart")),
                "H1WindowEnd": parse_opt_int(r.get("H1WindowEnd")),
                "H2WindowStart": parse_opt_int(r.get("H2WindowStart")),
                "H2WindowEnd": parse_opt_int(r.get("H2WindowEnd")),
            })
    return out

def parse_opt_int(value):
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    try:
        return int(float(s))
    except ValueError as exc:
        raise SystemExit(f"Invalid integer value in optional window column: {value}") from exc

def run_dump(norm, region1, region2, outfile):
    cmd = [
        java_bin, "-Djava.awt.headless=true", f"-Xmx{java_xmx}", "-jar", juicer_jar,
        "dump", "observed", norm, hic_file,
        region1, region2, "BP", str(bin_size), outfile,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, (proc.stderr or "").strip(), (proc.stdout or "").strip()

def dense_from_dump(path, h1_start, h1_end, h2_start, h2_end):
    h1_bins = int(math.ceil((h1_end - h1_start + 1) / float(bin_size)))
    h2_bins = int(math.ceil((h2_end - h2_start + 1) / float(bin_size)))
    mat = np.zeros((h1_bins, h2_bins), dtype=np.float64)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return mat
    with open(path) as f:
        for ln in f:
            if not ln.strip():
                continue
            c = ln.split()
            if len(c) < 3:
                continue
            try:
                x = int(float(c[0]))
                y = int(float(c[1]))
                v = float(c[2])
            except ValueError:
                continue
            if not math.isfinite(v):
                continue
            i = int((x - h1_start) // bin_size)
            j = int((y - h2_start) // bin_size)
            if i < 0 or j < 0 or i >= h1_bins or j >= h2_bins:
                continue
            mat[i, j] += v
    return mat

def fmt_mb(bp):
    return f"{bp / 1e6:.2f}".rstrip("0").rstrip(".")

GUIDE_LINEWIDTH = 2.0
TITLE_FONTSIZE = 12
AXIS_LABEL_FONTSIZE = 11
TICK_LABEL_FONTSIZE = 9
COLORBAR_LABEL_FONTSIZE = 10
COLORBAR_TICK_FONTSIZE = 9

def plot_inter(ax, mat, title, h1_inv_start, h1_inv_end, h1_win_start, h1_win_end, h2_inv_start, h2_inv_end, h2_win_start, h2_win_end, vmax):
    # Enforce same x/y scale and center the matrix in a square plotting canvas.
    nbin_y, nbin_x = mat.shape
    n_common = max(1, max(nbin_y, nbin_x))
    x_off = 0.5 * (n_common - nbin_x)
    y_off = 0.5 * (n_common - nbin_y)

    im = ax.imshow(
        np.log1p(mat),
        origin="lower",
        cmap="Reds",
        interpolation="nearest",
        aspect="equal",
        vmin=0.0,
        vmax=vmax,
        extent=(x_off - 0.5, x_off + nbin_x - 0.5, y_off - 0.5, y_off + nbin_y - 0.5),
    )
    y1 = (h1_inv_start - h1_win_start) / float(bin_size) + y_off
    y2 = (h1_inv_end - h1_win_start) / float(bin_size) + y_off
    x1 = (h2_inv_start - h2_win_start) / float(bin_size) + x_off
    x2 = (h2_inv_end - h2_win_start) / float(bin_size) + x_off
    for x in (x1, x2):
        ax.axvline(x=x, color="deepskyblue", lw=GUIDE_LINEWIDTH, ls="--")
    for y in (y1, y2):
        ax.axhline(y=y, color="deepskyblue", lw=GUIDE_LINEWIDTH, ls="--")
    ax.set_title(title, fontsize=TITLE_FONTSIZE, pad=12)
    ax.set_xlim(-0.5, n_common - 0.5)
    ax.set_ylim(-0.5, n_common - 0.5)
    ax.set_anchor("C")

    if nbin_x > 1:
        xticks = np.linspace(x_off, x_off + nbin_x - 1, 5)
        xlabels = [fmt_mb(x) for x in np.linspace(h2_win_start, h2_win_end, 5)]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels, fontsize=TICK_LABEL_FONTSIZE)
    if nbin_y > 1:
        yticks = np.linspace(y_off, y_off + nbin_y - 1, 5)
        ylabels = [fmt_mb(y) for y in np.linspace(h1_win_start, h1_win_end, 5)]
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, fontsize=TICK_LABEL_FONTSIZE)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE)
    ax.set_xlabel("H2 Position (Mb)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("H1 Position (Mb)", fontsize=AXIS_LABEL_FONTSIZE)
    return im

fai = load_fai(joint_fai)
inv_all = parse_inv(inv_tsv)

if explicit_ids:
    selected = [r for r in inv_all if r["ID"] in explicit_ids]
    missing = sorted(explicit_ids - {r["ID"] for r in selected})
    if missing:
        raise SystemExit(f"Requested IDs not found: {','.join(missing)}")
else:
    selected = [
        r for r in inv_all
        if r["MinLen"] >= min_inv_len
        and (r["RefChr"] in priority_scaffolds or r["QryChr"] in priority_scaffolds)
    ]
    if not selected:
        selected = [r for r in inv_all if r["MinLen"] >= min_inv_len]
    selected.sort(key=lambda x: x["MinLen"], reverse=True)
    selected = selected[:max_inv]

if not selected:
    raise SystemExit("No inversion selected. Adjust --min-inv-len or --inv-ids.")

rows = []
mat_store = {}
all_logs = []

for r in selected:
    inv_id = r["ID"]
    h1_chr = f"H1|{r['RefChr']}"
    h2_chr = f"H2|{r['QryChr']}"
    if h1_chr not in fai:
        raise SystemExit(f"{h1_chr} not found in joint FAI.")
    if h2_chr not in fai:
        raise SystemExit(f"{h2_chr} not found in joint FAI.")

    custom_h1 = [r["H1WindowStart"], r["H1WindowEnd"]]
    custom_h2 = [r["H2WindowStart"], r["H2WindowEnd"]]
    if any(v is not None for v in custom_h1) and not all(v is not None for v in custom_h1):
        raise SystemExit(f"{inv_id}: both H1WindowStart and H1WindowEnd are required when overriding the H1 window.")
    if any(v is not None for v in custom_h2) and not all(v is not None for v in custom_h2):
        raise SystemExit(f"{inv_id}: both H2WindowStart and H2WindowEnd are required when overriding the H2 window.")

    if all(v is not None for v in custom_h1):
        h1_start = max(1, min(custom_h1))
        h1_end = min(fai[h1_chr], max(custom_h1))
    else:
        h1_start = max(1, min(r["RefStart"], r["RefEnd"]) - flank_bp)
        h1_end = min(fai[h1_chr], max(r["RefStart"], r["RefEnd"]) + flank_bp)

    if all(v is not None for v in custom_h2):
        h2_start = max(1, min(custom_h2))
        h2_end = min(fai[h2_chr], max(custom_h2))
    else:
        h2_start = max(1, min(r["QryStart"], r["QryEnd"]) - flank_bp)
        h2_end = min(fai[h2_chr], max(r["QryStart"], r["QryEnd"]) + flank_bp)

    if not (h1_start <= min(r["RefStart"], r["RefEnd"]) <= h1_end and h1_start <= max(r["RefStart"], r["RefEnd"]) <= h1_end):
        raise SystemExit(f"{inv_id}: H1 window does not fully contain the inversion coordinates.")
    if not (h2_start <= min(r["QryStart"], r["QryEnd"]) <= h2_end and h2_start <= max(r["QryStart"], r["QryEnd"]) <= h2_end):
        raise SystemExit(f"{inv_id}: H2 window does not fully contain the inversion coordinates.")

    reg1 = f"{h1_chr}:{h1_start}:{h1_end}"
    reg2 = f"{h2_chr}:{h2_start}:{h2_end}"

    dump_tsv = os.path.join(dump_dir, f"{inv_id}.H1vsH2.{r['RefChr']}_{r['QryChr']}.{h1_start}_{h1_end}.{h2_start}_{h2_end}.{primary_norm}.tsv")
    rc, err, _ = run_dump(primary_norm, reg1, reg2, dump_tsv)
    used_norm = primary_norm

    if rc != 0 and fallback_norm and fallback_norm != primary_norm:
        dump_tsv = os.path.join(dump_dir, f"{inv_id}.H1vsH2.{r['RefChr']}_{r['QryChr']}.{h1_start}_{h1_end}.{h2_start}_{h2_end}.{fallback_norm}.tsv")
        rc2, err2, _ = run_dump(fallback_norm, reg1, reg2, dump_tsv)
        if rc2 == 0:
            rc = 0
            used_norm = fallback_norm
            err = ""
        else:
            err = (err + "\n" + err2).strip()

    if rc != 0:
        raise SystemExit(
            f"juicer_tools dump failed for {inv_id} H1-vs-H2\n"
            f"{reg1} x {reg2}\n"
            f"Primary={primary_norm} fallback={fallback_norm}\n{err}"
        )

    mat = dense_from_dump(dump_tsv, h1_start, h1_end, h2_start, h2_end)
    mat_store[inv_id] = {
        "h1_chr": h1_chr,
        "h2_chr": h2_chr,
        "h1_start": h1_start,
        "h1_end": h1_end,
        "h2_start": h2_start,
        "h2_end": h2_end,
        "h1_inv_start": min(r["RefStart"], r["RefEnd"]),
        "h1_inv_end": max(r["RefStart"], r["RefEnd"]),
        "h2_inv_start": min(r["QryStart"], r["QryEnd"]),
        "h2_inv_end": max(r["QryStart"], r["QryEnd"]),
        "norm": used_norm,
        "dump_tsv": dump_tsv,
        "mat": mat,
    }

    finite_pos = np.isfinite(mat) & (mat > 0)
    nz = np.log1p(mat[finite_pos])
    if nz.size:
        all_logs.append(nz)

# global scale across all panels for standardized comparison
if all_logs:
    all_vals = np.concatenate(all_logs)
    global_vmax = float(np.quantile(all_vals, 0.99))
else:
    global_vmax = 1.0

global_vmax = max(global_vmax, 1e-6)

for r in selected:
    inv_id = r["ID"]
    m = mat_store[inv_id]
    rows.append({
        "ID": inv_id,
        "RefChr": r["RefChr"], "RefStart": r["RefStart"], "RefEnd": r["RefEnd"], "RefLen": r["RefLen"],
        "QryChr": r["QryChr"], "QryStart": r["QryStart"], "QryEnd": r["QryEnd"], "QryLen": r["QryLen"],
        "MinLen": r["MinLen"],
        "H1_window_start": m["h1_start"], "H1_window_end": m["h1_end"],
        "H2_window_start": m["h2_start"], "H2_window_end": m["h2_end"],
        "Inter_contacts": int(np.nansum(m["mat"])),
        "Norm": m["norm"],
        "Inter_dump": m["dump_tsv"],
        "hic_file": hic_file,
        "bin_size": bin_size,
        "flank_bp": flank_bp,
    })

rows.sort(key=lambda x: x["MinLen"], reverse=True)

header = [
    "ID", "RefChr", "RefStart", "RefEnd", "RefLen",
    "QryChr", "QryStart", "QryEnd", "QryLen", "MinLen",
    "H1_window_start", "H1_window_end",
    "H2_window_start", "H2_window_end",
    "Inter_contacts", "Norm", "Inter_dump",
    "hic_file", "bin_size", "flank_bp",
]
with open(summary_tsv, "w", newline="") as f:
    w = csv.DictWriter(f, delimiter="\t", fieldnames=header)
    w.writeheader()
    for row in rows:
        w.writerow(row)

# Save matrices
for r in rows:
    inv_id = r["ID"]
    m = mat_store[inv_id]
    np.savez_compressed(
        os.path.join(mat_dir, f"{inv_id}.hic_context.juicer.npz"),
        H1_vs_H2=m["mat"],
        H1_chr=m["h1_chr"], H1_start=m["h1_start"], H1_end=m["h1_end"],
        H2_chr=m["h2_chr"], H2_start=m["h2_start"], H2_end=m["h2_end"],
        H1_inv_start=m["h1_inv_start"], H1_inv_end=m["h1_inv_end"],
        H2_inv_start=m["h2_inv_start"], H2_inv_end=m["h2_inv_end"],
        norm=m["norm"],
        bin_size=bin_size,
    )

# Combined panel (cap panel count to keep figure tractable for large candidate sets)
n_total = len(rows)
combined_max = 20
panel_rows = rows[:combined_max]
n = len(panel_rows)

fig_h = max(7.0 * n + 2.0, 12.0)
fig, axes = plt.subplots(nrows=n, ncols=1, figsize=(19.0, fig_h), constrained_layout=False)
if n == 1:
    axes = np.array([axes])

for i, r in enumerate(panel_rows):
    inv_id = r["ID"]
    m = mat_store[inv_id]
    title = (
        f"{inv_id}  H1 {r['RefChr']}:{r['RefStart']}-{r['RefEnd']}  vs  "
        f"H2 {r['QryChr']}:{r['QryStart']}-{r['QryEnd']} ({m['norm']})"
    )
    im = plot_inter(
        axes[i],
        m["mat"],
        title,
        m["h1_inv_start"],
        m["h1_inv_end"],
        m["h1_start"],
        m["h1_end"],
        m["h2_inv_start"],
        m["h2_inv_end"],
        m["h2_start"],
        m["h2_end"],
        global_vmax,
    )
    cb = fig.colorbar(im, ax=axes[i], fraction=0.025, pad=0.02)
    cb.set_label("log1p(normalized contact)", fontsize=COLORBAR_LABEL_FONTSIZE)
    cb.ax.tick_params(labelsize=COLORBAR_TICK_FONTSIZE)

if n_total > n:
    fig.suptitle(
        f"Inter-haplotype Hi-C context maps (first {n} of {n_total}; H1 window vs H2 window)",
        fontsize=15,
        y=0.99,
    )
else:
    fig.suptitle("Inter-haplotype Hi-C context maps (H1 window vs H2 window)", fontsize=15, y=0.99)
fig.subplots_adjust(left=0.07, right=0.90, top=0.93, bottom=0.05, hspace=0.72)
fig.savefig(panel_png, dpi=400, bbox_inches="tight", pad_inches=0.25)
fig.savefig(panel_pdf, bbox_inches="tight", pad_inches=0.25)
plt.close(fig)

# Per inversion panels
for r in rows:
    inv_id = r["ID"]
    m = mat_store[inv_id]
    fig, ax = plt.subplots(1, 1, figsize=(13.5, 10.0), constrained_layout=False)
    im = plot_inter(
        ax,
        m["mat"],
        f"{inv_id} H1 {r['RefChr']} vs H2 {r['QryChr']} ({m['norm']})",
        m["h1_inv_start"],
        m["h1_inv_end"],
        m["h1_start"],
        m["h1_end"],
        m["h2_inv_start"],
        m["h2_inv_end"],
        m["h2_start"],
        m["h2_end"],
        global_vmax,
    )
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("log1p(normalized contact)", fontsize=COLORBAR_LABEL_FONTSIZE)
    cb.ax.tick_params(labelsize=COLORBAR_TICK_FONTSIZE)
    fig.subplots_adjust(left=0.09, right=0.89, top=0.90, bottom=0.10)
    fig.savefig(os.path.join(png_dir, f"{inv_id}.hic_context.juicer.png"), dpi=400, bbox_inches="tight", pad_inches=0.25)
    fig.savefig(os.path.join(pdf_dir, f"{inv_id}.hic_context.juicer.pdf"), bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)

print(f"[INFO] wrote summary: {summary_tsv}", file=sys.stderr)
print(f"[INFO] wrote panel png: {panel_png}", file=sys.stderr)
print(f"[INFO] wrote panel pdf: {panel_pdf}", file=sys.stderr)
PY

echo "[INFO] Done: $(date)"
echo "[INFO] HIC file : $HIC_FILE"
echo "[INFO] Summary  : $SUMMARY_TSV"
echo "[INFO] Panels   : $PANEL_PNG"
echo "[INFO] Panels   : $PANEL_PDF"
echo "[INFO] Per-INV  : $OUTDIR/per_inversion_png"
echo "[INFO] Matrices : $OUTDIR/matrices"
