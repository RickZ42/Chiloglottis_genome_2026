#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Create publication-grade Hi-C inversion context maps using Juicer-compatible processing.

This script:
1) Builds a joint H1+H2 reference with prefixed chromosome IDs.
2) Maps Hi-C reads with BWA-MEM (`-SP5M`) as in Juicer single-CPU workflow.
3) Runs Juicer core text processing (`chimeric_blacklist.awk`, fragment assignment, `dups.awk`).
4) Builds normalized `.hic` with `juicer_tools pre`.
5) Extracts local windows with `juicer_tools dump observed` and plots standardized panels.

Usage:
  bash run_joint_hic_inversion_context_map.sh [options]

Options:
  --h1-fa PATH              H1 FASTA
  --h2-fa PATH              H2 FASTA
  --inv-tsv PATH            inversion TSV (columns: ID,RefChr,RefStart,RefEnd,QryChr,QryStart,QryEnd)
  --r1 PATH                 Hi-C read1 FASTQ(.gz)
  --r2 PATH                 Hi-C read2 FASTQ(.gz)
  --outdir PATH             output directory
  --threads INT             threads for BWA/sort (default: 32)
  --mapq INT                MAPQ cutoff used by juicer_tools pre (default: 30)
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
MAPQ=30
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

if [[ ! -s "${JOINT_FA}.bwt" || ! -s "${JOINT_FA}.sa" || "$FORCE_REBUILD" -eq 1 ]]; then
  echo "[INFO] bwa index"
  "$BWA" index "$JOINT_FA"
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
            })
    return out

def run_dump(norm, chrom, start, end, outfile):
    region = f"{chrom}:{start}:{end}"
    cmd = [
        java_bin, "-Djava.awt.headless=true", f"-Xmx{java_xmx}", "-jar", juicer_jar,
        "dump", "observed", norm, hic_file,
        region, region, "BP", str(bin_size), outfile,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, (proc.stderr or "").strip(), (proc.stdout or "").strip()

def dense_from_dump(path, start, end):
    nbin = int(math.ceil((end - start + 1) / float(bin_size)))
    mat = np.zeros((nbin, nbin), dtype=np.float64)
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
            i = int((x - start) // bin_size)
            j = int((y - start) // bin_size)
            if i < 0 or j < 0 or i >= nbin or j >= nbin:
                continue
            mat[i, j] += v
            if i != j:
                mat[j, i] += v
    return mat

def plot_one(ax, mat, title, inv_start, inv_end, win_start, vmax):
    im = ax.imshow(np.log1p(mat), origin="lower", cmap="Reds", interpolation="nearest", vmin=0.0, vmax=vmax)
    sbin = (inv_start - win_start) / float(bin_size)
    ebin = (inv_end - win_start) / float(bin_size)
    for x in (sbin, ebin):
        ax.axvline(x=x, color="deepskyblue", lw=1.1, ls="--")
        ax.axhline(y=x, color="deepskyblue", lw=1.1, ls="--")
    ax.set_title(title, fontsize=9)
    nbin = mat.shape[0]
    if nbin > 1:
        ticks = np.linspace(0, nbin - 1, 5)
        labels = [f"{(win_start + t * bin_size) / 1e6:.1f}" for t in ticks]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Position (Mb)")
    ax.set_ylabel("Position (Mb)")
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
    for hap, chr_raw, s, e in (
        ("H1", r["RefChr"], r["RefStart"], r["RefEnd"]),
        ("H2", r["QryChr"], r["QryStart"], r["QryEnd"]),
    ):
        chrn = f"{hap}|{chr_raw}"
        if chrn not in fai:
            raise SystemExit(f"{chrn} not found in joint FAI.")
        clen = fai[chrn]
        a = max(1, min(s, e) - flank_bp)
        b = min(clen, max(s, e) + flank_bp)

        dump_tsv = os.path.join(dump_dir, f"{inv_id}.{hap}.{chr_raw}.{a}_{b}.{primary_norm}.tsv")
        rc, err, _ = run_dump(primary_norm, chrn, a, b, dump_tsv)
        used_norm = primary_norm

        if rc != 0 and fallback_norm and fallback_norm != primary_norm:
            dump_tsv = os.path.join(dump_dir, f"{inv_id}.{hap}.{chr_raw}.{a}_{b}.{fallback_norm}.tsv")
            rc2, err2, _ = run_dump(fallback_norm, chrn, a, b, dump_tsv)
            if rc2 == 0:
                rc = 0
                used_norm = fallback_norm
                err = ""
            else:
                err = (err + "\n" + err2).strip()

        if rc != 0:
            raise SystemExit(
                f"juicer_tools dump failed for {inv_id} {hap} {chrn}:{a}-{b}\n"
                f"Primary={primary_norm} fallback={fallback_norm}\n{err}"
            )

        mat = dense_from_dump(dump_tsv, a, b)
        mat_store[(inv_id, hap)] = {
            "chr": chrn,
            "start": a,
            "end": b,
            "inv_start": min(s, e),
            "inv_end": max(s, e),
            "norm": used_norm,
            "dump_tsv": dump_tsv,
            "mat": mat,
        }

        nz = np.log1p(mat[mat > 0])
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
    h1 = mat_store[(inv_id, "H1")]
    h2 = mat_store[(inv_id, "H2")]
    rows.append({
        "ID": inv_id,
        "RefChr": r["RefChr"], "RefStart": r["RefStart"], "RefEnd": r["RefEnd"], "RefLen": r["RefLen"],
        "QryChr": r["QryChr"], "QryStart": r["QryStart"], "QryEnd": r["QryEnd"], "QryLen": r["QryLen"],
        "MinLen": r["MinLen"],
        "H1_window_start": h1["start"], "H1_window_end": h1["end"], "H1_contacts": int(h1["mat"].sum() // 2), "H1_norm": h1["norm"], "H1_dump": h1["dump_tsv"],
        "H2_window_start": h2["start"], "H2_window_end": h2["end"], "H2_contacts": int(h2["mat"].sum() // 2), "H2_norm": h2["norm"], "H2_dump": h2["dump_tsv"],
        "hic_file": hic_file,
        "bin_size": bin_size,
        "flank_bp": flank_bp,
    })

rows.sort(key=lambda x: x["MinLen"], reverse=True)

header = [
    "ID", "RefChr", "RefStart", "RefEnd", "RefLen",
    "QryChr", "QryStart", "QryEnd", "QryLen", "MinLen",
    "H1_window_start", "H1_window_end", "H1_contacts", "H1_norm", "H1_dump",
    "H2_window_start", "H2_window_end", "H2_contacts", "H2_norm", "H2_dump",
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
    h1 = mat_store[(inv_id, "H1")]
    h2 = mat_store[(inv_id, "H2")]
    np.savez_compressed(
        os.path.join(mat_dir, f"{inv_id}.hic_context.juicer.npz"),
        H1=h1["mat"], H2=h2["mat"],
        H1_chr=h1["chr"], H1_start=h1["start"], H1_end=h1["end"], H1_norm=h1["norm"],
        H2_chr=h2["chr"], H2_start=h2["start"], H2_end=h2["end"], H2_norm=h2["norm"],
        bin_size=bin_size,
    )

# Combined panel
n = len(rows)
fig_h = max(3.2 * n, 4.0)
fig, axes = plt.subplots(nrows=n, ncols=2, figsize=(11, fig_h), constrained_layout=True)
if n == 1:
    axes = np.array([axes])

for i, r in enumerate(rows):
    inv_id = r["ID"]
    h1 = mat_store[(inv_id, "H1")]
    h2 = mat_store[(inv_id, "H2")]
    title1 = f"{inv_id} H1 {r['RefChr']}:{r['RefStart']}-{r['RefEnd']} ({h1['norm']})"
    title2 = f"{inv_id} H2 {r['QryChr']}:{r['QryStart']}-{r['QryEnd']} ({h2['norm']})"
    im = plot_one(axes[i, 0], h1["mat"], title1, r["RefStart"], r["RefEnd"], h1["start"], global_vmax)
    plot_one(axes[i, 1], h2["mat"], title2, r["QryStart"], r["QryEnd"], h2["start"], global_vmax)
    cb = fig.colorbar(im, ax=axes[i, :], fraction=0.02, pad=0.01)
    cb.set_label("log1p(normalized contact)")

fig.suptitle("Hi-C inversion context maps (Juicer standardized normalization)", fontsize=12)
fig.savefig(panel_png, dpi=400)
fig.savefig(panel_pdf)
plt.close(fig)

# Per inversion panels
for r in rows:
    inv_id = r["ID"]
    h1 = mat_store[(inv_id, "H1")]
    h2 = mat_store[(inv_id, "H2")]
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    im = plot_one(ax[0], h1["mat"], f"{inv_id} H1 {r['RefChr']} ({h1['norm']})", r["RefStart"], r["RefEnd"], h1["start"], global_vmax)
    plot_one(ax[1], h2["mat"], f"{inv_id} H2 {r['QryChr']} ({h2['norm']})", r["QryStart"], r["QryEnd"], h2["start"], global_vmax)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("log1p(normalized contact)")
    fig.savefig(os.path.join(png_dir, f"{inv_id}.hic_context.juicer.png"), dpi=400)
    fig.savefig(os.path.join(pdf_dir, f"{inv_id}.hic_context.juicer.pdf"))
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
