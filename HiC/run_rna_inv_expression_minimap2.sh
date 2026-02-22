#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Map RNA reads with minimap2 and summarize expression for selected inversions.

Pipeline:
1) Build minimap2 index for reference genome (if missing)
2) Align paired-end RNA reads with minimap2 and produce sorted/indexed BAMs
3) Quantify gene counts with featureCounts
4) Map gene expression to inversion-overlapping genes

Usage:
  bash run_rna_inv_expression_minimap2.sh [options]

Options:
  --reads-dir PATH        Directory with paired FASTQ files (default dataset path)
  --r1-suffix STR         Read1 suffix (default: _1.fastq.gz)
  --r2-suffix STR         Read2 suffix (default: _2.fastq.gz)
  --sample-list PATH      Optional file with one sample prefix per line
  --ref-fa PATH           Reference FASTA (default: 2026 H1 top20)
  --annotation-gtf PATH   GTF for featureCounts (default: 2026 H1 braker.gtf)
  --counts-tsv PATH       Gene count matrix path (default: outdir/rna_gene_counts.tsv)
  --inv-overlaps PATH     Inversion-gene overlap table (default: H1_INV_gene_overlaps.tsv)
  --inv-ids STR|PATH      IDs list or file (default: 1855,1963,1915,1926,1903,1820,1936)
  --outdir PATH           Output directory
  --threads INT           Threads (default: 32)
  --mm2-preset STR        minimap2 preset for -ax (default: splice)
  --minimap2 PATH         minimap2 binary (default: /g/data/xf3/miniconda/bin/minimap2)
  --samtools PATH         samtools binary (default: /g/data/xf3/miniconda/envs/common-tools/bin/samtools)
  --featurecounts PATH    featureCounts binary (default: featureCounts)
  --python PATH           Python interpreter (default: /g/data/xf3/miniconda/bin/python3)
  --skip-align            Skip minimap2 alignment; reuse BAMs in outdir/bam
  --skip-counts           Skip featureCounts; reuse --counts-tsv
  --skip-map              Skip inversion mapping stage
  --help                  Show this help message

Notes:
- minimap2 with RNA short reads is configured as: -ax splice --secondary=no
- For short-read RNA-seq, HISAT2/STAR is typically preferred, but this script
  is minimap2-based as requested.
EOF
}

# Defaults for current dataset
READS_DIR="/g/data/xf3/zz3507/RawData/ctrapeziformis_RNA/chiloglottis"
R1_SUFFIX="_1.fastq.gz"
R2_SUFFIX="_2.fastq.gz"
SAMPLE_LIST=""
REF_FA="/g/data/xf3/zz3507/Output/20260127Genome/H1/H1_20260127.FINAL.top20.fa"
ANNOT_GTF="/g/data/xf3/zz3507/Output/20260127Genome/H1/breaker/braker.gtf"
COUNTS_TSV=""
INV_OVERLAPS="/g/data/xf3/zz3507/Output/20260127Genome/H1/breaker/H1_INV_gene_overlaps.tsv"
INV_IDS="1855,1963,1915,1926,1903,1820,1936"
OUTDIR="/g/data/xf3/zz3507/Output/20260127Genome/H1/inversion_rna_minimap2"
THREADS=32
MM2_PRESET="splice"
MINIMAP2="${MINIMAP2:-/g/data/xf3/miniconda/bin/minimap2}"
SAMTOOLS="${SAMTOOLS:-/g/data/xf3/miniconda/envs/common-tools/bin/samtools}"
FEATURECOUNTS="${FEATURECOUNTS:-featureCounts}"
PYTHON="${PYTHON:-/g/data/xf3/miniconda/bin/python3}"

SKIP_ALIGN=0
SKIP_COUNTS=0
SKIP_MAP=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reads-dir) READS_DIR="$2"; shift 2 ;;
    --r1-suffix) R1_SUFFIX="$2"; shift 2 ;;
    --r2-suffix) R2_SUFFIX="$2"; shift 2 ;;
    --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
    --ref-fa) REF_FA="$2"; shift 2 ;;
    --annotation-gtf) ANNOT_GTF="$2"; shift 2 ;;
    --counts-tsv) COUNTS_TSV="$2"; shift 2 ;;
    --inv-overlaps) INV_OVERLAPS="$2"; shift 2 ;;
    --inv-ids) INV_IDS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --mm2-preset) MM2_PRESET="$2"; shift 2 ;;
    --minimap2) MINIMAP2="$2"; shift 2 ;;
    --samtools) SAMTOOLS="$2"; shift 2 ;;
    --featurecounts) FEATURECOUNTS="$2"; shift 2 ;;
    --python) PYTHON="$2"; shift 2 ;;
    --skip-align) SKIP_ALIGN=1; shift ;;
    --skip-counts) SKIP_COUNTS=1; shift ;;
    --skip-map) SKIP_MAP=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

MAP_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/map_rna_expression_to_inversions.py"

for f in "$REF_FA" "$INV_OVERLAPS" "$MAP_SCRIPT"; do
  [[ -s "$f" ]] || { echo "ERROR: Missing input file: $f" >&2; exit 1; }
done

for b in "$MINIMAP2" "$SAMTOOLS" "$PYTHON"; do
  [[ -x "$b" ]] || { echo "ERROR: Missing executable: $b" >&2; exit 1; }
done

mkdir -p "$OUTDIR"/{log,bam,tmp}
LOG="$OUTDIR/log/run.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "[INFO] Start: $(date)"
echo "[INFO] READS_DIR=$READS_DIR"
echo "[INFO] REF_FA=$REF_FA"
echo "[INFO] ANNOT_GTF=$ANNOT_GTF"
echo "[INFO] INV_OVERLAPS=$INV_OVERLAPS"
echo "[INFO] INV_IDS=$INV_IDS"
echo "[INFO] OUTDIR=$OUTDIR"
echo "[INFO] THREADS=$THREADS MM2_PRESET=$MM2_PRESET"
echo "[INFO] SKIP_ALIGN=$SKIP_ALIGN SKIP_COUNTS=$SKIP_COUNTS SKIP_MAP=$SKIP_MAP"

REF_MMI="$OUTDIR/reference.mmi"
BAM_DIR="$OUTDIR/bam"
if [[ -z "$COUNTS_TSV" ]]; then
  COUNTS_TSV="$OUTDIR/rna_gene_counts.tsv"
fi
OUT_PREFIX="$OUTDIR/inv_rna_expr_minimap2"

if [[ $SKIP_ALIGN -eq 0 ]]; then
  if [[ ! -s "$REF_MMI" ]]; then
    echo "[INFO] Building minimap2 index: $REF_MMI"
    "$MINIMAP2" -d "$REF_MMI" "$REF_FA"
  else
    echo "[INFO] Reuse minimap2 index: $REF_MMI"
  fi

  declare -a SAMPLES
  if [[ -n "$SAMPLE_LIST" ]]; then
    [[ -s "$SAMPLE_LIST" ]] || { echo "ERROR: sample list not found: $SAMPLE_LIST" >&2; exit 1; }
    mapfile -t SAMPLES < <(grep -v '^[[:space:]]*#' "$SAMPLE_LIST" | awk 'NF{print $1}')
  else
    mapfile -t SAMPLES < <(
      find "$READS_DIR" -maxdepth 1 -type f -name "*${R1_SUFFIX}" -printf '%f\n' \
      | sed "s/${R1_SUFFIX//\//\\/}$//" \
      | sort -u
    )
  fi
  [[ ${#SAMPLES[@]} -gt 0 ]] || { echo "ERROR: No samples found." >&2; exit 1; }
  echo "[INFO] Samples to process: ${#SAMPLES[@]}"

  for s in "${SAMPLES[@]}"; do
    r1="$READS_DIR/${s}${R1_SUFFIX}"
    r2="$READS_DIR/${s}${R2_SUFFIX}"
    bam="$BAM_DIR/${s}.sorted.bam"
    [[ -s "$r1" ]] || { echo "[WARN] Missing R1 for $s (skip): $r1"; continue; }
    [[ -s "$r2" ]] || { echo "[WARN] Missing R2 for $s (skip): $r2"; continue; }

    if [[ -s "$bam" && -s "${bam}.bai" ]]; then
      echo "[INFO] Reuse BAM: $bam"
      continue
    fi

    echo "[INFO] Align: $s"
    "$MINIMAP2" -t "$THREADS" -ax "$MM2_PRESET" --secondary=no "$REF_MMI" "$r1" "$r2" \
      | "$SAMTOOLS" sort -@ "$THREADS" -o "$bam" -
    "$SAMTOOLS" index "$bam"
  done
fi

mapfile -t BAMS < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*.sorted.bam" | sort)
[[ ${#BAMS[@]} -gt 0 ]] || { echo "ERROR: No BAM files found in $BAM_DIR" >&2; exit 1; }
echo "[INFO] BAM files: ${#BAMS[@]}"

if [[ $SKIP_COUNTS -eq 0 ]]; then
  [[ -s "$ANNOT_GTF" ]] || { echo "ERROR: Missing annotation GTF: $ANNOT_GTF" >&2; exit 1; }
  if [[ "$FEATURECOUNTS" == */* ]]; then
    [[ -x "$FEATURECOUNTS" ]] || { echo "ERROR: Missing featureCounts executable: $FEATURECOUNTS" >&2; exit 1; }
  else
    command -v "$FEATURECOUNTS" >/dev/null 2>&1 || {
      echo "ERROR: featureCounts not found in PATH; use --featurecounts /path/to/featureCounts" >&2
      exit 1
    }
  fi
  # Auto-detect whether BAMs retain paired-end FLAGs; minimap2 splice output in this
  # workflow may not mark records as paired, which breaks `featureCounts -p`.
  pe_flag_count="$("$SAMTOOLS" view -c -f 1 "${BAMS[0]}" | tr -d '[:space:]')"
  if [[ "${pe_flag_count:-0}" =~ ^[0-9]+$ ]] && [[ "$pe_flag_count" -gt 0 ]]; then
    echo "[INFO] Running featureCounts in paired-end mode (-p)"
    "$FEATURECOUNTS" -T "$THREADS" -a "$ANNOT_GTF" -o "$COUNTS_TSV" -p "${BAMS[@]}"
  else
    echo "[WARN] No paired-end FLAGs detected in BAMs; running featureCounts without -p (single-end counting mode)."
    "$FEATURECOUNTS" -T "$THREADS" -a "$ANNOT_GTF" -o "$COUNTS_TSV" "${BAMS[@]}"
  fi
else
  [[ -s "$COUNTS_TSV" ]] || { echo "ERROR: --skip-counts set but missing $COUNTS_TSV" >&2; exit 1; }
  echo "[INFO] Reuse counts: $COUNTS_TSV"
fi

if [[ $SKIP_MAP -eq 0 ]]; then
  echo "[INFO] Mapping expression to inversions"
  "$PYTHON" "$MAP_SCRIPT" \
    --inv-overlaps "$INV_OVERLAPS" \
    --expr "$COUNTS_TSV" \
    --inv-ids "$INV_IDS" \
    --out-prefix "$OUT_PREFIX"
else
  echo "[INFO] Skip inversion mapping stage"
fi

echo "[INFO] Done: $(date)"
echo "[INFO] Outputs:"
echo "  BAM dir: $BAM_DIR"
echo "  Counts : $COUNTS_TSV"
echo "  Prefix : $OUT_PREFIX"
