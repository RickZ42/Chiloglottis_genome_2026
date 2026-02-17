#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Re-run Hi-C inversion validation using joint mapping to both haplotypes.

This script:
1) Builds a joint reference (H1+H2 with prefixed scaffold names).
2) Maps Hi-C PE reads to the joint reference.
3) Evaluates SyRI inversion breakpoints with breakpoint-crossing Hi-C pairs.
4) Writes a per-inversion "independent support" conclusion.

Usage:
  bash run_joint_hic_inv_validation.sh [options]

Options:
  --h1-fa PATH              H1 FASTA
  --h2-fa PATH              H2 FASTA
  --inv-tsv PATH            SyRI inversion table (with Ref/Qry coords)
  --r1 PATH                 Hi-C read1 FASTQ(.gz)
  --r2 PATH                 Hi-C read2 FASTQ(.gz)
  --outdir PATH             Output directory
  --threads INT             Threads for minimap2/samtools (default: 32)
  --mapq INT                Mapping quality cutoff (default: 20)
  --flank-bp INT            Breakpoint flank window size in bp (default: 200000)
  --min-bp-pairs INT        Minimum pairs on each breakpoint side (default: 25)
  --min-bp-to-ctrl-frac F   Min breakpoint/control ratio (default: 0.15)
  --priority-scaffolds CSV  Priority scaffolds (default: scaffold_5,scaffold_10,scaffold_19,scaffold_20)
  --min-inv-len INT         Minimum min(ref_len,qry_len) to analyze (default: 0)
  --help                    Show this message

Environment overrides:
  MINIMAP2, SAMTOOLS, PYTHON

Outputs:
  inv_joint_hic_support.summary.tsv
  inv_joint_hic_support.priority.tsv
  inv_joint_hic_support.conclusion.txt
  H1H2_joint.fa(.fai/.mmi)
EOF
}

# Defaults for current 2026 Chiloglottis dataset
H1_FA="/g/data/xf3/zz3507/Output/20260127Genome/H1/H1_20260127.FINAL.top20.fa"
H2_FA="/g/data/xf3/zz3507/Output/20260127Genome/H2/H2_20260127.FINAL.top20.ordered.renamed.fa"
INV_TSV="/g/data/xf3/zz3507/Output/20260127Genome/syri/syri_asm10/H1_vs_H2syri.highconfINV10kb_coordinates.tsv"
R1="/g/data/xf3/zz3507/RawData/ctrapeziformis_hic/trimmed_reads/40769_R1_001_val_1.fq.gz"
R2="/g/data/xf3/zz3507/RawData/ctrapeziformis_hic/trimmed_reads/40769_R2_001_val_2.fq.gz"
OUTDIR="/g/data/xf3/zz3507/Output/20260127Genome/hic_inv_validation_joint"
THREADS=32
MAPQ=20
FLANK_BP=200000
MIN_BP_PAIRS=25
MIN_BP_TO_CTRL_FRAC=0.15
PRIORITY_SCAFFOLDS="scaffold_5,scaffold_10,scaffold_19,scaffold_20"
MIN_INV_LEN=0

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
    --flank-bp) FLANK_BP="$2"; shift 2 ;;
    --min-bp-pairs) MIN_BP_PAIRS="$2"; shift 2 ;;
    --min-bp-to-ctrl-frac) MIN_BP_TO_CTRL_FRAC="$2"; shift 2 ;;
    --priority-scaffolds) PRIORITY_SCAFFOLDS="$2"; shift 2 ;;
    --min-inv-len) MIN_INV_LEN="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

MINIMAP2="${MINIMAP2:-/g/data/xf3/miniconda/bin/minimap2}"
SAMTOOLS="${SAMTOOLS:-/g/data/xf3/miniconda/envs/common-tools/bin/samtools}"
PYTHON="${PYTHON:-/g/data/xf3/miniconda/bin/python}"

for f in "$H1_FA" "$H2_FA" "$INV_TSV" "$R1" "$R2"; do
  [[ -s "$f" ]] || { echo "ERROR: Missing input file: $f" >&2; exit 1; }
done
for b in "$MINIMAP2" "$SAMTOOLS" "$PYTHON"; do
  [[ -x "$b" ]] || { echo "ERROR: Missing executable: $b" >&2; exit 1; }
done

mkdir -p "$OUTDIR"/{log,tmp}
LOG="$OUTDIR/log/run.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "[INFO] Start: $(date)"
echo "[INFO] H1: $H1_FA"
echo "[INFO] H2: $H2_FA"
echo "[INFO] INV: $INV_TSV"
echo "[INFO] R1 : $R1"
echo "[INFO] R2 : $R2"
echo "[INFO] OUT: $OUTDIR"
echo "[INFO] THREADS=$THREADS MAPQ=$MAPQ FLANK_BP=$FLANK_BP MIN_BP_PAIRS=$MIN_BP_PAIRS MIN_BP_TO_CTRL_FRAC=$MIN_BP_TO_CTRL_FRAC"

JOINT_FA="$OUTDIR/H1H2_joint.fa"
JOINT_FAI="$JOINT_FA.fai"
JOINT_MMI="$OUTDIR/H1H2_joint.mmi"
SUMMARY_TSV="$OUTDIR/inv_joint_hic_support.summary.tsv"
PRIORITY_TSV="$OUTDIR/inv_joint_hic_support.priority.tsv"
CONCLUSION_TXT="$OUTDIR/inv_joint_hic_support.conclusion.txt"

if [[ ! -s "$JOINT_FA" ]]; then
  echo "[INFO] Building joint FASTA: $JOINT_FA"
  awk '/^>/{h=substr($0,2); split(h,a,/[\t ]+/); print ">H1|"a[1]; next} {print}' "$H1_FA" > "$JOINT_FA"
  awk '/^>/{h=substr($0,2); split(h,a,/[\t ]+/); print ">H2|"a[1]; next} {print}' "$H2_FA" >> "$JOINT_FA"
else
  echo "[INFO] Reuse joint FASTA: $JOINT_FA"
fi

if [[ ! -s "$JOINT_FAI" ]]; then
  echo "[INFO] samtools faidx"
  "$SAMTOOLS" faidx "$JOINT_FA"
else
  echo "[INFO] Reuse FAI: $JOINT_FAI"
fi

if [[ ! -s "$JOINT_MMI" ]]; then
  echo "[INFO] minimap2 index"
  "$MINIMAP2" -d "$JOINT_MMI" "$JOINT_FA"
else
  echo "[INFO] Reuse minimap2 index: $JOINT_MMI"
fi

echo "[INFO] Mapping + support counting ..."
"$MINIMAP2" -t "$THREADS" -x sr --secondary=no -a "$JOINT_MMI" "$R1" "$R2" \
  | "$SAMTOOLS" view -@ "$THREADS" -h -q "$MAPQ" -f 1 -F 3340 - \
  | "$SAMTOOLS" collate -@ "$THREADS" -O -u - \
  | "$SAMTOOLS" view -h - \
  | "$PYTHON" - "$INV_TSV" "$JOINT_FAI" "$SUMMARY_TSV" "$PRIORITY_TSV" "$CONCLUSION_TXT" "$FLANK_BP" "$MIN_BP_PAIRS" "$MIN_BP_TO_CTRL_FRAC" "$PRIORITY_SCAFFOLDS" "$MIN_INV_LEN" <<'PY'
import csv
import sys
from collections import defaultdict

if len(sys.argv) != 11:
    raise SystemExit(f"Expected 10 args, got {len(sys.argv)-1}")

inv_tsv = sys.argv[1]
joint_fai = sys.argv[2]
summary_tsv = sys.argv[3]
priority_tsv = sys.argv[4]
conclusion_txt = sys.argv[5]
flank_bp = int(sys.argv[6])
min_bp_pairs = int(sys.argv[7])
min_bp_to_ctrl_frac = float(sys.argv[8])
priority_scaffolds = set(x.strip() for x in sys.argv[9].split(",") if x.strip())
min_inv_len = int(sys.argv[10])

def load_fai(path):
    d = {}
    with open(path) as f:
        for ln in f:
            if not ln.strip():
                continue
            c = ln.rstrip("\n").split("\t")
            d[c[0]] = int(c[1])
    return d

def clip_interval(start, end, chrom_len):
    s = max(1, int(start))
    e = min(int(end), int(chrom_len))
    if s > e:
        return None
    return (s, e)

def in_iv(pos, iv):
    return iv is not None and iv[0] <= pos <= iv[1]

def cross(pos1, pos2, iv_a, iv_b):
    return (in_iv(pos1, iv_a) and in_iv(pos2, iv_b)) or (in_iv(pos1, iv_b) and in_iv(pos2, iv_a))

def ratio(bp, ctrl):
    if ctrl <= 0:
        return 1.0 if bp > 0 else 0.0
    return bp / ctrl

fai = load_fai(joint_fai)

records = []
with open(inv_tsv) as f:
    rd = csv.DictReader(f, delimiter="\t")
    needed = {"RefChr", "RefStart", "RefEnd", "QryChr", "QryStart", "QryEnd", "ID"}
    miss = needed - set(rd.fieldnames or [])
    if miss:
        raise SystemExit(f"Missing columns in inversion table: {sorted(miss)}")
    for row in rd:
        ref_s = int(row["RefStart"])
        ref_e = int(row["RefEnd"])
        qry_s = int(row["QryStart"])
        qry_e = int(row["QryEnd"])
        ref_len = abs(ref_e - ref_s) + 1
        qry_len = abs(qry_e - qry_s) + 1
        min_len = min(ref_len, qry_len)
        if min_len < min_inv_len:
            continue
        rec = {
            "ID": row["ID"],
            "RefChr": row["RefChr"],
            "RefStart": ref_s,
            "RefEnd": ref_e,
            "QryChr": row["QryChr"],
            "QryStart": qry_s,
            "QryEnd": qry_e,
            "RefLen": ref_len,
            "QryLen": qry_len,
            "MinLen": min_len,
            "Priority": (row["RefChr"] in priority_scaffolds) or (row["QryChr"] in priority_scaffolds),
        }
        records.append(rec)

if not records:
    raise SystemExit("No inversion records after filtering.")

windows = {}
regions_by_chr = defaultdict(list)
for rec in records:
    inv_id = rec["ID"]
    h1_chr = f"H1|{rec['RefChr']}"
    h2_chr = f"H2|{rec['QryChr']}"
    if h1_chr not in fai:
        raise SystemExit(f"Chromosome not found in joint reference: {h1_chr}")
    if h2_chr not in fai:
        raise SystemExit(f"Chromosome not found in joint reference: {h2_chr}")

    for hap, chrom, s, e in (
        ("H1", h1_chr, rec["RefStart"], rec["RefEnd"]),
        ("H2", h2_chr, rec["QryStart"], rec["QryEnd"]),
    ):
        clen = fai[chrom]
        left_out = clip_interval(s - flank_bp, s - 1, clen)
        left_in = clip_interval(s, s + flank_bp - 1, clen)
        right_in = clip_interval(e - flank_bp + 1, e, clen)
        right_out = clip_interval(e + 1, e + flank_bp, clen)
        left_ctrl_a = clip_interval(s - 2 * flank_bp, s - flank_bp - 1, clen)
        left_ctrl_b = clip_interval(s - flank_bp, s - 1, clen)
        right_ctrl_a = clip_interval(e + 1, e + flank_bp, clen)
        right_ctrl_b = clip_interval(e + flank_bp + 1, e + 2 * flank_bp, clen)
        windows[(inv_id, hap)] = {
            "chrom": chrom,
            "left_out": left_out,
            "left_in": left_in,
            "right_in": right_in,
            "right_out": right_out,
            "left_ctrl_a": left_ctrl_a,
            "left_ctrl_b": left_ctrl_b,
            "right_ctrl_a": right_ctrl_a,
            "right_ctrl_b": right_ctrl_b,
        }
        regions_by_chr[chrom].append((inv_id, hap))

counts = defaultdict(lambda: {"left_bp": 0, "right_bp": 0, "left_ctrl": 0, "right_ctrl": 0})

def process_group(group):
    if not group:
        return
    r1 = None
    r2 = None
    for flag, chrom, pos in group:
        if (flag & 64) and r1 is None:
            r1 = (chrom, pos)
        elif (flag & 128) and r2 is None:
            r2 = (chrom, pos)
    if r1 is None or r2 is None:
        if len(group) >= 2:
            r1 = (group[0][1], group[0][2])
            r2 = (group[1][1], group[1][2])
        else:
            return
    if r1[0] != r2[0]:
        return
    chrom = r1[0]
    if chrom not in regions_by_chr:
        return
    p1, p2 = r1[1], r2[1]
    for inv_id, hap in regions_by_chr[chrom]:
        w = windows[(inv_id, hap)]
        c = counts[(inv_id, hap)]
        if cross(p1, p2, w["left_out"], w["left_in"]):
            c["left_bp"] += 1
        if cross(p1, p2, w["right_in"], w["right_out"]):
            c["right_bp"] += 1
        if cross(p1, p2, w["left_ctrl_a"], w["left_ctrl_b"]):
            c["left_ctrl"] += 1
        if cross(p1, p2, w["right_ctrl_a"], w["right_ctrl_b"]):
            c["right_ctrl"] += 1

current_qname = None
group = []
pairs_scanned = 0

for ln in sys.stdin:
    if not ln or ln[0] == "@":
        continue
    c = ln.rstrip("\n").split("\t")
    if len(c) < 4:
        continue
    qname = c[0]
    flag = int(c[1])
    chrom = c[2]
    if chrom == "*":
        continue
    pos = int(c[3])
    if current_qname is None:
        current_qname = qname
    if qname != current_qname:
        process_group(group)
        pairs_scanned += 1
        group = []
        current_qname = qname
    group.append((flag, chrom, pos))

if group:
    process_group(group)
    pairs_scanned += 1

def hap_decision(c):
    left_ratio = ratio(c["left_bp"], c["left_ctrl"])
    right_ratio = ratio(c["right_bp"], c["right_ctrl"])
    left_ok = (c["left_bp"] >= min_bp_pairs) and (left_ratio >= min_bp_to_ctrl_frac)
    right_ok = (c["right_bp"] >= min_bp_pairs) and (right_ratio >= min_bp_to_ctrl_frac)
    return left_ok and right_ok, left_ratio, right_ratio, left_ok, right_ok

rows = []
for rec in records:
    inv_id = rec["ID"]
    c1 = counts[(inv_id, "H1")]
    c2 = counts[(inv_id, "H2")]
    h1_ok, h1_lr, h1_rr, h1_left_ok, h1_right_ok = hap_decision(c1)
    h2_ok, h2_lr, h2_rr, h2_left_ok, h2_right_ok = hap_decision(c2)
    independent = h1_ok and h2_ok
    if independent:
        conclusion = "SUPPORTED"
    elif h1_ok or h2_ok:
        conclusion = "PARTIAL"
    else:
        conclusion = "NOT_SUPPORTED"
    rows.append({
        **rec,
        "H1_left_bp_pairs": c1["left_bp"],
        "H1_right_bp_pairs": c1["right_bp"],
        "H1_left_ctrl_pairs": c1["left_ctrl"],
        "H1_right_ctrl_pairs": c1["right_ctrl"],
        "H1_left_ratio": round(h1_lr, 4),
        "H1_right_ratio": round(h1_rr, 4),
        "H1_left_ok": "YES" if h1_left_ok else "NO",
        "H1_right_ok": "YES" if h1_right_ok else "NO",
        "H1_support": "YES" if h1_ok else "NO",
        "H2_left_bp_pairs": c2["left_bp"],
        "H2_right_bp_pairs": c2["right_bp"],
        "H2_left_ctrl_pairs": c2["left_ctrl"],
        "H2_right_ctrl_pairs": c2["right_ctrl"],
        "H2_left_ratio": round(h2_lr, 4),
        "H2_right_ratio": round(h2_rr, 4),
        "H2_left_ok": "YES" if h2_left_ok else "NO",
        "H2_right_ok": "YES" if h2_right_ok else "NO",
        "H2_support": "YES" if h2_ok else "NO",
        "IndependentSupport": "YES" if independent else "NO",
        "Conclusion": conclusion,
        "Priority": "YES" if rec["Priority"] else "NO",
    })

rows.sort(key=lambda r: (r["Priority"] != "YES", -r["MinLen"], r["ID"]))

header = [
    "ID", "Priority", "RefChr", "RefStart", "RefEnd", "RefLen",
    "QryChr", "QryStart", "QryEnd", "QryLen", "MinLen",
    "H1_left_bp_pairs", "H1_right_bp_pairs", "H1_left_ctrl_pairs", "H1_right_ctrl_pairs",
    "H1_left_ratio", "H1_right_ratio", "H1_left_ok", "H1_right_ok", "H1_support",
    "H2_left_bp_pairs", "H2_right_bp_pairs", "H2_left_ctrl_pairs", "H2_right_ctrl_pairs",
    "H2_left_ratio", "H2_right_ratio", "H2_left_ok", "H2_right_ok", "H2_support",
    "IndependentSupport", "Conclusion",
]

with open(summary_tsv, "w", newline="") as f:
    w = csv.DictWriter(f, delimiter="\t", fieldnames=header)
    w.writeheader()
    for r in rows:
        w.writerow({k: r[k] for k in header})

priority_rows = [r for r in rows if r["Priority"] == "YES"]
with open(priority_tsv, "w", newline="") as f:
    w = csv.DictWriter(f, delimiter="\t", fieldnames=header)
    w.writeheader()
    for r in priority_rows:
        w.writerow({k: r[k] for k in header})

total = len(rows)
priority_n = len(priority_rows)
supported = sum(r["Conclusion"] == "SUPPORTED" for r in rows)
partial = sum(r["Conclusion"] == "PARTIAL" for r in rows)
not_supported = sum(r["Conclusion"] == "NOT_SUPPORTED" for r in rows)
pr_supported = sum(r["Conclusion"] == "SUPPORTED" for r in priority_rows)
pr_partial = sum(r["Conclusion"] == "PARTIAL" for r in priority_rows)
pr_not = sum(r["Conclusion"] == "NOT_SUPPORTED" for r in priority_rows)

with open(conclusion_txt, "w") as f:
    f.write("Joint Hi-C inversion support summary\n")
    f.write(f"Pairs scanned (name-collated SAM groups): {pairs_scanned}\n")
    f.write(f"Inversions analyzed: {total}\n")
    f.write(f"Priority scaffolds: {','.join(sorted(priority_scaffolds))}\n")
    f.write(f"Priority inversions: {priority_n}\n")
    f.write("\nOverall conclusion counts:\n")
    f.write(f"  SUPPORTED: {supported}\n")
    f.write(f"  PARTIAL: {partial}\n")
    f.write(f"  NOT_SUPPORTED: {not_supported}\n")
    f.write("\nPriority conclusion counts:\n")
    f.write(f"  SUPPORTED: {pr_supported}\n")
    f.write(f"  PARTIAL: {pr_partial}\n")
    f.write(f"  NOT_SUPPORTED: {pr_not}\n")
    f.write("\nPriority details:\n")
    for r in priority_rows:
        f.write(
            f"  {r['ID']} {r['RefChr']}:{r['RefStart']}-{r['RefEnd']} "
            f"<-> {r['QryChr']}:{r['QryStart']}-{r['QryEnd']} "
            f"=> {r['Conclusion']} (IndependentSupport={r['IndependentSupport']})\n"
        )

print(f"[INFO] pairs_scanned={pairs_scanned}", file=sys.stderr)
print(f"[INFO] wrote {summary_tsv}", file=sys.stderr)
print(f"[INFO] wrote {priority_tsv}", file=sys.stderr)
print(f"[INFO] wrote {conclusion_txt}", file=sys.stderr)
PY

echo "[INFO] Done: $(date)"
echo "[INFO] Summary : $SUMMARY_TSV"
echo "[INFO] Priority: $PRIORITY_TSV"
echo "[INFO] Report  : $CONCLUSION_TXT"

