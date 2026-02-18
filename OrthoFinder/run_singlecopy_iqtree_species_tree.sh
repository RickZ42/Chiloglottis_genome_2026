#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Build a species tree from single-copy OrthoFinder HOGs using MAFFT + IQ-TREE.

Default behavior:
1) Parse N0.tsv from a completed OrthoFinder run.
2) Keep loci where each species has <=1 copy and occupancy >= --min-occupancy.
3) Export per-locus FASTA with species names as sequence headers.
4) Align each locus with MAFFT.
5) Concatenate alignments and run partitioned IQ-TREE.

Usage:
  bash run_singlecopy_iqtree_species_tree.sh [options]

Options:
  --dataset-dir PATH     Dataset dir with *.faa and OrthoFinder folder.
                         default: /g/data/xf3/zz3507/database/GeneFamilyExpentionContractionDatasetH12026
  --n0-tsv PATH          Root HOG table (N0.tsv).
                         default: <dataset>/OrthoFinder/Results_Feb16/Phylogenetic_Hierarchical_Orthogroups/N0.tsv
  --outdir PATH          Output directory.
                         default: <dataset>/OrthoFinder/Results_Feb16/Species_Tree/singlecopy_iqtree_YYYYMMDD_HHMMSS
  --threads INT          Threads for MAFFT and IQ-TREE. default: 24
  --bootstrap INT        Ultrafast bootstrap replicates. default: 1000
  --min-occupancy FLOAT  Min fraction of species present per locus (0-1).
                         default: 1.0 (strict all-species single-copy)
  --seed INT             Random seed for IQ-TREE. default: 20260218
  --redo                 Remove prior IQ-TREE files for same prefix.
  --help                 Show this help.
USAGE
}

DATASET_DIR="/g/data/xf3/zz3507/database/GeneFamilyExpentionContractionDatasetH12026"
N0_TSV=""
OUTDIR=""
THREADS=24
BOOTSTRAP=1000
MIN_OCCUPANCY="1.0"
SEED=20260218
REDO=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset-dir) DATASET_DIR="$2"; shift 2 ;;
    --n0-tsv) N0_TSV="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --bootstrap) BOOTSTRAP="$2"; shift 2 ;;
    --min-occupancy) MIN_OCCUPANCY="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --redo) REDO=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$N0_TSV" ]]; then
  N0_TSV="${DATASET_DIR}/OrthoFinder/Results_Feb16/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
fi
if [[ -z "$OUTDIR" ]]; then
  OUTDIR="${DATASET_DIR}/OrthoFinder/Results_Feb16/Species_Tree/singlecopy_iqtree_$(date +%Y%m%d_%H%M%S)"
fi

ORTHO_ENV="/g/data/xf3/miniconda/envs/orthofinder"
MAFFT="${ORTHO_ENV}/bin/mafft"
IQTREE="${ORTHO_ENV}/bin/iqtree2"
PYTHON="${ORTHO_ENV}/bin/python"
[[ -x "$PYTHON" ]] || PYTHON="$(command -v python3)"

[[ -s "$N0_TSV" ]] || { echo "ERROR: missing N0.tsv: $N0_TSV" >&2; exit 1; }
[[ -d "$DATASET_DIR" ]] || { echo "ERROR: missing dataset dir: $DATASET_DIR" >&2; exit 1; }
[[ -x "$MAFFT" ]] || { echo "ERROR: MAFFT not executable: $MAFFT" >&2; exit 1; }
[[ -x "$IQTREE" ]] || { echo "ERROR: IQ-TREE not executable: $IQTREE" >&2; exit 1; }
[[ -x "$PYTHON" ]] || { echo "ERROR: python not found" >&2; exit 1; }

mkdir -p "$OUTDIR"/{logs,tmp,single_copy_fastas,alignments,iqtree}
LOG="$OUTDIR/logs/run_singlecopy_iqtree.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "[INFO] start: $(date)"
echo "[INFO] DATASET_DIR=$DATASET_DIR"
echo "[INFO] N0_TSV=$N0_TSV"
echo "[INFO] OUTDIR=$OUTDIR"
echo "[INFO] THREADS=$THREADS BOOTSTRAP=$BOOTSTRAP MIN_OCCUPANCY=$MIN_OCCUPANCY SEED=$SEED"

SELECTED_TSV="$OUTDIR/selected_hogs.tsv"
SPECIES_LIST="$OUTDIR/species.list.txt"
FASTA_DIR="$OUTDIR/single_copy_fastas"
ALIGN_DIR="$OUTDIR/alignments"

echo "[INFO] selecting single-copy loci and exporting FASTA"
"$PYTHON" - "$DATASET_DIR" "$N0_TSV" "$MIN_OCCUPANCY" "$SELECTED_TSV" "$SPECIES_LIST" "$FASTA_DIR" <<'PY'
import csv
import os
import re
import sys
from collections import defaultdict

dataset_dir, n0_tsv, min_occ_s, selected_tsv, species_list_fp, fasta_dir = sys.argv[1:7]
min_occ = float(min_occ_s)
if not (0.0 < min_occ <= 1.0):
    raise SystemExit(f"min-occupancy must be in (0,1], got {min_occ}")

def fasta_ids(path):
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                yield line[1:].strip().split()[0]

def parse_cell(cell, valid_ids):
    if not cell:
        return []
    hits = []
    for part in re.split(r",\s+", cell):
        token = part.strip().split()[0] if part.strip() else ""
        if token in valid_ids:
            hits.append(token)
    dedup = []
    seen = set()
    for x in hits:
        if x not in seen:
            dedup.append(x)
            seen.add(x)
    return dedup

with open(n0_tsv) as f:
    header = f.readline().rstrip("\n").split("\t")
species = header[3:]
if not species:
    raise SystemExit("No species columns found in N0.tsv")

os.makedirs(fasta_dir, exist_ok=True)
with open(species_list_fp, "w") as out:
    for sp in species:
        out.write(sp + "\n")

idsets = {}
for sp in species:
    fa = os.path.join(dataset_dir, f"{sp}.faa")
    if not os.path.exists(fa):
        raise SystemExit(f"Missing species FASTA: {fa}")
    idsets[sp] = set(fasta_ids(fa))

selected = []
needed_ids = defaultdict(set)
total_rows = 0
with open(n0_tsv) as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    for row in reader:
        total_rows += 1
        hog = row[0]
        og = row[1] if len(row) > 1 else ""
        chosen = {}
        multi_copy = False
        present = 0
        for i, sp in enumerate(species, start=3):
            ids = parse_cell(row[i].strip(), idsets[sp])
            if len(ids) > 1:
                multi_copy = True
                break
            if len(ids) == 1:
                chosen[sp] = ids[0]
                present += 1
        if multi_copy:
            continue
        occ = present / len(species)
        if present > 0 and occ >= min_occ:
            selected.append((hog, og, occ, chosen))
            for sp, sid in chosen.items():
                needed_ids[sp].add(sid)

if not selected:
    raise SystemExit("No loci passed filters; lower --min-occupancy or inspect N0.tsv")

seqs = {sp: {} for sp in species}
for sp in species:
    targets = needed_ids[sp]
    if not targets:
        continue
    fa = os.path.join(dataset_dir, f"{sp}.faa")
    with open(fa) as f:
        cur_id = None
        cur_seq = []
        for line in f:
            if line.startswith(">"):
                if cur_id in targets:
                    seqs[sp][cur_id] = "".join(cur_seq)
                cur_id = line[1:].strip().split()[0]
                cur_seq = []
            else:
                if cur_id in targets:
                    cur_seq.append(line.strip())
        if cur_id in targets:
            seqs[sp][cur_id] = "".join(cur_seq)
    missing = targets - set(seqs[sp])
    if missing:
        miss = ",".join(sorted(list(missing))[:5])
        raise SystemExit(f"Missing {len(missing)} target sequences for {sp}; first IDs: {miss}")

with open(selected_tsv, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["HOG", "OG", "present_species", "occupancy"])
    for hog, og, occ, chosen in selected:
        w.writerow([hog, og, len(chosen), f"{occ:.6f}"])
        fp = os.path.join(fasta_dir, f"{hog}.fa")
        with open(fp, "w") as fa_out:
            for sp in species:
                sid = chosen.get(sp)
                if sid is None:
                    continue
                fa_out.write(f">{sp}\n{seqs[sp][sid]}\n")

print(f"total_hogs\t{total_rows}")
print(f"selected_hogs\t{len(selected)}")
print(f"species\t{len(species)}")
PY

N_HOGS=$(awk 'NR>1{n++}END{print n+0}' "$SELECTED_TSV")
echo "[INFO] selected loci: $N_HOGS"
if [[ "$N_HOGS" -lt 1 ]]; then
  echo "ERROR: no selected loci" >&2
  exit 1
fi

echo "[INFO] running MAFFT per locus"
mapfile -t HOGS < <(awk 'NR>1{print $1}' "$SELECTED_TSV")
for HOG in "${HOGS[@]}"; do
  IN_FA="$FASTA_DIR/${HOG}.fa"
  OUT_ALN="$ALIGN_DIR/${HOG}.aln.fa"
  "$MAFFT" --quiet --auto --thread "$THREADS" "$IN_FA" > "$OUT_ALN"
done

SUPERMATRIX="$OUTDIR/singlecopy_supermatrix.fa"
PARTITIONS="$OUTDIR/singlecopy_partitions.nex"

echo "[INFO] concatenating alignments"
"$PYTHON" - "$SPECIES_LIST" "$SELECTED_TSV" "$ALIGN_DIR" "$SUPERMATRIX" "$PARTITIONS" <<'PY'
import csv
import os
import sys

species_list_fp, selected_tsv, align_dir, supermatrix_fp, partitions_fp = sys.argv[1:6]

with open(species_list_fp) as f:
    species = [x.strip() for x in f if x.strip()]
if not species:
    raise SystemExit("No species found in species list")

hogs = []
with open(selected_tsv) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        hogs.append(row["HOG"])
if not hogs:
    raise SystemExit("No HOGs in selected_hogs.tsv")

concat = {sp: [] for sp in species}
start = 1
parts = []

for hog in hogs:
    aln_fp = os.path.join(align_dir, f"{hog}.aln.fa")
    if not os.path.exists(aln_fp):
        raise SystemExit(f"Missing alignment: {aln_fp}")
    seqs = {}
    cur = None
    with open(aln_fp) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                cur = line[1:].split()[0]
                if cur in seqs:
                    raise SystemExit(f"Duplicate sequence '{cur}' in {aln_fp}")
                seqs[cur] = []
            else:
                if cur is None:
                    raise SystemExit(f"FASTA parse error in {aln_fp}")
                seqs[cur].append(line)
    seqs = {k: "".join(v) for k, v in seqs.items()}
    lengths = {len(v) for v in seqs.values()}
    if len(lengths) != 1:
        raise SystemExit(f"Non-uniform alignment length in {aln_fp}")
    if not lengths:
        raise SystemExit(f"Empty alignment: {aln_fp}")
    L = lengths.pop()

    for sp in species:
        concat[sp].append(seqs.get(sp, "-" * L))
    end = start + L - 1
    parts.append((hog, start, end))
    start = end + 1

with open(supermatrix_fp, "w") as out:
    for sp in species:
        out.write(f">{sp}\n{''.join(concat[sp])}\n")

with open(partitions_fp, "w") as out:
    out.write("#nexus\nbegin sets;\n")
    for hog, s, e in parts:
        out.write(f"  charset {hog} = {s}-{e};\n")
    out.write("end;\n")

print(f"concat_loci\t{len(parts)}")
print(f"concat_length\t{parts[-1][2] if parts else 0}")
PY

IQ_PREFIX="$OUTDIR/iqtree/singlecopy_species_tree"
if [[ "$REDO" -eq 1 ]]; then
  rm -f "${IQ_PREFIX}".*
fi

echo "[INFO] running IQ-TREE"
"$IQTREE" \
  -s "$SUPERMATRIX" \
  -p "$PARTITIONS" \
  -m MFP+MERGE \
  -B "$BOOTSTRAP" \
  --alrt 1000 \
  --bnni \
  -nt "$THREADS" \
  -seed "$SEED" \
  -pre "$IQ_PREFIX"

echo "[INFO] done: $(date)"
echo "[INFO] treefile: ${IQ_PREFIX}.treefile"
echo "[INFO] support tree: ${IQ_PREFIX}.contree"
echo "[INFO] report: ${IQ_PREFIX}.iqtree"
