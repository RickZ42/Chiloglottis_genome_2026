#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Build a supermatrix species tree from near-single-copy orthogroups using an existing
OrthoFinder -M msa result (reuses OrthoFinder per-OG alignments).

Near-single-copy definition (default):
- occupancy >= 0.80
- no species has >2 copies in an OG
- no more than 4 species are multi-copy in an OG
- duplicated species are reduced to one representative sequence (longest ungapped MSA seq)

Usage:
  bash run_near_singlecopy_iqtree_from_orthofinder_msa.sh [options]

Options:
  --results-dir PATH          OrthoFinder results dir (contains Orthogroups/, MultipleSequenceAlignments/).
                              default: /g/data/xf3/zz3507/database/GeneFamilyExpentionContractionDatasetH22026/OrthoFinder/Results_Feb15
  --outdir PATH               Output dir.
                              default: <results-dir>/Species_Tree/near_singlecopy_iqtree_<timestamp>
  --min-occupancy FLOAT       Min fraction of species present. default: 0.8
  --max-copy INT              Max copies allowed for any one species in an OG. default: 2
  --max-multicopy-species INT Max number of species with >1 copy in an OG. default: 4
  --drop-all-gap-cols         Drop all-gap columns per OG after duplicate pruning (default: on)
  --keep-all-gap-cols         Keep all columns exactly as in OrthoFinder MSA
  --threads INT               IQ-TREE threads. default: 4
  --iqtree-model STR          IQ-TREE model string. default: MFP+MERGE
  --bootstrap INT             UFBoot replicates. default: 1000
  --alrt INT                  SH-aLRT replicates. default: 1000
  --seed INT                  IQ-TREE seed. default: 20260223
  --ml-only                   Run IQ-TREE ML tree search only (no support tests)
  --compare-tree PATH         Compare final tree against a reference Newick tree (RF-like split comparison)
  --redo                      Remove prior IQ-TREE files for same prefix
  --help                      Show help
USAGE
}

RESULTS_DIR="/g/data/xf3/zz3507/database/GeneFamilyExpentionContractionDatasetH22026/OrthoFinder/Results_Feb15"
OUTDIR=""
MIN_OCCUPANCY="0.8"
MAX_COPY=2
MAX_MULTI=4
DROP_ALL_GAP_COLS=1
THREADS=4
IQTREE_MODEL="MFP+MERGE"
BOOTSTRAP=1000
ALRT=1000
SEED=20260223
ML_ONLY=0
COMPARE_TREE=""
REDO=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --min-occupancy) MIN_OCCUPANCY="$2"; shift 2 ;;
    --max-copy) MAX_COPY="$2"; shift 2 ;;
    --max-multicopy-species) MAX_MULTI="$2"; shift 2 ;;
    --drop-all-gap-cols) DROP_ALL_GAP_COLS=1; shift ;;
    --keep-all-gap-cols) DROP_ALL_GAP_COLS=0; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --iqtree-model) IQTREE_MODEL="$2"; shift 2 ;;
    --bootstrap) BOOTSTRAP="$2"; shift 2 ;;
    --alrt) ALRT="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --ml-only) ML_ONLY=1; shift ;;
    --compare-tree) COMPARE_TREE="$2"; shift 2 ;;
    --redo) REDO=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

OG_COUNT_TSV="$RESULTS_DIR/Orthogroups/Orthogroups.GeneCount.tsv"
MSA_DIR="$RESULTS_DIR/MultipleSequenceAlignments"
[[ -s "$OG_COUNT_TSV" ]] || { echo "ERROR: missing $OG_COUNT_TSV" >&2; exit 1; }
[[ -d "$MSA_DIR" ]] || { echo "ERROR: missing $MSA_DIR" >&2; exit 1; }

if [[ -z "$OUTDIR" ]]; then
  OUTDIR="$RESULTS_DIR/Species_Tree/near_singlecopy_iqtree_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$OUTDIR"/{logs,iqtree,tmp}

ORTHO_ENV="/g/data/xf3/miniconda/envs/orthofinder"
IQTREE="$ORTHO_ENV/bin/iqtree2"
PYTHON="$ORTHO_ENV/bin/python"
[[ -x "$PYTHON" ]] || PYTHON="$(command -v python3)"
[[ -x "$IQTREE" ]] || { echo "ERROR: missing IQ-TREE at $IQTREE" >&2; exit 1; }
[[ -x "$PYTHON" ]] || { echo "ERROR: python not found" >&2; exit 1; }

LOG="$OUTDIR/logs/run_near_singlecopy_iqtree.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "[INFO] start: $(date)"
echo "[INFO] RESULTS_DIR=$RESULTS_DIR"
echo "[INFO] OUTDIR=$OUTDIR"
echo "[INFO] MIN_OCCUPANCY=$MIN_OCCUPANCY MAX_COPY=$MAX_COPY MAX_MULTI=$MAX_MULTI DROP_ALL_GAP_COLS=$DROP_ALL_GAP_COLS"
echo "[INFO] THREADS=$THREADS IQTREE_MODEL=$IQTREE_MODEL BOOTSTRAP=$BOOTSTRAP ALRT=$ALRT SEED=$SEED ML_ONLY=$ML_ONLY"
if [[ -n "$COMPARE_TREE" ]]; then
  echo "[INFO] COMPARE_TREE=$COMPARE_TREE"
fi

SELECTED_STATS="$OUTDIR/selected_near_singlecopy_ogs.tsv"
SPECIES_LIST="$OUTDIR/species.list.txt"
SUPERMATRIX="$OUTDIR/near_singlecopy_supermatrix.fa"
PARTITIONS="$OUTDIR/near_singlecopy_partitions.nex"

echo "[INFO] selecting OGs and building concatenated supermatrix from OrthoFinder MSAs"
"$PYTHON" - "$OG_COUNT_TSV" "$MSA_DIR" "$MIN_OCCUPANCY" "$MAX_COPY" "$MAX_MULTI" "$DROP_ALL_GAP_COLS" \
  "$SELECTED_STATS" "$SPECIES_LIST" "$SUPERMATRIX" "$PARTITIONS" <<'PY'
import csv
import os
import sys

gene_count_tsv, msa_dir, min_occ_s, max_copy_s, max_multi_s, drop_gap_cols_s, selected_stats_fp, species_list_fp, supermatrix_fp, partitions_fp = sys.argv[1:11]
min_occ = float(min_occ_s)
max_copy = int(max_copy_s)
max_multi = int(max_multi_s)
drop_gap_cols = bool(int(drop_gap_cols_s))

if not (0 < min_occ <= 1):
    raise SystemExit(f"Invalid min occupancy: {min_occ}")
if max_copy < 1 or max_multi < 0:
    raise SystemExit("Invalid max-copy or max-multicopy-species")

def read_fasta(fp):
    seqs = []
    name = None
    buf = []
    with open(fp) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name = line[1:].strip()
                buf = []
            else:
                if name is None:
                    raise ValueError(f"FASTA parse error in {fp}")
                buf.append(line.strip())
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs

with open(gene_count_tsv) as f:
    r = csv.reader(f, delimiter="\t")
    header = next(r)
    species = header[1:-1]
    rows = [row for row in r]

species_sorted = sorted(species, key=len, reverse=True)
nsp = len(species)
if nsp == 0:
    raise SystemExit("No species in Orthogroups.GeneCount.tsv")

for sp in species:
    if not sp:
        raise SystemExit("Empty species name in header")

for sp in species:
    pass
with open(species_list_fp, "w") as out:
    for sp in species:
        out.write(sp + "\n")

concat = {sp: [] for sp in species}
parts = []
selected_rows = []
start = 1
total_ogs = 0
selected_ogs = 0

for row in rows:
    total_ogs += 1
    og = row[0]
    counts = [int(x) for x in row[1:-1]]
    present = sum(c > 0 for c in counts)
    if present == 0:
        continue
    occ = present / nsp
    if occ < min_occ:
        continue
    if any(c > max_copy for c in counts):
        continue
    multi_species = sum(c > 1 for c in counts)
    if multi_species > max_multi:
        continue

    msa_fp = os.path.join(msa_dir, f"{og}.fa")
    if not os.path.exists(msa_fp):
        continue
    records = read_fasta(msa_fp)
    if not records:
        continue

    per_sp = {sp: [] for sp in species}
    unknown = []
    for header_full, seq in records:
        matched = None
        for sp in species_sorted:
            prefix = sp + "_"
            if header_full.startswith(prefix):
                matched = sp
                break
        if matched is None:
            unknown.append(header_full)
            continue
        per_sp[matched].append((header_full, seq))
    if unknown:
        raise SystemExit(f"Could not map {len(unknown)} MSA headers to species in {msa_fp}; first={unknown[0]}")

    # Choose one representative per species by longest ungapped aligned sequence.
    chosen = {}
    aln_len = None
    dup_pruned = 0
    for sp in species:
        recs = per_sp[sp]
        if not recs:
            continue
        if aln_len is None:
            aln_len = len(recs[0][1])
        for _, s in recs:
            if len(s) != aln_len:
                raise SystemExit(f"Non-uniform alignment length in {msa_fp}")
        if len(recs) > 1:
            dup_pruned += len(recs) - 1
        recs_sorted = sorted(recs, key=lambda x: (sum(ch != "-" for ch in x[1]), x[0]), reverse=True)
        chosen[sp] = recs_sorted[0][1]

    if aln_len is None or aln_len == 0:
        continue

    # Optional: remove columns that are all-gap after duplicate pruning / missing taxa.
    keep_cols = None
    if drop_gap_cols:
        keep_cols = []
        for i in range(aln_len):
            all_gap = True
            for sp in species:
                s = chosen.get(sp)
                if s is not None and s[i] != "-":
                    all_gap = False
                    break
            if not all_gap:
                keep_cols.append(i)
        if not keep_cols:
            continue
        trimmed_len = len(keep_cols)
    else:
        trimmed_len = aln_len

    for sp in species:
        seq = chosen.get(sp)
        if seq is None:
            concat[sp].append("-" * trimmed_len)
        else:
            if keep_cols is None:
                concat[sp].append(seq)
            else:
                concat[sp].append("".join(seq[i] for i in keep_cols))

    end = start + trimmed_len - 1
    parts.append((og, start, end))
    selected_rows.append((og, present, occ, multi_species, dup_pruned, aln_len, trimmed_len))
    start = end + 1
    selected_ogs += 1

if selected_ogs == 0:
    raise SystemExit("No OGs selected. Relax thresholds.")

with open(selected_stats_fp, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["Orthogroup", "present_species", "occupancy", "multicopy_species", "duplicates_pruned", "msa_len", "used_len"])
    for og, present, occ, multi, pruned, raw_len, used_len in selected_rows:
        w.writerow([og, present, f"{occ:.6f}", multi, pruned, raw_len, used_len])

with open(supermatrix_fp, "w") as out:
    for sp in species:
        out.write(f">{sp}\n{''.join(concat[sp])}\n")

with open(partitions_fp, "w") as out:
    out.write("#nexus\nbegin sets;\n")
    for og, s, e in parts:
        out.write(f"  charset {og} = {s}-{e};\n")
    out.write("end;\n")

print(f"species\t{nsp}")
print(f"total_orthogroups\t{total_ogs}")
print(f"selected_orthogroups\t{selected_ogs}")
print(f"concat_length\t{parts[-1][2]}")
print(f"thresholds\tocc>={min_occ};maxcopy<={max_copy};maxmulti<={max_multi}")
PY

N_OGS=$(awk 'NR>1{n++}END{print n+0}' "$SELECTED_STATS")
CONCAT_LEN=$(awk 'NR==2{print length($0)}' "$SUPERMATRIX")
echo "[INFO] selected OGs: $N_OGS"
echo "[INFO] supermatrix length (first seq): $CONCAT_LEN"

IQDIR="$OUTDIR/iqtree"
mkdir -p "$IQDIR"
IQ_PREFIX="$IQDIR/near_singlecopy_species_tree"
if [[ "$REDO" -eq 1 ]]; then
  rm -f "${IQ_PREFIX}".*
fi

echo "[INFO] running IQ-TREE"
if [[ "$ML_ONLY" -eq 1 ]]; then
  "$IQTREE" \
    -s "$SUPERMATRIX" \
    -p "$PARTITIONS" \
    -m "$IQTREE_MODEL" \
    -nt "$THREADS" \
    -seed "$SEED" \
    -pre "$IQ_PREFIX"
else
  "$IQTREE" \
    -s "$SUPERMATRIX" \
    -p "$PARTITIONS" \
    -m "$IQTREE_MODEL" \
    -B "$BOOTSTRAP" \
    --alrt "$ALRT" \
    --bnni \
    -nt "$THREADS" \
    -seed "$SEED" \
    -pre "$IQ_PREFIX"
fi

TREE_FOR_COMPARE="${IQ_PREFIX}.treefile"
if [[ ! -s "$TREE_FOR_COMPARE" ]]; then
  echo "ERROR: IQ-TREE did not produce $TREE_FOR_COMPARE" >&2
  exit 1
fi

if [[ -n "$COMPARE_TREE" ]]; then
  [[ -s "$COMPARE_TREE" ]] || { echo "ERROR: compare tree missing: $COMPARE_TREE" >&2; exit 1; }
  CMP_OUT="$OUTDIR/tree_topology_comparison.tsv"
  echo "[INFO] comparing topology to reference tree"
  "$PYTHON" - "$TREE_FOR_COMPARE" "$COMPARE_TREE" "$CMP_OUT" <<'PY'
import sys
from pathlib import Path

new_fp, ref_fp, out_fp = sys.argv[1:4]

try:
    from Bio import Phylo
except Exception as e:
    raise SystemExit(f"Biopython required for topology comparison: {e}")

def get_term_names(tree):
    return [t.name for t in tree.get_terminals()]

def split_set(tree):
    taxa = get_term_names(tree)
    all_taxa = set(taxa)
    n = len(all_taxa)
    if n < 4:
        return set(), all_taxa
    splits = set()
    root = tree.root
    def rec(clade):
        if clade.is_terminal():
            return {clade.name}
        leaves = set()
        for ch in clade.clades:
            leaves |= rec(ch)
        if clade is not root:
            k = len(leaves)
            if 1 < k < n - 1:
                side = frozenset(leaves if k <= n - k else (all_taxa - leaves))
                splits.add(side)
        return leaves
    rec(root)
    return splits, all_taxa

new_tree = Phylo.read(new_fp, "newick")
ref_tree = Phylo.read(ref_fp, "newick")
new_splits, new_taxa = split_set(new_tree)
ref_splits, ref_taxa = split_set(ref_tree)

with open(out_fp, "w") as out:
    out.write("metric\tvalue\n")
    if new_taxa != ref_taxa:
        out.write(f"taxa_match\tFalse\n")
        out.write(f"new_taxa_n\t{len(new_taxa)}\n")
        out.write(f"ref_taxa_n\t{len(ref_taxa)}\n")
        only_new_taxa = sorted(new_taxa - ref_taxa)
        only_ref_taxa = sorted(ref_taxa - new_taxa)
        out.write(f"only_in_new_taxa\t{','.join(only_new_taxa)}\n")
        out.write(f"only_in_ref_taxa\t{','.join(only_ref_taxa)}\n")
        print("taxa_match\tFalse")
        sys.exit(0)

    shared = new_splits & ref_splits
    only_new = sorted(new_splits - ref_splits, key=lambda s: (len(s), sorted(s)))
    only_ref = sorted(ref_splits - new_splits, key=lambda s: (len(s), sorted(s)))
    rf = len(only_new) + len(only_ref)
    max_rf = len(new_splits) + len(ref_splits)
    norm_rf = (rf / max_rf) if max_rf else 0.0

    out.write("taxa_match\tTrue\n")
    out.write(f"taxa_n\t{len(new_taxa)}\n")
    out.write(f"new_splits\t{len(new_splits)}\n")
    out.write(f"ref_splits\t{len(ref_splits)}\n")
    out.write(f"shared_splits\t{len(shared)}\n")
    out.write(f"rf_like_distance\t{rf}\n")
    out.write(f"normalized_rf_like\t{norm_rf:.6f}\n")
    for i, s in enumerate(only_new[:10], start=1):
        out.write(f"only_new_split_{i}\t{'|'.join(sorted(s))}\n")
    for i, s in enumerate(only_ref[:10], start=1):
        out.write(f"only_ref_split_{i}\t{'|'.join(sorted(s))}\n")

    print("taxa_match\tTrue")
    print(f"taxa_n\t{len(new_taxa)}")
    print(f"shared_splits\t{len(shared)}")
    print(f"new_splits\t{len(new_splits)}")
    print(f"ref_splits\t{len(ref_splits)}")
    print(f"rf_like_distance\t{rf}")
    print(f"normalized_rf_like\t{norm_rf:.6f}")
PY
fi

echo "[INFO] done: $(date)"
echo "[INFO] selected OG stats: $SELECTED_STATS"
echo "[INFO] supermatrix: $SUPERMATRIX"
echo "[INFO] partitions: $PARTITIONS"
echo "[INFO] IQ-TREE report: ${IQ_PREFIX}.iqtree"
echo "[INFO] IQ-TREE treefile: ${IQ_PREFIX}.treefile"
if [[ -s "${IQ_PREFIX}.contree" ]]; then
  echo "[INFO] IQ-TREE contree: ${IQ_PREFIX}.contree"
fi
