#!/bin/bash
set -euo pipefail

ROOT="/g/data/xf3/zz3507"
OUT_BASE="${OUT_BASE:-${ROOT}/Output/20260127Genome/compleasm_orchidoideae_odb12}"
LINEAGE_NAME="${LINEAGE_NAME:-embryophyta_odb12}"
LINEAGE_ROOT="${LINEAGE_ROOT:-${ROOT}/database/busco_downloads/lineages}"
COMPLEASM_BIN="${COMPLEASM_BIN:-/g/data/xf3/miniconda/envs/compleasm/bin/compleasm}"
THREADS="${THREADS:-${PBS_NCPUS:-16}}"
MANIFEST="${MANIFEST:-${ROOT}/script/Chiloglottis_genome_2026/compleasm/orchidoideae_odb12.manifest.tsv}"

export PATH="$(dirname "${COMPLEASM_BIN}"):${PATH}"

mkdir -p "${OUT_BASE}/runs" "${OUT_BASE}/logs"

if [[ ! -x "${COMPLEASM_BIN}" ]]; then
  echo "ERROR: compleasm not executable: ${COMPLEASM_BIN}" >&2
  exit 1
fi
if [[ ! -d "${LINEAGE_ROOT}/${LINEAGE_NAME}" ]]; then
  echo "ERROR: lineage not found: ${LINEAGE_ROOT}/${LINEAGE_NAME}" >&2
  exit 1
fi
if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: manifest not found: ${MANIFEST}" >&2
  exit 1
fi

cp -f "${MANIFEST}" "${OUT_BASE}/manifest.tsv"

NOTE_FILE="${OUT_BASE}/notes.txt"
cat > "${NOTE_FILE}" <<'EOF'
Scope: Orchidoideae-focused compleasm completeness comparison (embryophyta_odb12).
Note: Platanthera guangdongensis genome FASTA was not found in local database paths, so it is not included in this run.
EOF

SUMMARY_TSV="${OUT_BASE}/compleasm_orchidoideae_odb12_summary.tsv"
printf "label\tfasta\tlineage\ttotal_buscos\tcomplete\tcomplete_single\tcomplete_duplicated\tfragmented\tmissing\tC_pct\tS_pct\tD_pct\tF_pct\tM_pct\trun_dir\tstatus\n" > "${SUMMARY_TSV}"

echo "[$(date)] compleasm version: $(${COMPLEASM_BIN} --version 2>/dev/null || echo unknown)"
echo "[$(date)] Output base: ${OUT_BASE}"
echo "[$(date)] Manifest: ${MANIFEST}"
echo "[$(date)] Lineage: ${LINEAGE_NAME} (${LINEAGE_ROOT})"
echo "[$(date)] Threads: ${THREADS}"

parse_full_table() {
  local full_table="$1"
  local label="$2"
  local fasta="$3"
  local run_dir="$4"

  python3 - "$full_table" "$label" "$fasta" "$LINEAGE_NAME" "$run_dir" <<'PY'
import csv, sys
from collections import Counter

full_table, label, fasta, lineage_name, run_dir = sys.argv[1:]

rank = {"Missing": 0, "Fragmented": 1, "Single": 2, "Complete": 2, "Duplicated": 3}
best = {}
with open(full_table, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        gene = row.get("Gene") or row.get("Busco id")
        status = row["Status"]
        prev = best.get(gene)
        if prev is None or rank.get(status, -1) > rank.get(prev, -1):
            best[gene] = status

cnt = Counter(best.values())
total = len(best)
single = cnt.get("Single", 0) + cnt.get("Complete", 0)
dup = cnt.get("Duplicated", 0)
frag = cnt.get("Fragmented", 0)
miss = cnt.get("Missing", 0)
comp = single + dup

def pct(x):
    return f"{(100.0 * x / total):.2f}" if total else "0.00"

print("\t".join([
    label, fasta, lineage_name, str(total), str(comp), str(single), str(dup),
    str(frag), str(miss), pct(comp), pct(single), pct(dup), pct(frag), pct(miss),
    run_dir, "ok"
]))
PY
}

while IFS=$'\t' read -r label fasta_rel; do
  [[ -z "${label}" ]] && continue
  [[ "${label}" == "label" ]] && continue

  fasta="${ROOT}/${fasta_rel}"
  run_dir="${OUT_BASE}/runs/${label}"
  full_table="${run_dir}/full_table.tsv"
  log_file="${OUT_BASE}/logs/${label}.compleasm.log"

  echo "[$(date)] ===== ${label} ====="
  echo "[$(date)] FASTA: ${fasta}"

  if [[ ! -f "${fasta}" ]]; then
    echo "[$(date)] MISSING FASTA: ${fasta}" | tee -a "${log_file}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tmissing_fasta\n" \
      "${label}" "${fasta_rel}" "${LINEAGE_NAME}" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "${run_dir}" >> "${SUMMARY_TSV}"
    continue
  fi

  if [[ -f "${full_table}" ]]; then
    echo "[$(date)] Reusing existing compleasm output: ${run_dir}" | tee -a "${log_file}"
  else
    rm -rf "${run_dir}"
    mkdir -p "${run_dir}"
    {
      echo "[$(date)] Running compleasm for ${label}"
      "${COMPLEASM_BIN}" run \
        -a "${fasta}" \
        -o "${run_dir}" \
        -l "${LINEAGE_NAME}" \
        -L "${LINEAGE_ROOT}" \
        -t "${THREADS}"
      echo "[$(date)] compleasm finished for ${label}"
    } 2>&1 | tee "${log_file}"
  fi

  if [[ -f "${full_table}" ]]; then
    parse_full_table "${full_table}" "${label}" "${fasta_rel}" "${run_dir}" >> "${SUMMARY_TSV}"
  else
    echo "[$(date)] ERROR: expected ${full_table} not found" | tee -a "${log_file}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tno_full_table\n" \
      "${label}" "${fasta_rel}" "${LINEAGE_NAME}" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "${run_dir}" >> "${SUMMARY_TSV}"
  fi
done < "${MANIFEST}"

echo "[$(date)] Summary written: ${SUMMARY_TSV}"
column -t -s $'\t' "${SUMMARY_TSV}" | tee "${OUT_BASE}/compleasm_orchidoideae_odb12_summary.pretty.txt" >/dev/null || true

