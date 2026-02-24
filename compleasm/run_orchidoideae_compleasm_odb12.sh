#!/bin/bash
set -euo pipefail

ROOT="/g/data/xf3/zz3507"
OUT_BASE="${OUT_BASE:-${ROOT}/Output/20260127Genome/compleasm_orchidoideae_odb12}"
LINEAGE_NAME="${LINEAGE_NAME:-embryophyta_odb12}"
LINEAGE_ROOT="${LINEAGE_ROOT:-${ROOT}/database/busco_downloads/lineages}"
COMPLEASM_BIN="${COMPLEASM_BIN:-${ROOT}/script/Chiloglottis_genome_2026/compleasm/bin/compleasm_v027}"
MINIPROT_BIN="${MINIPROT_BIN:-/g/data/xf3/miniconda/envs/compleasm/bin/miniprot}"
HMMSEARCH_BIN="${HMMSEARCH_BIN:-/g/data/xf3/miniconda/envs/compleasm/bin/hmmsearch}"
THREADS="${THREADS:-${PBS_NCPUS:-16}}"
MANIFEST="${MANIFEST:-${ROOT}/script/Chiloglottis_genome_2026/compleasm/orchidoideae_odb12.manifest.tsv}"
COMPLEASM_DB_WORKROOT="${COMPLEASM_DB_WORKROOT:-${OUT_BASE}/offline_mb_downloads}"
COMPLEASM_VERSION_STR=""

export PATH="$(dirname "${COMPLEASM_BIN}"):$(dirname "${MINIPROT_BIN}"):$(dirname "${HMMSEARCH_BIN}"):${PATH}"

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
if [[ ! -x "${MINIPROT_BIN}" ]]; then
  echo "ERROR: miniprot not executable: ${MINIPROT_BIN}" >&2
  exit 1
fi
if [[ ! -x "${HMMSEARCH_BIN}" ]]; then
  echo "ERROR: hmmsearch not executable: ${HMMSEARCH_BIN}" >&2
  exit 1
fi

cp -f "${MANIFEST}" "${OUT_BASE}/manifest.tsv"

prepare_compleasm_offline_db() {
  # compleasm always initializes a BUSCO downloader and tries to check file_versions
  # online unless *.done markers exist. We prepare an isolated local "mb_downloads"
  # style directory to force offline execution.
  #
  # v0.2.6 also hardcodes *_odb10 names (legacy mode), while v0.2.7 supports odb12.
  local src_lineage_dir="${LINEAGE_ROOT}/${LINEAGE_NAME}"
  local legacy_mode=0
  local lineage_target_name="${LINEAGE_NAME}"
  local default_lineage_name="eukaryota_odb12"

  if [[ "${COMPLEASM_VERSION_STR}" == *"0.2.6"* ]]; then
    legacy_mode=1
    lineage_target_name="${LINEAGE_NAME}_odb10"
    default_lineage_name="eukaryota_odb10"
  fi

  local alias_dir="${COMPLEASM_DB_WORKROOT}/${lineage_target_name}"
  local default_lineage_dir="${COMPLEASM_DB_WORKROOT}/${default_lineage_name}"

  mkdir -p "${COMPLEASM_DB_WORKROOT}" "${COMPLEASM_DB_WORKROOT}/placement_files"

  # Prevent online version checks: mark file_versions and placement files as already available.
  : > "${COMPLEASM_DB_WORKROOT}/file_versions.tsv"
  : > "${COMPLEASM_DB_WORKROOT}/file_versions.tsv.done"
  : > "${COMPLEASM_DB_WORKROOT}/placement_files.done"
  rm -f "${COMPLEASM_DB_WORKROOT}/file_versions.tsv.tmp" "${COMPLEASM_DB_WORKROOT}/placement_files.tmp"

  # compleasm Downloader.__init__ always requests a default eukaryota lineage first.
  # Mark it as already available to bypass any network attempt; autolineage is not used.
  mkdir -p "${default_lineage_dir}"
  : > "${COMPLEASM_DB_WORKROOT}/${default_lineage_name}.done"

  # Create/refresh lineage directory (or alias) expected by compleasm.
  if [[ -L "${alias_dir}" || -e "${alias_dir}" ]]; then
    rm -rf "${alias_dir}"
  fi
  ln -s "${src_lineage_dir}" "${alias_dir}"
  : > "${COMPLEASM_DB_WORKROOT}/${lineage_target_name}.done"

  echo "[$(date)] Prepared offline compleasm DB workroot: ${COMPLEASM_DB_WORKROOT}"
  if [[ "${legacy_mode}" -eq 1 ]]; then
    echo "[$(date)] Legacy compleasm mode detected (${COMPLEASM_VERSION_STR}); alias ${lineage_target_name} -> ${src_lineage_dir}"
  else
    echo "[$(date)] Native odb12 mode detected (${COMPLEASM_VERSION_STR}); lineage ${lineage_target_name} -> ${src_lineage_dir}"
  fi
}

NOTE_FILE="${OUT_BASE}/notes.txt"
cat > "${NOTE_FILE}" <<'EOF'
Scope: Orchidoideae-focused compleasm completeness comparison (embryophyta_odb12).
Note: Platanthera guangdongensis genome FASTA was not found in local database paths, so it is not included in this run.
Compatibility: workflow prepares a local offline compleasm mb_downloads cache to avoid internet calls on NCI; native odb12 is used with compleasm v0.2.7+ (legacy alias only if older compleasm is selected).
EOF

SUMMARY_TSV="${OUT_BASE}/compleasm_orchidoideae_odb12_summary.tsv"
printf "label\tfasta\tlineage\ttotal_buscos\tcomplete\tcomplete_single\tcomplete_duplicated\tfragmented\tmissing\tC_pct\tS_pct\tD_pct\tF_pct\tM_pct\trun_dir\tstatus\n" > "${SUMMARY_TSV}"

COMPLEASM_VERSION_STR="$(${COMPLEASM_BIN} --version 2>/dev/null || echo unknown)"
echo "[$(date)] compleasm version: ${COMPLEASM_VERSION_STR}"
echo "[$(date)] Output base: ${OUT_BASE}"
echo "[$(date)] Manifest: ${MANIFEST}"
echo "[$(date)] Lineage: ${LINEAGE_NAME} (${LINEAGE_ROOT})"
echo "[$(date)] miniprot: ${MINIPROT_BIN}"
echo "[$(date)] hmmsearch: ${HMMSEARCH_BIN}"
echo "[$(date)] Threads: ${THREADS}"

prepare_compleasm_offline_db

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
  full_table_nested="${run_dir}/${LINEAGE_NAME}/full_table.tsv"
  log_file="${OUT_BASE}/logs/${label}.compleasm.log"

  echo "[$(date)] ===== ${label} ====="
  echo "[$(date)] FASTA: ${fasta}"

  if [[ ! -f "${fasta}" ]]; then
    echo "[$(date)] MISSING FASTA: ${fasta}" | tee -a "${log_file}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tmissing_fasta\n" \
      "${label}" "${fasta_rel}" "${LINEAGE_NAME}" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "${run_dir}" >> "${SUMMARY_TSV}"
    continue
  fi

  if [[ -f "${full_table}" || -f "${full_table_nested}" ]]; then
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
        -L "${COMPLEASM_DB_WORKROOT}" \
        -t "${THREADS}" \
        --miniprot_execute_path "${MINIPROT_BIN}" \
        --hmmsearch_execute_path "${HMMSEARCH_BIN}"
      echo "[$(date)] compleasm finished for ${label}"
    } 2>&1 | tee "${log_file}"
  fi

  if [[ -f "${full_table}" ]]; then
    parse_full_table "${full_table}" "${label}" "${fasta_rel}" "${run_dir}" >> "${SUMMARY_TSV}"
  elif [[ -f "${full_table_nested}" ]]; then
    parse_full_table "${full_table_nested}" "${label}" "${fasta_rel}" "${run_dir}" >> "${SUMMARY_TSV}"
  else
    echo "[$(date)] ERROR: expected ${full_table} or ${full_table_nested} not found" | tee -a "${log_file}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tno_full_table\n" \
      "${label}" "${fasta_rel}" "${LINEAGE_NAME}" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "${run_dir}" >> "${SUMMARY_TSV}"
  fi
done < "${MANIFEST}"

echo "[$(date)] Summary written: ${SUMMARY_TSV}"
column -t -s $'\t' "${SUMMARY_TSV}" | tee "${OUT_BASE}/compleasm_orchidoideae_odb12_summary.pretty.txt" >/dev/null || true
