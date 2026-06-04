#!/usr/bin/env bash
# Fetch SRA FASTQs for accessions listed in samples/*.csv and align to SARS-CoV-2
# with minimap2. Intended to run inside micromamba env freyja-sc2 (sra-tools, minimap2, samtools).
set -euo pipefail

ENV_NAME="${ENV_NAME:-freyja-sc2}"
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 4)}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

resolve_workdir() {
  if [[ -n "${WORKDIR:-}" ]]; then
    cd "${WORKDIR}" && pwd
    return
  fi
  # Prefer data directory next to scripts/ when fastq/bam live there.
  if [[ -d "${SCRIPT_DIR}/../../validate-detections" ]]; then
    cd "${SCRIPT_DIR}/../../validate-detections" && pwd
    return
  fi
  if [[ -d "${SCRIPT_DIR}/../samples" ]]; then
    cd "${SCRIPT_DIR}/.." && pwd
    return
  fi
  cd "${SCRIPT_DIR}/.." && pwd
}

resolve_repo_root() {
  local d="${1}"
  while [[ "${d}" != "/" ]]; do
    if [[ -f "${d}/Assets/NC_045512_Hu-1.fasta" ]]; then
      echo "${d}"
      return 0
    fi
    d="$(dirname "${d}")"
  done
  return 1
}

WORKDIR="$(resolve_workdir)"
SAMPLES_DIR="${WORKDIR}/samples"
FASTQ_DIR="${WORKDIR}/fastq"
BAM_DIR="${WORKDIR}/bam"
FASTERQ_TMP="${WORKDIR}/tmp/fasterq"
FAILED_LOG="${WORKDIR}/failed_accessions.log"

if ! REPO_ROOT="$(resolve_repo_root "${WORKDIR}")"; then
  echo "error: could not find freyja-global repo root (Assets/NC_045512_Hu-1.fasta)" >&2
  exit 1
fi
REF_FASTA="${REF_FASTA:-${REPO_ROOT}/Assets/NC_045512_Hu-1.fasta}"

# Sample CSVs may live next to the script when not copied into WORKDIR.
if [[ ! -d "${SAMPLES_DIR}" && -d "${SCRIPT_DIR}/../samples" ]]; then
  SAMPLES_DIR="${SCRIPT_DIR}/../samples"
fi

# Per-sample failures are handled in the main loop; do not abort the full run.
declare -a FAILED_ACCESSIONS=()
declare -a SKIPPED_ACCESSIONS=()
declare -a OK_ACCESSIONS=()

activate_env() {
  if [[ -n "${CONDA_PREFIX:-}" && "$(basename "${CONDA_PREFIX}")" == "${ENV_NAME}" ]]; then
    return 0
  fi
  local mamba_bin=""
  if command -v micromamba &>/dev/null; then
    mamba_bin="$(command -v micromamba)"
  elif [[ -x "${HOME}/micromamba/bin/micromamba" ]]; then
    mamba_bin="${HOME}/micromamba/bin/micromamba"
  elif [[ -x "${HOME}/.local/bin/micromamba" ]]; then
    mamba_bin="${HOME}/.local/bin/micromamba"
  fi
  if [[ -z "${mamba_bin}" ]]; then
    echo "error: micromamba not found; run: micromamba activate ${ENV_NAME}" >&2
    exit 1
  fi
  # shellcheck source=/dev/null
  eval "$("${mamba_bin}" shell hook -s bash)"
  micromamba activate "${ENV_NAME}"
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [options] [sample_csv ...]

Download FASTQ from SRA for SRR accessions in sample CSVs (column 'accession'),
then align reads to ${REF_FASTA} with minimap2 (-ax sr) into ${BAM_DIR}.

Options:
  -n, --dry-run     Print planned commands without running them
  -h, --help        Show this help

Environment:
  WORKDIR           Analysis root (default: auto-detect from script location)
  ENV_NAME          Micromamba env (default: freyja-sc2)
  THREADS           CPU threads (default: nproc)
  REF_FASTA         Reference FASTA path

Samples with an existing ${BAM_DIR}/<accession>.sorted.bam and .bai are skipped.
Failures are logged to ${FAILED_LOG}; the script continues with remaining samples.

If no CSV paths are given, all files in ${SAMPLES_DIR}/*.csv are used.
EOF
}

DRY_RUN=0
CSV_FILES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    -*) echo "error: unknown option: $1" >&2; usage >&2; exit 1 ;;
    *) CSV_FILES+=("$1"); shift ;;
  esac
done

if [[ ${#CSV_FILES[@]} -eq 0 ]]; then
  mapfile -t CSV_FILES < <(find "${SAMPLES_DIR}" -maxdepth 1 -name '*.csv' -print | sort)
fi

if [[ ${#CSV_FILES[@]} -eq 0 ]]; then
  echo "error: no sample CSV files found under ${SAMPLES_DIR}" >&2
  exit 1
fi

if [[ ! -f "${REF_FASTA}" ]]; then
  echo "error: reference not found: ${REF_FASTA}" >&2
  exit 1
fi

mkdir -p "${FASTQ_DIR}" "${BAM_DIR}" "${FASTERQ_TMP}"

run() {
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    printf '[dry-run]'; printf ' %q' "$@"; printf '\n'
    return 0
  fi
  "$@"
}

log_failure() {
  local acc="$1" step="$2" msg="$3"
  local ts
  ts="$(date -Iseconds)"
  printf '%s\t%s\t%s\n' "${ts}" "${acc}" "${step}: ${msg}" >> "${FAILED_LOG}"
  echo "[${acc}] FAILED (${step}): ${msg}" >&2
}

collect_accessions() {
  local f
  for f in "${CSV_FILES[@]}"; do
    if [[ ! -f "${f}" ]]; then
      echo "error: file not found: ${f}" >&2
      exit 1
    fi
    awk -F',' '
      NR == 1 {
        col = 0
        for (i = 1; i <= NF; i++) {
          gsub(/^[ \t\r]+|[ \t\r]+$/, "", $i)
          if ($i == "accession") col = i
        }
        if (col == 0) {
          print "error: no accession column in " FILENAME > "/dev/stderr"
          exit 1
        }
        next
      }
      {
        gsub(/^[ \t\r]+|[ \t\r]+$/, "", $col)
        if ($col ~ /^SRR[0-9]+$/) print $col
      }
    ' "${f}"
  done | sort -u
}

mapfile -t ACCESSIONS < <(collect_accessions)

if [[ ${#ACCESSIONS[@]} -eq 0 ]]; then
  echo "error: no SRR accessions found in sample CSVs" >&2
  exit 1
fi

echo "Using env: ${ENV_NAME}"
echo "WORKDIR:   ${WORKDIR}"
echo "Samples:   ${#ACCESSIONS[@]} unique accessions from ${#CSV_FILES[@]} CSV file(s)"
echo "FASTQ dir: ${FASTQ_DIR}"
echo "BAM dir:   ${BAM_DIR}"
echo "SRA temp:  ${FASTERQ_TMP}"
echo "Reference: ${REF_FASTA}"
echo "Fail log:  ${FAILED_LOG}"

activate_env

for cmd in prefetch fasterq-dump minimap2 samtools; do
  if ! command -v "${cmd}" &>/dev/null; then
    echo "error: ${cmd} not on PATH after activating ${ENV_NAME}" >&2
    exit 1
  fi
done

sample_complete() {
  local acc="$1"
  local bam="${BAM_DIR}/${acc}.sorted.bam"
  local bai="${bam}.bai"
  [[ -s "${bam}" && -f "${bai}" ]]
}

fastq_ready() {
  local acc="$1"
  [[ -f "${FASTQ_DIR}/${acc}_1.fastq.gz" && -f "${FASTQ_DIR}/${acc}_2.fastq.gz" ]] \
    || [[ -s "${FASTQ_DIR}/${acc}.fastq.gz" ]] \
    || [[ -f "${FASTQ_DIR}/${acc}_1.fastq" && -f "${FASTQ_DIR}/${acc}_2.fastq" ]] \
    || [[ -s "${FASTQ_DIR}/${acc}.fastq" ]]
}

cleanup_partial_fastq() {
  local acc="$1" f
  shopt -s nullglob
  for f in "${FASTQ_DIR}/${acc}"*.fastq "${FASTQ_DIR}/${acc}"*.fastq.gz; do
    [[ -e "${f}" ]] || continue
    if [[ "${DRY_RUN}" -eq 0 ]]; then
      rm -f "${f}"
    else
      echo "[dry-run] rm -f ${f}"
    fi
  done
  shopt -u nullglob
}

cleanup_fasterq_staging() {
  local d
  shopt -s nullglob
  for d in "${WORKDIR}"/fasterq.tmp.*; do
    if [[ -d "${d}" && "${DRY_RUN}" -eq 0 ]]; then
      rm -rf "${d}"
    fi
  done
  shopt -u nullglob
}

gzip_fastqs() {
  local acc="$1" f
  shopt -s nullglob
  for f in "${FASTQ_DIR}/${acc}"*.fastq; do
    if [[ -f "${f}" && ! -f "${f}.gz" ]]; then
      run gzip -f "${f}" || return 1
    fi
  done
  shopt -u nullglob
  return 0
}

fetch_fastq() {
  local acc="$1"
  if fastq_ready "${acc}"; then
    echo "[${acc}] FASTQ already present, skipping download"
    gzip_fastqs "${acc}" || return 1
    return 0
  fi

  cleanup_partial_fastq "${acc}"
  cleanup_fasterq_staging

  echo "[${acc}] prefetch"
  if ! run prefetch "${acc}"; then
    log_failure "${acc}" "prefetch" "prefetch exited with status $?"
    return 1
  fi

  echo "[${acc}] fasterq-dump"
  if ! run fasterq-dump "${acc}" \
    -O "${FASTQ_DIR}" \
    --split-files \
    -e "${THREADS}" \
    -t "${FASTERQ_TMP}" \
    -p; then
    log_failure "${acc}" "fasterq-dump" "fasterq-dump exited with status $?"
    cleanup_partial_fastq "${acc}"
    cleanup_fasterq_staging
    return 1
  fi

  if ! gzip_fastqs "${acc}"; then
    log_failure "${acc}" "gzip" "failed to compress FASTQ files"
    return 1
  fi

  cleanup_fasterq_staging

  if [[ "${DRY_RUN}" -eq 0 ]] && ! fastq_ready "${acc}"; then
    log_failure "${acc}" "fasterq-dump" "no FASTQ output found in ${FASTQ_DIR}"
    cleanup_partial_fastq "${acc}"
    return 1
  fi
  return 0
}

resolve_reads() {
  local acc="$1"
  READS=()
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    READS=("${FASTQ_DIR}/${acc}_1.fastq.gz" "${FASTQ_DIR}/${acc}_2.fastq.gz")
    return 0
  fi
  if [[ -f "${FASTQ_DIR}/${acc}_1.fastq.gz" && -f "${FASTQ_DIR}/${acc}_2.fastq.gz" ]]; then
    READS=("${FASTQ_DIR}/${acc}_1.fastq.gz" "${FASTQ_DIR}/${acc}_2.fastq.gz")
  elif [[ -f "${FASTQ_DIR}/${acc}_1.fastq" && -f "${FASTQ_DIR}/${acc}_2.fastq" ]]; then
    READS=("${FASTQ_DIR}/${acc}_1.fastq" "${FASTQ_DIR}/${acc}_2.fastq")
  elif [[ -f "${FASTQ_DIR}/${acc}.fastq.gz" ]]; then
    READS=("${FASTQ_DIR}/${acc}.fastq.gz")
  elif [[ -f "${FASTQ_DIR}/${acc}.fastq" ]]; then
    READS=("${FASTQ_DIR}/${acc}.fastq")
  else
    return 1
  fi
  return 0
}

align_sample() {
  local acc="$1"
  local bam="${BAM_DIR}/${acc}.sorted.bam"
  local bai="${bam}.bai"

  if sample_complete "${acc}"; then
    echo "[${acc}] already processed (sorted BAM + index present), skipping"
    return 2
  fi

  if ! fetch_fastq "${acc}"; then
    return 1
  fi

  if ! resolve_reads "${acc}"; then
    log_failure "${acc}" "align" "cannot determine read layout in ${FASTQ_DIR}"
    return 1
  fi

  echo "[${acc}] minimap2 align (${#READS[@]} file(s))"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    run minimap2 -ax sr -t "${THREADS}" "${REF_FASTA}" "${READS[@]}"
    run samtools sort -@ "${THREADS}" -o "${bam}"
    run samtools index "${bam}"
    return 0
  fi

  local sort_tmp="${BAM_DIR}/.${acc}.sorted.bam.tmp"
  if ! minimap2 -ax sr -t "${THREADS}" "${REF_FASTA}" "${READS[@]}" \
    | samtools sort -@ "${THREADS}" -o "${sort_tmp}"; then
    log_failure "${acc}" "align" "minimap2 or samtools sort failed"
    rm -f "${sort_tmp}"
    return 1
  fi

  if ! samtools index "${sort_tmp}"; then
    log_failure "${acc}" "align" "samtools index failed"
    rm -f "${sort_tmp}" "${sort_tmp}.bai"
    return 1
  fi

  mv -f "${sort_tmp}" "${bam}"
  mv -f "${sort_tmp}.bai" "${bai}"
  return 0
}

for acc in "${ACCESSIONS[@]}"; do
  set +e
  align_sample "${acc}"
  rc=$?
  set -e
  case "${rc}" in
    0) OK_ACCESSIONS+=("${acc}") ;;
    2) SKIPPED_ACCESSIONS+=("${acc}") ;;
    *)
      FAILED_ACCESSIONS+=("${acc}")
      ;;
  esac
done

echo ""
echo "========== Summary =========="
echo "Total accessions: ${#ACCESSIONS[@]}"
echo "Completed now:    ${#OK_ACCESSIONS[@]}"
echo "Skipped (done):   ${#SKIPPED_ACCESSIONS[@]}"
echo "Failed:           ${#FAILED_ACCESSIONS[@]}"

if [[ ${#FAILED_ACCESSIONS[@]} -gt 0 ]]; then
  echo "Failed accessions: ${FAILED_ACCESSIONS[*]}"
  echo "Details: ${FAILED_LOG}"
  exit 1
fi

echo "All samples succeeded or were already processed."
