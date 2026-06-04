#!/usr/bin/env bash
# Run covar on sorted/indexed BAMs produced by fetch_and_align_sra.sh.
# Requires micromamba env freyja-sc2 (or covar on PATH after activation).
set -euo pipefail

ENV_NAME="${ENV_NAME:-freyja-sc2}"
THREADS="16"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

resolve_workdir() {
  if [[ -n "${WORKDIR:-}" ]]; then
    cd "${WORKDIR}" && pwd
    return
  fi
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
BAM_DIR="${WORKDIR}/bam"
COVAR_DIR="${WORKDIR}/covar"
FAILED_LOG="${WORKDIR}/failed_covar.log"

if ! REPO_ROOT="$(resolve_repo_root "${WORKDIR}")"; then
  echo "error: could not find freyja-global repo root (Assets/NC_045512_Hu-1.fasta)" >&2
  exit 1
fi

REF_FASTA="${REF_FASTA:-${REPO_ROOT}/Assets/NC_045512_Hu-1.fasta}"
GFF="${GFF:-${REPO_ROOT}/Analyses/Bangalore/scripts/NC_045512_Hu-1.gff}"

declare -a FAILED_SAMPLES=()
declare -a SKIPPED_SAMPLES=()
declare -a OK_SAMPLES=()

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
Usage: $(basename "$0") [options] [bam_file ...]

Run covar on BAM files in ${BAM_DIR} (default) or on the paths given as arguments.

Options:
  -n, --dry-run     Print planned commands without running them
  -f, --force       Re-run covar even if output TSV already exists
  -h, --help        Show this help

Environment:
  WORKDIR           Analysis root (default: auto-detect)
  ENV_NAME          Micromamba env (default: freyja-sc2)
  THREADS           Threads passed to covar (default: nproc)
  REF_FASTA         Reference FASTA (default: Assets/NC_045512_Hu-1.fasta)
  GFF               Annotation GFF3 (default: Analyses/Bangalore/scripts/NC_045512_Hu-1.gff)

Output: ${COVAR_DIR}/<accession>.covar.tsv

covar expects primer-trimmed, sorted, indexed BAMs. Failures are logged to
${FAILED_LOG}; processing continues for remaining samples.
EOF
}

DRY_RUN=0
FORCE=0
BAM_FILES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--dry-run) DRY_RUN=1; shift ;;
    -f|--force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    -*) echo "error: unknown option: $1" >&2; usage >&2; exit 1 ;;
    *) BAM_FILES+=("$1"); shift ;;
  esac
done

if [[ ! -f "${REF_FASTA}" ]]; then
  echo "error: reference not found: ${REF_FASTA}" >&2
  exit 1
fi

if [[ ! -f "${GFF}" ]]; then
  echo "error: GFF not found: ${GFF}" >&2
  exit 1
fi

mkdir -p "${COVAR_DIR}"

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
  mapfile -t BAM_FILES < <(find "${BAM_DIR}" -maxdepth 1 -name '*.sorted.bam' -print | sort)
fi

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
  echo "error: no BAM files found (looked in ${BAM_DIR} for *.sorted.bam)" >&2
  exit 1
fi

echo "Using env: ${ENV_NAME}"
echo "WORKDIR:   ${WORKDIR}"
echo "BAM dir:   ${BAM_DIR}"
echo "Output:    ${COVAR_DIR}"
echo "Reference: ${REF_FASTA}"
echo "GFF:       ${GFF}"
echo "BAMs:      ${#BAM_FILES[@]} file(s)"
echo "Fail log:  ${FAILED_LOG}"

activate_env

if ! command -v covar &>/dev/null; then
  echo "error: covar not on PATH after activating ${ENV_NAME}" >&2
  exit 1
fi

run() {
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    printf '[dry-run]'; printf ' %q' "$@"; printf '\n'
    return 0
  fi
  "$@"
}

log_failure() {
  local sample="$1" step="$2" msg="$3"
  local ts
  ts="$(date -Iseconds)"
  printf '%s\t%s\t%s\n' "${ts}" "${sample}" "${step}: ${msg}" >> "${FAILED_LOG}"
  echo "[${sample}] FAILED (${step}): ${msg}" >&2
}

accession_from_bam() {
  local bam="$1"
  local base
  base="$(basename "${bam}")"
  base="${base%.sorted.bam}"
  base="${base%.bam}"
  echo "${base}"
}

covar_output_path() {
  local acc="$1"
  echo "${COVAR_DIR}/${acc}.covar.tsv"
}

covar_done() {
  local out="$1"
  [[ -s "${out}" ]]
}

run_covar() {
  local bam="$1"
  local acc out bai
  acc="$(accession_from_bam "${bam}")"
  out="$(covar_output_path "${acc}")"
  bai="${bam}.bai"

  if [[ ! -f "${bam}" ]]; then
    log_failure "${acc}" "input" "BAM not found: ${bam}"
    return 1
  fi

  if [[ ! -f "${bai}" ]]; then
    echo "[${acc}] indexing BAM"
    if ! run samtools index "${bam}"; then
      log_failure "${acc}" "index" "samtools index failed"
      return 1
    fi
  fi

  if [[ "${FORCE}" -eq 0 ]] && covar_done "${out}"; then
    echo "[${acc}] covar output present, skipping (${out})"
    return 2
  fi

  echo "[${acc}] covar"
  if ! run covar \
    --input "${bam}" \
    --reference "${REF_FASTA}" \
    --annotation "${GFF}" \
    --output "${out}" \
    --threads "${THREADS}"; then
    log_failure "${acc}" "covar" "covar exited with status $?"
    rm -f "${out}"
    return 1
  fi

  if [[ "${DRY_RUN}" -eq 0 ]] && ! covar_done "${out}"; then
    log_failure "${acc}" "covar" "no output written to ${out}"
    return 1
  fi

  return 0
}

for bam in "${BAM_FILES[@]}"; do
  set +e
  run_covar "${bam}"
  rc=$?
  set -e
  acc="$(accession_from_bam "${bam}")"
  case "${rc}" in
    0) OK_SAMPLES+=("${acc}") ;;
    2) SKIPPED_SAMPLES+=("${acc}") ;;
    *) FAILED_SAMPLES+=("${acc}") ;;
  esac
done

echo ""
echo "========== Summary =========="
echo "BAMs processed:   ${#BAM_FILES[@]}"
echo "Completed now:    ${#OK_SAMPLES[@]}"
echo "Skipped (done):   ${#SKIPPED_SAMPLES[@]}"
echo "Failed:           ${#FAILED_SAMPLES[@]}"

if [[ ${#FAILED_SAMPLES[@]} -gt 0 ]]; then
  echo "Failed samples: ${FAILED_SAMPLES[*]}"
  echo "Details: ${FAILED_LOG}"
  exit 1
fi

echo "All samples succeeded or were already processed."
