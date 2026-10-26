#!/usr/bin/env bash
# Shared helpers for scripts/*.sh — no framework, just idempotent file checks.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FORCE="${FORCE:-0}"
BOOTSTRAP_ITERS="${BOOTSTRAP_ITERS:-999}"

# load_cohort_config <p212|p2521>: sources config/<cohort>.env into the current shell
load_cohort_config() {
  local cohort="$1"
  local env_file="${REPO_ROOT}/config/${cohort}.env"
  [[ -f "$env_file" ]] || { echo "No such cohort config: $env_file" >&2; exit 1; }
  # shellcheck disable=SC1090
  source "$env_file"
}

# king_keep_file <p212|p2521> <0deg|2deg|3deg>: path of the KING retained-ID list
king_keep_file() {
  if [[ "$1" == "p2521" ]]; then
    echo "/ess/p2521/cluster/users/espehage/BiobankNorge_AD/results/husk/all_chr.cutoff_${2}.king.cutoff.in.id"
  else
    echo "data/interim/relatedness/p212_${2}.king.cutoff.in.id"
  fi
}

# should_run <output_file> [<output_file> ...]
# Prints a skip message and returns false if every given output already exists
# (unless FORCE=1). Long-running/one-off steps use this to avoid recomputation;
# this replaces Nextflow's storeDir/-resume caching with plain file checks.
should_run() {
  if [[ "$FORCE" == "1" ]]; then
    return 0
  fi
  local f
  for f in "$@"; do
    [[ -e "$f" ]] || return 0
  done
  echo "  [skip] already present: $* (rerun with FORCE=1 to regenerate)"
  return 1
}

# run_rscript <script.R> [args...]
# Runs an Rscript and fails loudly (rather than silently continuing) if it
# exits non-zero — belt-and-suspenders on top of `set -e`.
run_rscript() {
  if ! Rscript "$@"; then
    echo "ERROR: Rscript $* exited with a non-zero status; aborting before moving any output." >&2
    exit 1
  fi
}

# run_rscript_docker <image> <script.R> [args...]
# Same as run_rscript, but inside the given docker image (repo mounted at the
# same path so relative sourcing/reads/writes in the R scripts still work).
run_rscript_docker() {
  local image="$1"; shift
  if ! docker run --rm --platform linux/amd64 \
        -v "${REPO_ROOT}:${REPO_ROOT}" -w "$REPO_ROOT" \
        "$image" Rscript "$@"; then
    echo "ERROR: Rscript $* (docker image: $image) exited with a non-zero status; aborting before moving any output." >&2
    exit 1
  fi
}

# require_files <file> [<file> ...]
# Fails loudly if any expected R output is missing, e.g. because the script
# exited 0 but errored/warned partway through without writing that file.
require_files() {
  local missing=() f
  for f in "$@"; do
    [[ -f "$f" ]] || missing+=("$f")
  done
  if (( ${#missing[@]} > 0 )); then
    echo "ERROR: expected output(s) not found (check R output above for errors/warnings): ${missing[*]}" >&2
    exit 1
  fi
}
