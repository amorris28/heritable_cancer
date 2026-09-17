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
