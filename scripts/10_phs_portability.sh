#!/usr/bin/env bash
# Portability of the published PHS601 weights: Cox slope on the unstandardized PHS.
# Usage: scripts/10_phs_portability.sh <p212|p2521> [--force]
# Needs phenotypes.tsv with a phs_raw column (rerun scripts/02_preprocess_metadata.sh).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

COHORT="${1:?Usage: $0 <p212|p2521> [--force]}"
[[ "${2:-}" == "--force" ]] && FORCE=1
load_cohort_config "$COHORT"
cd "$REPO_ROOT"

PHENOTYPES=data/processed/phenotypes.tsv
[[ -f "$PHENOTYPES" ]] || { echo "Missing $PHENOTYPES — run scripts/02_preprocess_metadata.sh $CLUSTER first" >&2; exit 1; }

OUTPUT="data/public/${CLUSTER}_phs_portability.tsv"
should_run "$OUTPUT" || exit 0

run_rscript src/phs-portability.R "$PHENOTYPES" "$CLUSTER"
require_files "${CLUSTER}_phs_portability.tsv"
mv "${CLUSTER}_phs_portability.tsv" "$OUTPUT"

echo "Done. Output: $OUTPUT"
