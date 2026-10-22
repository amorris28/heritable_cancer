#!/usr/bin/env bash
# Absolute-risk calibration of the frozen PHS601 weights against registry incidence.
# Usage: scripts/10_calibration_frozen.sh <p212|p2521> [--force]
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

OUTPUT="data/public/${CLUSTER}_calibration_frozen.tsv"
should_run "$OUTPUT" || exit 0

run_rscript src/calibration-frozen.R "$PHENOTYPES" "$CLUSTER"
require_files "${CLUSTER}_calibration_frozen.tsv"
mv "${CLUSTER}_calibration_frozen.tsv" "$OUTPUT"
mv "${CLUSTER}_calibration_frozen.png" "data/public/${CLUSTER}_calibration_frozen.png"

echo "Done. Output: $OUTPUT"
