#!/usr/bin/env bash
# Build the p212 KING retained-ID lists (0deg/2deg/3deg) from the externally
# computed KING text kinship matrix. Light enough to run on a login node.
# Usage: scripts/08_p212_relatedness.sh [--force]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"
load_cohort_config p212

[[ "${1:-}" == "--force" ]] && FORCE=1
cd "$REPO_ROOT"

KING_TEXT="${KING_TEXT:-/ess/p212/cluster/king/kinship_results25.king}"
KING_IDS="${KING_IDS:-${KING_TEXT}.id}"
BFILE="${KING_BFILE:-$BED_FILE}"
PLINK2="../bin/plink2"
OUT="data/interim/relatedness"

outputs=("$OUT"/p212_{0,2,3}deg.king.cutoff.in.id)
should_run "${outputs[@]}" || exit 0

mkdir -p "$OUT"
python3 "${SCRIPT_DIR}/lib/king_text_to_bin.py" "$KING_TEXT" "$KING_IDS" "$OUT/p212"

# 0deg: duplicates/MZ; 2deg: second degree; 3deg: third degree
for spec in 0deg:0.35333 2deg:0.0883883 3deg:0.0441942; do
  "$PLINK2" --bfile "$BFILE" \
    --king-cutoff "$OUT/p212" "${spec#*:}" \
    --out "$OUT/p212_${spec%%:*}"
done
