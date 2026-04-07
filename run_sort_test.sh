#!/usr/bin/env bash
# HomSort benchmark: 6005, 60175, 600572 records. Output tee'd to timestamped log.
set -e
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -x "${DIR}/build/bin/sort_test" ]]; then
    BIN="${DIR}/build/bin/sort_test"
elif [[ -x "${DIR}/bin/sort_test" ]]; then
    BIN="${DIR}/bin/sort_test"
else
    echo "Binary not found. Run: cmake . && make sort_test -j4 (or build/)"
    exit 1
fi
TS=$(date +%Y%m%d_%H%M%S)
LOG="${DIR}/sort_test_${TS}.log"

echo "===== Sort test (6005,60175,600572) -> $LOG ====="
for N in 6005 60175 600572; do
    echo ""
    echo "========== $N records =========="
    "$BIN" "$N" 2>&1
done | tee "$LOG"
echo ""
echo "Done. Log: $LOG"
