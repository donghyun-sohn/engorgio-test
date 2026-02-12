#!/usr/bin/env bash
# HomSort benchmark: 4, 16, 64, 256, 1024, 4096, 16384 records. Output tee'd to timestamped log.
set -e
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${DIR}/build/bin/sort_test"
TS=$(date +%Y%m%d_%H%M%S)
LOG="${DIR}/sort_test_${TS}.log"

if [[ ! -x "$BIN" ]]; then
    echo "Binary not found: $BIN (run: cd build && make sort_test -j4)" >&2
    exit 1
fi

echo "===== Sort test (4,16,64,256,1024,4096,16384) -> $LOG ====="
for N in 4 16 64 256 1024 4096 16384; do
    echo ""
    echo "========== $N records =========="
    "$BIN" "$N" 2>&1
done | tee "$LOG"
echo ""
echo "Done. Log: $LOG"
