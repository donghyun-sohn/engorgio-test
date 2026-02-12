#!/usr/bin/env bash
# Q1, Q6, Q12 benchmark: 2K, 4K, 8K, 16K rows each. Output tee'd to timestamped log.
set -e
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${DIR}/build/bin/relational_query_test"
TS=$(date +%Y%m%d_%H%M%S)
LOG="${DIR}/relational_bench_${TS}.log"

if [[ ! -x "$BIN" ]]; then
    echo "Binary not found: $BIN (run: cd build && make relational_query_test -j4)" >&2
    exit 1
fi

echo "===== Relational benchmark (Q1, Q6, Q12 @ 2K,4K,8K,16K) -> $LOG ====="
for N in 2048 4096 8192 16384; do
    echo ""
    echo "========================================  $N rows  ========================================"
    echo "--- Q1 ($N rows) ---"
    "$BIN" --q1 "$N" 2>&1
    echo ""
    echo "--- Q6 ($N rows) ---"
    "$BIN" --q6 "$N" 2>&1
    echo ""
    echo "--- Q12 ($N rows) ---"
    "$BIN" --q12 "$N" 2>&1
done | tee "$LOG"

echo ""
echo "Done. Log: $LOG"
