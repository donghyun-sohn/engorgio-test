#!/usr/bin/env bash
# E2: Filter+Agg benchmark (Q1, Q6, Q12) — records_num=60175 (SF=0.01)
set -e
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -x "${DIR}/build/bin/bench_e2" ]]; then
    BIN="${DIR}/build/bin/bench_e2"
elif [[ -x "${DIR}/bin/bench_e2" ]]; then
    BIN="${DIR}/bin/bench_e2"
else
    echo "Binary not found. Run: cmake . && make bench_e2 -j4" >&2
    exit 1
fi
RECORDS="${1:-60175}"
TS=$(date +%Y%m%d_%H%M%S)
LOG="${DIR}/bench_e2_${TS}.log"

echo "===== E2: Q1, Q6, Q12 @ ${RECORDS} rows -> $LOG ====="
stdbuf -oL "$BIN" "$RECORDS" 2>&1 | tee "$LOG"
echo ""
echo "Done. Log: $LOG"
