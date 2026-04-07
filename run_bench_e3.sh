#!/usr/bin/env bash
# E3: HomSort + Sync benchmark — n={64,256,1024,4096,8192}, c=8 columns
set -e
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -x "${DIR}/build/bin/bench_e3" ]]; then
    BIN="${DIR}/build/bin/bench_e3"
elif [[ -x "${DIR}/bin/bench_e3" ]]; then
    BIN="${DIR}/bin/bench_e3"
else
    echo "Binary not found. Run: cmake . && make bench_e3 -j4" >&2
    exit 1
fi
TS=$(date +%Y%m%d_%H%M%S)
LOG="${DIR}/bench_e3_${TS}.log"

echo "===== E3: HomSort+Sync sweep -> $LOG ====="
stdbuf -oL "$BIN" 2>&1 | tee "$LOG"
echo ""
echo "Done. Log: $LOG"
echo "CSV:  engorgio_e3_results.csv"
