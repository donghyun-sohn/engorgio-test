#!/usr/bin/env bash
# Run E2 (Filter+Agg) then E3 (HomSort+Sync) in order
set -e
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RECORDS="${1:-60175}"

echo "============================================"
echo "Engorgio Full Benchmark: E2 then E3"
echo "Started: $(date)"
echo "============================================"
echo ""

# --- E2 ---
"${DIR}/run_bench_e2.sh" "$RECORDS"

echo ""
echo "============================================"
echo "E2 finished: $(date)"
echo "Starting E3..."
echo "============================================"
echo ""

# --- E3 ---
"${DIR}/run_bench_e3.sh"

echo ""
echo "============================================"
echo "All benchmarks complete: $(date)"
echo "============================================"
