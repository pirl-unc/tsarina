#!/usr/bin/env bash
# Run the tsarina test suite with a memory-aware pytest-xdist worker
# count: ``-n auto`` spawns one worker per core, but each worker
# loads its own DataFrame state and can OOM a 32 GB Mac when other
# pytest suites are running concurrently. We pick the smaller of
# (cores, available_RAM / per_worker_gb).

set -e

# Per-worker memory budget. A soft heuristic, not a hard cap.
readonly PER_WORKER_GB=2.0

# macOS available-RAM heuristic: free + inactive + speculative pages
# are reclaimable on demand, so they count as headroom.
macos_available_bytes() {
    local page_size pages
    page_size=$(sysctl -n hw.pagesize)
    pages=$(vm_stat | awk '
        /Pages free/        { gsub(/\./, "", $3); free     = $3 }
        /Pages inactive/    { gsub(/\./, "", $3); inactive = $3 }
        /Pages speculative/ { gsub(/\./, "", $3); spec     = $3 }
        END                 { print free + inactive + spec }
    ')
    echo $(( pages * page_size ))
}

# Pick a pytest -n value that respects both CPU and available RAM.
pytest_workers() {
    local cpus avail_bytes
    cpus=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.logicalcpu)
    if [[ "$(uname)" == "Darwin" ]]; then
        avail_bytes=$(macos_available_bytes)
    else
        avail_bytes=$(awk '/MemAvailable/ { print $2 * 1024 }' /proc/meminfo)
    fi
    awk -v cpus="$cpus" -v bytes="$avail_bytes" -v budget="$PER_WORKER_GB" '
        BEGIN {
            by_memory = int(bytes / 1024^3 / budget)
            n = (cpus < by_memory ? cpus : by_memory)
            print (n < 1 ? 1 : n)
        }
    '
}

# pytest-xdist is optional — fall back to serial pytest if missing.
xdist_flag=()
if python -c "import xdist" 2>/dev/null; then
    workers=$(pytest_workers)
    xdist_flag=(-n "$workers")
    echo "Running pytest with -n ${workers} (per-worker budget ≈ ${PER_WORKER_GB} GB)"
fi

exec pytest "${xdist_flag[@]}" --cov=tsarina/ --cov-report=term-missing tests "$@"
