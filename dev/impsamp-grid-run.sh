#!/usr/bin/env bash
# Run one arm of the impsamp grid: bash dev/impsamp-grid-run.sh <arm> [lib]
# One R process per chunk, at most 8 at once. Run from the worktree root.
arm="$1"
lib="${2:-}"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
n=8
for c in $(seq 0 $((n - 1))); do
  "$R" dev/impsamp-grid.R "$arm" "$c" "$n" "$lib" \
    > "dev/impsamp-out/grid-$arm-$c.log" 2>&1 &
done
wait
echo "grid $arm finished"
