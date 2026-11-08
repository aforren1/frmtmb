#!/bin/bash
# Per-family memory grid, three processes at a time.
#   bash dev/nestquad-mem-grid.sh <outfile> "<family> <rows>" ...
cd C:/Users/adf44/source/r/frmtmb-wt-nestquad
export PATH="$PATH:/c/rtools45/usr/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
out=$1; shift
i=0
for cfg in "$@"; do
  ("$R" dev/nestquad-mem.R $cfg 2>&1 | grep "^MEM" >> "$out") &
  i=$((i + 1))
  if [ $((i % 3)) -eq 0 ]; then wait; fi
done
wait
