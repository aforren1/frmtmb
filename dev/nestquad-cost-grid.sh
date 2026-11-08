#!/bin/bash
# Memory and time grid for the nested tape, three processes at a time.
#   bash dev/nestquad-cost-grid.sh <outfile> <config>...
# Each config is "family n_outer inner_per_outer rows_per_inner ko ki".
cd C:/Users/adf44/source/r/frmtmb-wt-nestquad
export PATH="$PATH:/c/rtools45/usr/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
out=$1; shift
i=0
for cfg in "$@"; do
  ("$R" dev/nestquad-cost.R lib $cfg 2>&1 | grep COST >> "$out") &
  i=$((i + 1))
  if [ $((i % 3)) -eq 0 ]; then wait; fi
done
wait
