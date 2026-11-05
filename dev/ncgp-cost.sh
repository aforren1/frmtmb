#!/usr/bin/env bash
# Wall clock alone, with nothing else of the lane's running: the reviewer's script with
# 3 repeats, base and lane, one process at a time, then the warm start's own parts.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
O=dev/ncgp-out/cost.txt
: > "$O"
for rep in 1 2; do
  for n in 60 200 500 1000; do
    for arm in base lane; do
      "$RS" dev/review-ncgp-cost.R $arm $n 3 2>&1 | grep -a -E '^(lib|RCO)' >> "$O"
    done
  done
done
"$RS" dev/ncgp-warmcost.R > dev/ncgp-out/warmcost.txt 2>&1
echo COST-DONE
