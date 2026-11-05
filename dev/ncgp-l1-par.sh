#!/usr/bin/env bash
# dev/ncgp-l1.R over every spelling, on brms's own program and then the attribution
# arm, 6 at a time, one log per job in dev/ncgp-out/l1/.
set -u
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
mkdir -p dev/ncgp-out/l1
for arm in "" attr; do
  for i in $(seq 1 20); do
    while [ "$(jobs -rp | wc -l)" -ge 6 ]; do sleep 1; done
    ( "$RS" dev/ncgp-l1.R "$i" $arm > "dev/ncgp-out/l1/l1-$i-${arm:-brms}.txt" 2>&1 ) &
  done
done
wait
echo L1-ALL-DONE
