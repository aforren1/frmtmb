#!/usr/bin/env bash
# The rest of the low-noise read, apart from brms's own program, whose chains at
# depth 14 run long enough to hold the others back.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
L=C:/Users/adf44/source/r/wt-ncgp-lib
s=1e-2
"$RS" dev/ncgp-sampler.R $s 30 1000 0.99 > dev/ncgp-out/post/sampler-$s.txt 2>&1 &
for c in pois12 pois16 binom1e7; do
  "$RS" dev/ncgp-nongauss2.R $c $L 1500 > dev/ncgp-out/post/ngz-lane-$c.txt 2>&1
done
wait
echo POST5-DONE
