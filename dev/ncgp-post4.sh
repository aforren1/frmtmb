#!/usr/bin/env bash
# The low-noise posterior read, with converged chains, and brms's own program at
# the lowest noise it mixes at; then the large-mean designs from frm()'s z.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
L=C:/Users/adf44/source/r/wt-ncgp-lib
"$RS" dev/ncgp-marg.R base > dev/ncgp-out/post/marg-base.txt 2>&1
"$RS" dev/ncgp-marg.R $L > dev/ncgp-out/post/marg-lane.txt 2>&1
s=1e-2
"$RS" dev/ncgp-post.R $s 30 base > dev/ncgp-out/post/base-$s.txt 2>&1
"$RS" dev/ncgp-post.R $s 30 $L 2000 0.99 > dev/ncgp-out/post/lane-$s.txt 2>&1
"$RS" dev/ncgp-sampler.R $s 30 2000 0.99 > dev/ncgp-out/post/sampler-$s.txt 2>&1
for c in pois12 pois16 binom1e7; do
  "$RS" dev/ncgp-nongauss2.R $c $L 1500 > dev/ncgp-out/post/ngz-lane-$c.txt 2>&1
done
echo POST4-DONE
