#!/usr/bin/env bash
# Evidence against brms's posterior, one 4-chain run at a time. Logs in
# dev/ncgp-out/post/.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
L=C:/Users/adf44/source/r/wt-ncgp-lib
mkdir -p dev/ncgp-out/post
for c in pois4 pois12 pois16 binom1e5 binom1e7 bern; do
  "$RS" dev/ncgp-nongauss.R $c base > dev/ncgp-out/post/ng-base-$c.txt 2>&1
  "$RS" dev/ncgp-nongauss.R $c $L 1500 > dev/ncgp-out/post/ng-lane-$c.txt 2>&1
done
for s in 3e-3 1e-3; do
  "$RS" dev/ncgp-post.R $s 20 base > dev/ncgp-out/post/base-$s.txt 2>&1
  "$RS" dev/ncgp-post.R $s 20 $L 3000 0.99 > dev/ncgp-out/post/lane-$s.txt 2>&1
  "$RS" dev/ncgp-sampler.R $s 20 3000 0.99 > dev/ncgp-out/post/sampler-$s.txt 2>&1
done
echo POST-DONE
