#!/usr/bin/env bash
# What dev/ncgp-post.sh had left when the build moved, then the large-mean designs
# with chains started at frm()'s z.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
L=C:/Users/adf44/source/r/wt-ncgp-lib
"$RS" dev/ncgp-sampler.R 3e-3 20 3000 0.99 > dev/ncgp-out/post/sampler-3e-3.txt 2>&1
s=1e-3
"$RS" dev/ncgp-post.R $s 20 base > dev/ncgp-out/post/base-$s.txt 2>&1
"$RS" dev/ncgp-post.R $s 20 $L 3000 0.99 > dev/ncgp-out/post/lane-$s.txt 2>&1
"$RS" dev/ncgp-sampler.R $s 20 3000 0.99 > dev/ncgp-out/post/sampler-$s.txt 2>&1
for c in pois12 pois16 binom1e7; do
  "$RS" dev/ncgp-nongauss2.R $c $L 1500 > dev/ncgp-out/post/ngz-lane-$c.txt 2>&1
done
echo POST3-DONE
