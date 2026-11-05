#!/usr/bin/env bash
# After dev/ncgp-post.sh: the large-mean designs with chains started at frm()'s z.
until grep -q POST-DONE dev/ncgp-out/post.log; do sleep 30; done
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
L=C:/Users/adf44/source/r/wt-ncgp-lib
for c in pois12 pois16 binom1e7; do
  "$RS" dev/ncgp-nongauss2.R $c $L 1500 > dev/ncgp-out/post/ngz-lane-$c.txt 2>&1
done
echo POST2-DONE
