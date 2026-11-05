#!/usr/bin/env bash
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
mkdir -p dev/ncgp-out/permsamp
for s in 101 102 103 104; do
  for a in sorted perm; do
    "$RS" dev/ncgp-permsamp.R $a $s 2000 12 15 > dev/ncgp-out/permsamp/$a-15-$s.txt 2>&1 &
  done
  while [ "$(jobs -rp | wc -l)" -ge 6 ]; do sleep 2; done
done
wait
echo PERMSAMP-DONE
