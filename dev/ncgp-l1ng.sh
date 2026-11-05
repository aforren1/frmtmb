#!/usr/bin/env bash
# dev/ncgp-l1ng.R on each large-mean case on the final build, and bernoulli.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
O=dev/ncgp-out/l1ng.txt
: > "$O"
for c in pois4 pois12 pois16 binom1e5 binom1e7 bern; do
  "$RS" dev/ncgp-l1ng.R "$c" 2>&1 | grep -a -E 'NCGPLN|^Error' >> "$O"
done
echo L1NG-DONE
