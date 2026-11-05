#!/usr/bin/env bash
# All of dev/ncgp-l1.R again on the final build, 6 at a time, then the mixture by
# hand, and the equivalence pairs that hold an exact gp() against the baseline.
set -u
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
rm -rf dev/ncgp-out/l1; mkdir -p dev/ncgp-out/l1
for arm in "" attr; do
  for i in $(seq 1 19); do
    while [ "$(jobs -rp | wc -l)" -ge 6 ]; do sleep 1; done
    ( "$RS" dev/ncgp-l1.R "$i" $arm > "dev/ncgp-out/l1/l1-$i-${arm:-brms}.txt" 2>&1 ) &
  done
done
wait
"$RS" dev/ncgp-mix.R > dev/ncgp-out/mix.txt 2>&1
echo L1-ALL-DONE
mkdir -p dev/ncgp-equiv-out
cp -n ../frmtmb-wt-r16base/dev/r16base-equiv-out/*.tsv dev/ncgp-equiv-out/
rm -f dev/ncgp-equiv-out/p022.tsv dev/ncgp-equiv-out/p023.tsv dev/ncgp-equiv-out/p094.tsv \
      dev/ncgp-equiv-out/p109.tsv dev/ncgp-equiv-out/p136.tsv
FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-ncgp-lib;C:/Users/adf44/source/r/rellib-r16" \
  bash dev/equiv-run.sh 5 dev/ncgp-equiv-out > dev/ncgp-out/equiv-run.txt 2>&1
echo EQUIV-DONE
