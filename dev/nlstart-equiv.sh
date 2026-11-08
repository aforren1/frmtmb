#!/usr/bin/env bash
# Remeasure the equivalence survey's pairs that hold an nl body, on the
# lane's build, three at a time. Output dev/nlstart-equiv-out/<id>.tsv.
cd /c/Users/adf44/source/r/frmtmb-wt-nlstart || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-nlstart-lib;C:/Users/adf44/source/r/rellib-r18"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlstart/dev/stan-cache
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
ids="p001 p014 p031 p032 p042 p047 p049 p119 p127 p136 p143 p158 p166 p171"
mkdir -p dev/nlstart-equiv-out
running=0
for id in $ids; do
  "$RS" dev/equiv-measure.R "$id" dev/nlstart-equiv-out --force \
    > "dev/nlstart-equiv-out/$id.log" 2>&1 &
  running=$((running + 1))
  if [ "$running" -ge 3 ]; then
    wait -n
    running=$((running - 1))
  fi
done
wait
echo "EQUIV DONE"
