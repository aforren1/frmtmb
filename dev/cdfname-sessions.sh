#!/usr/bin/env bash
# Each suite in ONE session, alphabetical order, on the lane's build:
# core in one slot, the eight extensions one after another in the other,
# so at most two R processes run.
#   bash dev/cdfname-sessions.sh
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-cdfname
cd $WT || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export NOT_CRAN=true FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-cdfname-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=${OUT:-dev/cdfname-sessions}
mkdir -p $OUT
( "$RS" dev/release/run-session.R frmtmb forward > $OUT/frmtmb.txt 2>&1 ) &
(
  for e in ${EXTS:-eam survival coupling latent learn ode sample spatial}; do
    "$RS" dev/release/run-session.R frmtmb.$e forward \
      > $OUT/frmtmb.$e.txt 2>&1
  done
) &
wait
for f in $OUT/*.txt; do
  echo "$(basename $f .txt): $(grep -a '^SESSION' $f | tail -1)"
  grep -a '^LEAK' $f | grep -v 'extension session state' | sed 's/^/  /'
done
echo "SESSIONS END"
