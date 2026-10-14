#!/usr/bin/env bash
# Lane illse's SE-check batteries (its dev/illse-sf-run.sh) on a build of
# the release tree: the false-alarm battery (5 groups, 820 fits per
# BLAS) and the 30-fit bitwise set, under the reference BLAS and
# OpenBLAS 0.3.32 with LAPACK (4 threads), as the lane ran them. The
# lane's scripts read and write its own worktree, so they run here as
# copies with the lane library replaced by the arm's library and the
# lane's dev/ by the copies' directory. Arms: "merged" (<lib>) and
# "lane14" (wt-illse-lib, the lane's final install 14; its own battery
# logs, sf-g2, are from install 13). The base arm (rellib-r8, frmtmb
# 0.70.0) is the lane's sf-g2 run, copied in; rellib-r8 has not changed.
#
#   bash dev/rel071-rr-illse.sh <lib> <tag>
set -u
LIB="$1"; TAG="$2"
REL=/c/Users/adf44/source/r/frmtmb-wt-release
LANE=/c/Users/adf44/source/r/frmtmb-wt-illse
W=$REL/dev/rel071-log/reruns/$TAG/illse
WW=${W/\/c\//C:/}
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$REL/dev/ciharden-out/Rob0.3.32-lapack/bin/x64/Rscript.exe"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-release/dev/stan-cache
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
rm -rf "$W"; mkdir -p "$W/falarm" "$W/base/falarm"
for arm in merged lane14; do
  mkdir -p "$W/s-$arm"
  l=$LIB; [ $arm = lane14 ] && l=C:/Users/adf44/source/r/wt-illse-lib
  for f in "$REL"/dev/illse-sf-*.R; do
    sed -e "s|C:/Users/adf44/source/r/wt-illse-lib|$l|g" \
        -e "s|C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-sf-|$WW/s-$arm/illse-sf-|g" \
        "$f" > "$W/s-$arm/$(basename "$f")"
  done
done
cp "$LANE"/dev/illse-log/sf-g2/falarm/*-base.tsv "$W/base/falarm/"
cp "$LANE"/dev/illse-log/sf-g2/bitwise-*-base.rds "$W/base/"
N=8
jobs_wait() { while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done; }
run() { # <cfg> <logname> <script> args...
  local c=$1 lg=$2; shift 2
  local R="$RS"
  if [ $c = ob ]; then R="$OB"; export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
  else unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS; fi
  jobs_wait
  echo "$(date +%T) START $c $lg :: $*" >> "$W/driver.txt"
  "$R" "$@" < /dev/null > "$W/$lg.txt" 2>&1 &
}
cd "$REL"
for c in ref ob; do
  for arm in merged lane14; do
    for g in nl lin poly mixed bounded; do
      run $c falarm/$g-$c-$arm "$W/s-$arm/illse-sf-rev3-falarm.R" lane $g \
        "$W/falarm/$g-$c-$arm.tsv" "1:20"
    done
    run $c bitwise-$c-$arm "$W/s-$arm/illse-sf-rev4-bitwise.R" lane \
      "$W/bitwise-$c-$arm.rds"
  done
done
wait
cp "$W"/base/falarm/*.tsv "$W/falarm/"
for c in ref ob; do
  for arm in merged lane14; do
    "$RS" "$REL/dev/seflags-rev4-bitcmp.R" "$W/bitwise-$c-$arm.rds" \
      "$W/base/bitwise-$c-base.rds" > "$W/bitcmp-$c-$arm-vs-base.txt" 2>&1
  done
done
for pair in "merged base" "lane14 base" "merged lane14"; do
  set -- $pair
  ILLSE_CFGS="ref ob" "$RS" "$REL/dev/illse-sf-falarm-cmp.R" "$W/falarm" "$1" \
    "$W/falarm" "$2" > "$W/falarm-cmp-$1-vs-$2.txt" 2>&1
done
echo "$(date +%T) all done" >> "$W/driver.txt"
