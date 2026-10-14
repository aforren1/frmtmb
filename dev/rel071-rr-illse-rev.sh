#!/usr/bin/env bash
# The illse reviewer's exact-reference scripts (dev/illse-rev-*.R,
# dev/illse-rev2-*.R; 256-bit references written without the lane's
# helpers) on a build of the release tree, as copies with the lane
# library replaced by the arm's and the lane's dev/ by the copies'
# directory. Arms "merged" (<lib>) and "lane14" (wt-illse-lib, the
# lane's final install) run the same scripts in the same session, so
# the two logs compare line by line after their first (arm) line; the
# lane's own logs (dev/illse-log/fx/) are from installs 13 and 14.
#   bash dev/rel071-rr-illse-rev.sh <lib> <tag>
set -u
LIB="$1"; TAG="$2"
REL=/c/Users/adf44/source/r/frmtmb-wt-release
W=$REL/dev/rel071-log/reruns/$TAG/illse-rev
WW=${W/\/c\//C:/}
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$REL/dev/ciharden-out/Rob0.3.32-lapack/bin/x64/Rscript.exe"
export NOT_CRAN=true
rm -rf "$W"; mkdir -p "$W"
for arm in merged lane14; do
  mkdir -p "$W/s-$arm"
  l=$LIB; [ $arm = lane14 ] && l=C:/Users/adf44/source/r/wt-illse-lib
  for f in "$REL"/dev/illse-rev-*.R "$REL"/dev/illse-rev2-*.R; do
    sed -e "s|C:/Users/adf44/source/r/wt-illse-lib|$l|g" \
        -e "s|C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev|$WW/s-$arm/illse-rev|g" \
        "$f" > "$W/s-$arm/$(basename "$f")"
  done
done
N=8
jobs_wait() { while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done; }
run() { # <cfg> <arm> <logname> <script> args...
  local c=$1 arm=$2 lg=$3 s=$4; shift 4
  local R="$RS"
  if [ $c = ob ]; then R="$OB"; export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
  else unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS; fi
  jobs_wait
  REV_ARM=lane "$R" "$W/s-$arm/$s" "$@" < /dev/null > "$W/$lg-$c-$arm.txt" 2>&1 &
}
cd "$REL"
for arm in merged lane14; do
  for c in ref ob; do
    run $c $arm lmm8 illse-rev-lmm.R 101:103 8
    run $c $arm glm-binom illse-rev-glm.R binom 201:204
    run $c $arm edge4 illse-rev-edge4.R
  done
  run ref $arm lmm7 illse-rev-lmm.R 101:103 7
  run ref $arm insight illse-rev-insight.R
  run ref $arm edge8 illse-rev2-edge8.R
done
wait
for f in "$W"/*-merged.txt; do
  g=${f%-merged.txt}-lane14.txt
  if diff -q <(tail -n +2 "$f") <(tail -n +2 "$g") > /dev/null; then
    echo "SAME merged and lane14 after line 1: $(basename "$f") ($(grep -c '' "$f") lines)"
  else
    echo "DIFFERENT: $(basename "$f") ($(diff <(tail -n +2 "$f") <(tail -n +2 "$g") | grep -c '^[<>]') lines)"
  fi
done
