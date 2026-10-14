#!/usr/bin/env bash
# Release reviewer of 0.71.0: lane illse's false-alarm battery (the
# lane seflags reviewer's script, 5 groups) on four builds under the
# reference BLAS and OpenBLAS 0.3.26 with LAPACK (4 threads, the
# ubuntu-24.04 runner), which the consolidation did not run (it ran
# 0.3.32). Arms: merged (rellib-r9), base (rellib-r8, 0.70.0), optstart
# (wt-optstart-lib), illse (wt-illse-lib).
#   bash dev/relrev071-falarm.sh <tag> <seeds R expr>
set -u
TAG="$1"; SEEDS="$2"
REL=/c/Users/adf44/source/r/frmtmb-wt-release
W=$REL/dev/relrev071-log/falarm-$TAG
WW=${W/\/c\//C:/}
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$REL/dev/relrev071-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
rm -rf "$W"; mkdir -p "$W/falarm"
declare -A LIBS=(
  [merged]=C:/Users/adf44/source/r/rellib-r9
  [base]=C:/Users/adf44/source/r/rellib-r8
  [optstart]=C:/Users/adf44/source/r/wt-optstart-lib
  [illse]=C:/Users/adf44/source/r/wt-illse-lib)
ARMS="${ARMS:-merged base optstart illse}"
for arm in $ARMS; do
  mkdir -p "$W/s-$arm"
  for f in "$REL"/dev/illse-sf-rev3-*.R; do
    sed -e "s|C:/Users/adf44/source/r/wt-illse-lib|${LIBS[$arm]}|g" \
        -e "s|C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-sf-|$WW/s-$arm/illse-sf-|g" \
        "$f" > "$W/s-$arm/$(basename "$f")"
  done
done
N=8
jobs_wait() { while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done; }
run() { # <cfg> <logname> <script> args...
  local c=$1 lg=$2; shift 2
  local R="$RS"
  if [ $c = ob ]; then R="$OB"; export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
  else unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS; fi
  jobs_wait
  echo "$(date +%T) START $c $lg" >> "$W/driver.txt"
  "$R" "$@" < /dev/null > "$W/$lg.txt" 2>&1 &
}
cd "$REL"
for c in ${CFGS:-ref ob}; do
  for arm in $ARMS; do
    for g in nl lin poly mixed bounded; do
      run $c falarm/$g-$c-$arm "$W/s-$arm/illse-sf-rev3-falarm.R" lane $g \
        "$W/falarm/$g-$c-$arm.tsv" "$SEEDS"
    done
  done
done
wait
echo "$(date +%T) all done" >> "$W/driver.txt"
