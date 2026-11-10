#!/usr/bin/env bash
# Measure the ten pairs of the compat lane, three R processes at once.
#   bash dev/compat-equiv-run.sh <stock|live> <lib;lib> <outdir>
# stock runs dev/equiv-measure.R (scores against the frozen sample
# status); live runs dev/compat-equiv-live.R (scores against the build's
# own frm_compat()). One log per pair in <outdir>.
set -u
mode=$1
export FRMTMB_EQUIV_LIB=$2
out=$3
cd "$(dirname "$0")/.."
export FRMTMB_STAN_CACHE="$PWD/dev/stan-cache"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p "$out"
ids="p167 p158 p171 p169 p166 p014 p146 p010 p011 p059"
run_one() {
  local id=$1
  if [ "$mode" = stock ]; then
    "$RS" dev/equiv-measure.R "$id" "$out" --force > "$out/$id.log" 2>&1
  else
    "$RS" dev/compat-equiv-live.R "$id" "$out" > "$out/$id.log" 2>&1
  fi
  echo "done $id exit $?"
}
n=0
for id in $ids; do
  run_one "$id" &
  n=$((n + 1))
  if [ $((n % 3)) -eq 0 ]; then wait; fi
done
wait
grep -h "^verdict" "$out"/*.log
