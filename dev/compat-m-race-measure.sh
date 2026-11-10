#!/usr/bin/env bash
# Start one dev/equiv-measure.R per id at the same instant, R rounds,
# each round into its own folder; record each parent's exit status.
# Each job is 2 R processes while its live-row child runs.
#   bash dev/compat-m-race-measure.sh <rounds> <outdir> <id> ...
set -u
R=$1; out=$2; shift 2
cd "$(dirname "$0")/.."
export FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-compat-lib;C:/Users/adf44/source/r/rellib-r19"
export FRMTMB_STAN_CACHE="$PWD/dev/stan-cache"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for r in $(seq 1 "$R"); do
  d="$out/r$r"; mkdir -p "$d"
  for id in "$@"; do
    ( "$RS" dev/equiv-measure.R "$id" "$d" --force > "$d/$id.log" 2>&1
      echo "rc=$?" >> "$d/$id.log" ) &
  done
  wait
  for id in "$@"; do
    echo "r$r $id $(grep -a -E '^verdict|^Error|rc=' "$d/$id.log" | tr '\n' ' ')"
  done
done
