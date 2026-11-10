#!/usr/bin/env bash
# Start K copies of dev/equiv-live-row.R at the same instant, R rounds,
# and record each child's exit status, stderr, live status and the
# extensions it loaded. No parent R process: K R processes at once.
#   bash dev/compat-m-race.sh <K> <rounds> <outdir> <id,id,...>
set -u
K=$1; R=$2; out=$3; ids=$4
cd "$(dirname "$0")/.."
export FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-compat-lib;C:/Users/adf44/source/r/rellib-r19"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p "$out"
for r in $(seq 1 "$R"); do
  for k in $(seq 1 "$K"); do
    ( "$RS" dev/equiv-live-row.R "$ids" "$out/r$r-k$k.tsv" \
        > "$out/r$r-k$k.out" 2> "$out/r$r-k$k.err"
      echo "rc=$?" >> "$out/r$r-k$k.err" ) &
  done
  wait
done
