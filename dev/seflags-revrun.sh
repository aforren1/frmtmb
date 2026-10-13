#!/usr/bin/env bash
# Lane seflags, punch round 1: the reviewer's scripts
# (dev/seflags-rev-*.R) on the lane build and the base, under one R
# build, into dev/seflags-log/rev-<tag>/; the false-alarm battery with
# its outputs moved to dev/seflags-log/falarm-<tag>/.
#   bash dev/seflags-revrun.sh <tag> <ref|ob26>
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
tag="$1"; cfg="$2"
case "$cfg" in
  ref) RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe" ;;
  ob26) RS="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
        export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 ;;
esac
L="$W/dev/seflags-log/rev-$tag"
F="$W/dev/seflags-log/falarm-$tag"
mkdir -p "$L" "$F"
for s in pred pred2 pred3 join sd sd2 fx1theta deg7 neg year; do
  for a in lane base; do
    "$RS" dev/seflags-rev-$s.R $a > "$L/$s-$a.txt" 2>&1 &
  done
done
wait
for g in nl lin poly mixed bounded; do
  for a in lane base; do
    "$RS" dev/seflags-rev-falarm.R $a $g "$F/$g-$a.tsv" "1:20" \
      > "$F/$g-$a.txt" 2>&1 &
  done
done
wait
echo REVRUN DONE
