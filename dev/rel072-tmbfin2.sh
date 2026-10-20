#!/usr/bin/env bash
# Consolidation item 2, with lane illse2's own driver: dev/illse2-p4wide.R
# 150 2000 1:2 (the reviewer's reproduction, rc-wide-lane-b.txt), its
# setup (dev/illse2-common.R) pointed at the library named here. Each
# run is a fresh session; the count is of "promise already under
# evaluation" lines on stderr.
#   bash dev/rel072-tmbfin2.sh <lib> <tag> <runs>
set -u
lib="$1"; tag="$2"; runs="${3:-2}"
REL=C:/Users/adf44/source/r/frmtmb-wt-release
W=$REL/dev/rel072-log/tmbfin2-$tag; mkdir -p "$W"
sed "s#^lib_lane <- .*#lib_lane <- \"$lib\"#" $REL/dev/illse2-common.R > "$W/common.R"
sed "s#^source(.*illse2-common.R\")#source(\"$W/common.R\")#" $REL/dev/illse2-p4wide.R > "$W/p4wide.R"
for i in $(seq 1 "$runs"); do
  ILLSE2_ARM=lane "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$W/p4wide.R" 150 2000 1:2 \
    > "$W/run$i.out" 2> "$W/run$i.err"
  echo "$tag run $i: $(head -1 "$W/run$i.out") | promise lines $(grep -c 'promise already under evaluation' "$W/run$i.err")"
done
