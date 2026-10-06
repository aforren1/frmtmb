#!/bin/sh
# Reviewer round 2: sum the gated core tier's RESULT lines and name every
# file with a failure, an error, a skip, a warning or no RESULT line.
cd /c/Users/adf44/source/r/frmtmb-wt-gpby/dev/gpby-rev2-gated
echo "files $(ls *.txt | wc -l), with RESULT $(grep -l '^RESULT' *.txt | wc -l), lib lines naming wt-gpby-lib $(grep -l '^lib: C:/Users/adf44/source/r/wt-gpby-lib/frmtmb' *.txt | wc -l)"
grep -h '^RESULT' *.txt | sed 's/.*pass=\([0-9]*\) fail=\([0-9]*\) error=\([0-9]*\) skip=\([0-9]*\) warn=\([0-9]*\)/\1 \2 \3 \4 \5/' | awk '{p+=$1;f+=$2;e+=$3;s+=$4;w+=$5} END {print "pass="p" fail="f" error="e" skip="s" warn="w}'
grep -h '^RESULT' *.txt | grep -v 'fail=0 error=0 skip=0 warn=0'
for f in *.txt; do grep -q '^RESULT' $f || echo "NO RESULT $f"; done
