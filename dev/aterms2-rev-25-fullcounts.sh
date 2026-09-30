#!/usr/bin/env bash
# Reviewer, re-check: the whole ungated suite of core and every
# extension on the lane build, one file per process
# (dev/aterms2-rev-list-full.txt). Totals generated from the logs.
# Output: dev/aterms2-rev-log-25-fullcounts.txt
cd /c/Users/adf44/source/r/frmtmb-wt-aterms2/dev
for d in aterms2-rev-suite-full-*; do
  pkg=${d#aterms2-rev-suite-full-}
  nf=$(ls $d/*.txt | wc -l)
  nr=$(grep -l '^RESULT' $d/*.txt | wc -l)
  grep -h '^RESULT' $d/*.txt | sed 's/[a-z]*=//g' | awk -v p=$pkg -v nf=$nf -v nr=$nr '
    {pa+=$4; fa+=$5; er+=$6; sk+=$7; wa+=$8}
    END {printf "%-16s files %3d with RESULT %3d  pass %6d fail %d error %d skip %d warn %d\n", p, nf, nr, pa, fa, er, sk, wa}'
done
echo "files with fail, error or warn > 0:"
grep -h '^RESULT' aterms2-rev-suite-full-*/*.txt | awk '{split($5,a,"=");split($6,b,"=");split($8,c,"=");if (a[2]+b[2]+c[2]>0) print}'
echo "files with skip > 0:"
grep -H '^RESULT' aterms2-rev-suite-full-*/*.txt | awk '{split($7,a,"="); if (a[2]>0) print $0}' | sed 's#aterms2-rev-suite-full-##; s#/lane__[^:]*:RESULT lane# #'
echo "logs without RESULT:"
for f in aterms2-rev-suite-full-*/*.txt; do grep -q '^RESULT' $f || echo "  $f"; done
