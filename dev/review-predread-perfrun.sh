#!/usr/bin/env bash
# Item 3: the timing arms in alternating fresh processes, three rounds, and the layered
# route's own cost. Runs alone on the review's share of the machine.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-predread
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
O=dev/review-predread-out
for r in 1 2 3; do
  REV_LIB=before "$RS" dev/review-predread-perf.R $O/perf-before-$r.tsv > $O/perf-before-$r.log 2>&1
  "$RS" dev/review-predread-perf.R $O/perf-after-$r.tsv > $O/perf-after-$r.log 2>&1
done
"$RS" dev/review-predread-perf-cs.R > $O/perf-cs.log 2>&1
echo finished > $O/perf.done
