#!/usr/bin/env bash
# The fix round's mutants, one process each, in sequence (dev/fitedge-mutrun.R).
WT=/c/Users/adf44/source/r/frmtmb-wt-fitedge
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
for job in "prefix test-gp-start.R" "gn1 test-gp-start.R" "noothers test-gp-start.R" \
           "ridgekept test-nl-family.R" "control test-gp-start.R"; do
  set -- $job
  "$RS" "$WT/dev/fitedge-mutrun.R" "$1" "$2" > "$WT/dev/fitedge-out/r2/mut-$1-$2.txt" 2>&1
  grep -h "^RESULT\|^DETAIL" "$WT/dev/fitedge-out/r2/mut-$1-$2.txt" | cut -c1-260
done
