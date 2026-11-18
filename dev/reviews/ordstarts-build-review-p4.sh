#!/usr/bin/env bash
# Reviewer of lane ordstarts, build round: probe 4 queue (one R process).
#   bash dev/reviews/ordstarts-build-review-p4.sh <base|lane>
cd "$(dirname "$0")/../.." || exit 1
export NOT_CRAN=true
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/reviews
q=$1
f=$O/ordstarts-build-review-p4-$q.tsv
rm -f "$f"
"$R" $O/ordstarts-build-review-p4.R "$q" compat2 400 1 "$f" > $O/ordstarts-build-review-p4-$q.log 2>&1
for cfg in "cc 300" "cs 300" "cs 1000"; do
  set -- $cfg
  "$R" $O/ordstarts-build-review-p4.R "$q" $1 $2 1:40 "$f" >> $O/ordstarts-build-review-p4-$q.log 2>&1
done
echo "P4 $q done" >> $O/ordstarts-build-review-p4-$q.log
