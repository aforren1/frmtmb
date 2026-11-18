#!/usr/bin/env bash
# Reviewer of lane ordstarts, build round: the two reproduction queues
# (one R process each).
#   bash dev/reviews/ordstarts-build-review-rep.sh <base|lane>
cd "$(dirname "$0")/../.." || exit 1
export NOT_CRAN=true
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/reviews
q=$1
for des in cs1000 cc300; do
  f=$O/ordstarts-build-review-rep-$q-$des.tsv
  rm -f "$f"
  t0=$(date +%s)
  "$R" $O/ordstarts-build-review-rep.R "$q" "$des" 1:10 "$f" \
    > $O/ordstarts-build-review-rep-$q-$des.log 2>&1
  echo "JOB $q $des exit $? wall $(( $(date +%s) - t0 ))s" >> $O/ordstarts-build-review-rep-jobs.txt
done
