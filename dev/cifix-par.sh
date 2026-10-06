#!/usr/bin/env bash
# Run test files one per process. Each job: <tag> <ref|ob> <lib> <pkg>
# <file>. Usage: dev/cifix-par.sh <jobfile>; logs go to
# dev/cifix-log/t-<tag>.log. dev/cifix-pardrv.sh is the throttled
# driver for whole suites.
cd /c/Users/adf44/source/r/frmtmb-wt-cifix
while read -r tag rbin lib pkg file; do
  [ -z "$tag" ] && continue
  if [ "$rbin" = ob ]; then
    R=./dev/cifix-out/Rob/bin/x64/Rscript.exe
    export OPENBLAS_NUM_THREADS=4
  else
    R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
  fi
  "$R" dev/cifix-run1.R "$lib" "$pkg" "$file" \
    > "dev/cifix-log/t-$tag.log" 2>&1 &
done < "$1"
wait
