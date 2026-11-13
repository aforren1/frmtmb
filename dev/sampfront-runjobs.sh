#!/usr/bin/env bash
# Lane sampfront: run test files one per R process, at most N at once.
#   bash dev/sampfront-runjobs.sh <job list> <out dir> [N]
# A job line is "<package> <test file> [base]"; "base" drops the lane
# library (dev/sampfront-run-test.R). Each process writes its own log;
# nothing tails them. The summary is <out dir>/summary.txt.
set -u
list="$1"; out="$2"; N="${3:-3}"
ROOT="C:/Users/adf44/source/r/frmtmb-wt-sampfront"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
: "${TMP:?TMP is unset}"
mkdir -p "$out"
while read -r p f arm; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  o="$out/$p--$(basename "$f" .R)${arm:+-$arm}.txt"
  # no stdin for R: it would read the job list this loop reads
  ( "$RS" "$ROOT/dev/sampfront-run-test.R" "$p" "$ROOT/$f" $arm > "$o" 2>&1 \
      < /dev/null ) &
done < "$list"
wait
: > "$out/summary.txt"
while read -r p f arm; do
  [ -n "$p" ] || continue
  o="$out/$p--$(basename "$f" .R)${arm:+-$arm}.txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  echo "${arm:-lane} ${hit:-NO RESULT $f}" >> "$out/summary.txt"
done < "$list"
cat "$out/summary.txt"
