#!/bin/sh
# Lane wt-defects: dev/brmsport-record.sh with every file in its own R
# process AT ONCE, which the 63 GB machine allows; the sequential
# original takes the sum of the file times. Same outputs:
# dev/brmsport-log/rec-<pkg>-<topic>.tsv and run-<pkg>-<topic>.txt.
#   FRMTMB_PORT_LIB=<lib> sh dev/defects-record.sh [topic ...]
cd "$(dirname "$0")/.." || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
: "${FRMTMB_PORT_LIB:?set FRMTMB_PORT_LIB to the library under test}"
export FRMTMB_PORT_LIB
mkdir -p dev/brmsport-log
jobs=""
for f in tests/testthat/test-brms-suite-*.R extensions/frmtmb.sample/tests/testthat/test-brms-suite-*.R; do
  topic=$(basename "$f" .R | sed 's/^test-brms-suite-//')
  [ "$topic" = helper-copy ] && continue
  if [ $# -gt 0 ]; then
    case " $* " in *" $topic "*) ;; *) continue ;; esac
  fi
  case "$f" in extensions/*) pkg=frmtmb.sample ;; *) pkg=frmtmb ;; esac
  "$R" dev/brmsport-run.R "$pkg" "$f" "dev/brmsport-log/rec-$pkg-$topic.tsv" \
    > "dev/brmsport-log/run-$pkg-$topic.txt" 2>&1 &
  jobs="$jobs $pkg-$topic"
done
wait
for j in $jobs; do
  grep -h "^RESULT" "dev/brmsport-log/run-$j.txt" || echo "NO RESULT $j"
done
