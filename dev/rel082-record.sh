#!/usr/bin/env bash
# Bin 1 of the ported brms suite in record mode on the release build:
# dev/defects-record.sh's loop (every file in its own R process, 3 at once,
# the same outputs, dev/brmsport-log/rec-<pkg>-<topic>.tsv and
# run-<pkg>-<topic>.txt), run under bash with Rtools LAST on PATH, which
# the 0.78.0 consolidation found keeps an R child from segfaulting at
# exit (dev/round-handoff.md, 0.78.0). Every extension is loaded.
#   bash dev/rel082-record.sh <topic ...>
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-rel082 || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r21
export FRMTMB_PORT_EXTENSIONS=all
mkdir -p dev/brmsport-log
jobs=""
# at most three R processes at once (the round's limit)
running=0
for f in tests/testthat/test-brms-suite-*.R extensions/frmtmb.sample/tests/testthat/test-brms-suite-*.R; do
  topic=$(basename "$f" .R | sed 's/^test-brms-suite-//')
  [ "$topic" = helper-copy ] && continue
  case " $* " in *" $topic "*) ;; *) continue ;; esac
  case "$f" in extensions/*) pkg=frmtmb.sample ;; *) pkg=frmtmb ;; esac
  "$R" dev/brmsport-run.R "$pkg" "$f" "dev/brmsport-log/rec-$pkg-$topic.tsv" \
    > "dev/brmsport-log/run-$pkg-$topic.txt" 2>&1 &
  jobs="$jobs $pkg-$topic"
  running=$((running + 1))
  if [ "$running" -ge "${REC_JOBS:-3}" ]; then wait -n; running=$((running - 1)); fi
done
wait
for j in $jobs; do
  grep -h "^RESULT" "dev/brmsport-log/run-$j.txt" || echo "NO RESULT $j"
  grep -h "extensions loaded" "dev/brmsport-log/run-$j.txt" | head -1
done
