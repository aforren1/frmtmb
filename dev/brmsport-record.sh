#!/bin/sh
# Record every generated brms-suite file, one R process each, into
# dev/brmsport-log/rec-<pkg>-<topic>.tsv. Usage: sh dev/brmsport-record.sh [topic ...]
# the tree is the caller's: run from the worktree root
[ -f dev/brmsport-run.R ] || { echo "run me from the worktree root"; exit 1; }
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
mkdir -p dev/brmsport-log
for f in tests/testthat/test-brms-suite-*.R extensions/frmtmb.sample/tests/testthat/test-brms-suite-*.R; do
  topic=$(basename "$f" .R | sed 's/^test-brms-suite-//')
  if [ $# -gt 0 ]; then
    case " $* " in *" $topic "*) ;; *) continue ;; esac
  fi
  case "$f" in extensions/*) pkg=frmtmb.sample ;; *) pkg=frmtmb ;; esac
  "$R" dev/brmsport-run.R "$pkg" "$f" "dev/brmsport-log/rec-$pkg-$topic.tsv" \
    > "dev/brmsport-log/run-$pkg-$topic.txt" 2>&1
  grep -h "^RESULT" "dev/brmsport-log/run-$pkg-$topic.txt" || echo "NO RESULT $pkg $topic"
done
