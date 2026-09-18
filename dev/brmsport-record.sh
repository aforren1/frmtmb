#!/bin/sh
# Record every generated brms-suite file, one R process each, into
# dev/brmsport-log/rec-<pkg>-<topic>.tsv. Usage: sh dev/brmsport-record.sh [topic ...]
# The tree is a REQUIRED argument, FRMTMB_PORT_ROOT. The first spelling
# named the brmsport worktree, which no longer exists, so the script ran
# nothing and said so only by failing. A DEFAULT is not the fix: the
# obvious default is the main checkout, and dev/organizer-rules.md
# forbids a lane agent from writing there, so an unset variable would
# have put this lane's logs in main. Refusing is the guard failing
# closed.
if [ -z "$FRMTMB_PORT_ROOT" ]; then
  echo "FRMTMB_PORT_ROOT is unset: set it to the tree to record, never"
  echo "to the main checkout, which lane agents must not write to."
  exit 1
fi
cd "$FRMTMB_PORT_ROOT" || exit 1
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
