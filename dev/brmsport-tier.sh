#!/bin/sh
# Run the generated brms-suite tier as the gated tier runs it (verdicts
# asserted), one R process per file, against the lane's private library.
# Output: dev/brmsport-log/tier-<pkg>-<topic>.txt and a summary on stdout.
# the tree is the caller's: run from the worktree root
[ -f dev/brmsport-run.R ] || { echo "run me from the worktree root"; exit 1; }
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
n=0; ok=0
for f in tests/testthat/test-brms-suite-*.R extensions/frmtmb.sample/tests/testthat/test-brms-suite-*.R; do
  topic=$(basename "$f" .R | sed 's/^test-brms-suite-//')
  case "$f" in extensions/*) pkg=frmtmb.sample ;; *) pkg=frmtmb ;; esac
  n=$((n + 1))
  "$R" dev/brmsport-run.R "$pkg" "$f" > "dev/brmsport-log/tier-$pkg-$topic.txt" 2>&1
  line=$(grep -h "^RESULT" "dev/brmsport-log/tier-$pkg-$topic.txt")
  if [ -n "$line" ]; then ok=$((ok + 1)); echo "$line"; else echo "NO RESULT $pkg $topic"; fi
done
echo "TIER ran $ok of $n files"
