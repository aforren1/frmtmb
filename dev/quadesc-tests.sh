#!/usr/bin/env bash
# Lane quadesc: one test file per R process under a chosen R build, the
# lane library ahead of rellib-r22 (or rellib-r22 alone for "base").
#   bash dev/quadesc-tests.sh <ref|0.3.26|0.3.32> <lib|base> <tag> <pkg> <file>...
# Logs dev/quadesc-log/t-<tag>-<file>.txt, then the RESULT lines.
set -u
VER="$1"; ARM="$2"; TAG="$3"; PKG="$4"; shift 4
cd C:/Users/adf44/source/r/frmtmb-wt-quadesc
export FRMTMB_BRMS_FIT_TESTS=true
if [ "$ARM" = lib ]; then
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-quadesc-lib
else
  unset FRMTMB_RELLIB
fi
for f in "$@"; do
  bash dev/quadesc-run.sh "$VER" default "t-$TAG-$(basename "$f" .R)" \
    dev/release/run-tests.R "$PKG" "$f"
  grep -h "^lib\|^RESULT" "dev/quadesc-log/t-$TAG-$(basename "$f" .R).txt"
done
