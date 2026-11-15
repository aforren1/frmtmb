#!/usr/bin/env bash
# Every generated bin-2 file, tier mode (verdicts asserted), one at a
# time, against the lane library over rellib-r21, with this tree's
# helpers and adapters; FRMTMB_PORT_ROOT unset on purpose (item 4448).
cd "$(dirname "$0")/.."
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
unset FRMTMB_PORT_ROOT
export FRMTMB_PORT_LIB="C:/Users/adf44/source/r/wt-hygiene-lib;C:/Users/adf44/source/r/rellib-r21"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/hygiene-out/bin2
for f in dev/brmsport-bin2/*/test-brms-suite-bin2-*.R; do
  pkg=$(basename "$(dirname "$f")")
  t=$(basename "$f" .R)
  HYG_STAGE=dev/hygiene-out/bin2-stage "$RS" dev/hygiene-bin2-run.R "$pkg" "$f" \
    > "dev/hygiene-out/bin2/$pkg-$t.txt" 2>&1
  grep -h "^lib:\|^RESULT\|^FAILING\|^  " "dev/hygiene-out/bin2/$pkg-$t.txt"
done
rm -rf dev/hygiene-out/bin2-stage
