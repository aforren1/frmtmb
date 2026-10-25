#!/usr/bin/env bash
# The ported brms-suite tier in record mode (dev/brmsport-run.R), on one
# arm, one R process per file, N at a time. dev/brmsport-record.sh does
# the same serially but writes into dev/brmsport-log/, which is tracked,
# so this lane writes to its own ignored directory instead.
#
#   bash dev/brmsnames-1024-port.sh base|lane [N]
#
# base reads rellib-r12; lane reads the lane library, which holds core
# and every extension built from this worktree. Output:
# dev/brmsnames-1024-port-out/<arm>/rec-<pkg>-<topic>.tsv and run-*.txt.
set -u
arm="$1"; N="${2:-4}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
if [ "$arm" = lane ]; then
  export FRMTMB_PORT_LIB=C:/Users/adf44/source/r/wt-brmsnames-lib
else
  export FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r12
fi
OUT="dev/brmsnames-1024-port-out/$arm"
mkdir -p "$OUT"
for f in tests/testthat/test-brms-suite-*.R \
         extensions/frmtmb.sample/tests/testthat/test-brms-suite-*.R; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  topic=$(basename "$f" .R | sed 's/^test-brms-suite-//')
  case "$f" in extensions/*) pkg=frmtmb.sample ;; *) pkg=frmtmb ;; esac
  "$RS" dev/brmsport-run.R "$pkg" "$f" "$OUT/rec-$pkg-$topic.tsv" \
    > "$OUT/run-$pkg-$topic.txt" 2>&1 &
done
wait
grep -h "^RESULT" "$OUT"/run-*.txt
