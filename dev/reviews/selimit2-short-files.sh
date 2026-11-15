#!/usr/bin/env bash
# Short review of lane selimit2's fix round: the six SE files on reference
# BLAS (lane build), then test-id-kron.R five times on each build through
# the lane's traced runner, one R process at a time (waits for any of
# this lane's runs already going).
#   bash selimit2-short-files.sh <out dir> [file ...]
set -u
OUT="$1"; shift
L=/c/Users/adf44/source/r/frmtmb-wt-selimit2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p "$OUT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$L/dev/stan-cache"
LANE=C:/Users/adf44/source/r/wt-selimit2-lib
mine() { # R processes of this review running through the lane's runners
  powershell -NoProfile -Command "(Get-CimInstance Win32_Process -Filter \"Name='Rscript.exe'\" | Where-Object { \$_.CommandLine -match 'frmtmb-wt-selimit2/dev/(release/run-tests|selimit2-run-traced)' -and \$_.CommandLine -notmatch 'PROGRA~1' } | Measure-Object).Count"
}
slot() { while [ "$(mine | tr -d '\r')" -ge 1 ]; do sleep 5; done; }
run() { # tag lib pkg file
  slot
  FRMTMB_RELLIB="$2" "$RS" "$L/dev/release/run-tests.R" "$3" "$L/$4" \
    > "$OUT/$1.txt" 2>&1 &
  sleep 5
}
runt() { # tag lib: test-id-kron.R through the traced runner
  slot
  SELIMIT2_TRACE="$OUT/$1.trace" FRMTMB_RELLIB="$2" \
    "$RS" "$L/dev/selimit2-run-traced.R" frmtmb \
    "$L/tests/testthat/test-id-kron.R" > "$OUT/$1.txt" 2>&1 &
  sleep 5
}
for f in "$@"; do
  run "lane--$f" "$LANE" frmtmb "tests/testthat/$f.R"
done
run "lane--eam-test-defects" "$LANE" frmtmb.eam \
  extensions/frmtmb.eam/tests/testthat/test-defects.R
for k in 1 2 3 4 5; do
  runt "idkron-lane-$k" "$LANE"
  runt "idkron-base-$k" ""
done
wait
echo done
