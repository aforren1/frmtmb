#!/usr/bin/env bash
# Punch round 2 final pass: both extensions' suites, plain and gated,
# base (rellib-r3) and this lane's final build (phase3b-lib), one file
# per process, two streams at a time. Labels plain-r2 and gated-r2 keep
# the round-1 logs. Marker: dev/phase3b-log/final-suites3.done
cd "$(dirname "$0")/.." || exit 1
B="$BASH"
gated() {
  ( export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true \
      FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
      FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache" \
      R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
    export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
    "$B" dev/phase3b-suite.sh "$1" "$2" C:/Users/adf44/source/r/phase3b-lib \
      gated-r2 )
}
plain() {
  "$B" dev/phase3b-suite.sh "$1" "$2" C:/Users/adf44/source/r/phase3b-lib \
    plain-r2
}
( plain frmtmb.eam new; plain frmtmb.learn new
  gated frmtmb.eam new; gated frmtmb.learn new ) &
( plain frmtmb.eam base; plain frmtmb.learn base
  gated frmtmb.eam base; gated frmtmb.learn base ) &
wait
date -Iseconds > dev/phase3b-log/final-suites3.done
