#!/usr/bin/env bash
# Punch round 1 final pass: both extensions' suites, plain and gated,
# base and this lane's final build (phase3b-lib), one file per process.
# TWO streams at a time, after the 09:24 crash. Marker:
# dev/phase3b-log/final-suites.done
cd "$(dirname "$0")/.." || exit 1
B="$BASH"
( "$B" dev/phase3b-suite.sh frmtmb.eam new
  "$B" dev/phase3b-suite.sh frmtmb.learn new
  "$B" dev/phase3b-gated.sh frmtmb.eam new
  "$B" dev/phase3b-gated.sh frmtmb.learn new ) &
( "$B" dev/phase3b-suite.sh frmtmb.eam base
  "$B" dev/phase3b-suite.sh frmtmb.learn base
  "$B" dev/phase3b-gated.sh frmtmb.eam base
  "$B" dev/phase3b-gated.sh frmtmb.learn base ) &
wait
date -Iseconds > dev/phase3b-log/final-suites.done
