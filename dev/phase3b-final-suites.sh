#!/usr/bin/env bash
# The final pass: both extensions' suites, plain and gated, on the base
# build and on this lane's final build, one file per process. Four
# streams in parallel; each writes its own files, so no two writers
# share a log. Marker: dev/phase3b-log/final-suites.done
cd "$(dirname "$0")/.." || exit 1
rm -f dev/phase3b-log/final-suites.done
B="$BASH"
( "$B" dev/phase3b-suite.sh frmtmb.eam base
  "$B" dev/phase3b-suite.sh frmtmb.learn base ) &
( "$B" dev/phase3b-suite.sh frmtmb.eam new
  "$B" dev/phase3b-suite.sh frmtmb.learn new ) &
( "$B" dev/phase3b-gated.sh frmtmb.eam base
  "$B" dev/phase3b-gated.sh frmtmb.learn base ) &
( "$B" dev/phase3b-gated.sh frmtmb.eam new
  "$B" dev/phase3b-gated.sh frmtmb.learn new ) &
wait
date -Iseconds > dev/phase3b-log/final-suites.done
