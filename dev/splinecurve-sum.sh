#!/bin/sh
# Lane splinecurve: sum the RESULT lines of a suite log, so a total in
# the findings is generated rather than typed. Strips the BOM and CRs
# PowerShell writes. Usage: sh dev/splinecurve-sum.sh <suite log>
tr -cd '\11\12\40-\176' < "$1" | grep "RESULT" |
  sed -E 's/.*pass=([0-9]+) fail=([0-9]+) err=([0-9]+) skip=([0-9]+)/\1 \2 \3 \4/' |
  awk -v f="$1" '{p+=$1; fl+=$2; e+=$3; s+=$4; n++}
    END {print f, "files", n, "pass", p, "fail", fl, "err", e, "skip", s}'
