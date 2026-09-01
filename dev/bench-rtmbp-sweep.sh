#!/bin/sh
## Interleaved calls-only sweep. Configurations are visited in rotation and
## the whole rotation is repeated, so thermal drift on this laptop spreads
## across configurations instead of landing on whichever one ran last.
## usage: sh dev/bench-rtmbp-sweep.sh <shape> <pass>
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
SHAPE=$1
PASS=$2
run() { "$R" dev/bench-rtmbp-core.R "$1" "$2" "$SHAPE" "$3" calls "$PASS" 2>&1 |
          grep -Ev '^\[1\]|^attr|logarithm'; }
run RTMB  1 FALSE
run RTMBp 1 FALSE
run RTMBp 1 TRUE
run RTMBp 2 TRUE
run RTMBp 4 TRUE
run RTMBp 8 TRUE
run RTMBp 16 TRUE
