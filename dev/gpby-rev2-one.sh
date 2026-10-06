#!/bin/sh
# Reviewer round 2: run one core test file gated; $1 is the file path.
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
"$R" dev/gpby-rev2-runtest-g.R frmtmb "$1" lane > "dev/gpby-rev2-gated/$(basename "$1" .R).txt" 2>&1
