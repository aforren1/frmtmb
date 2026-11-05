#!/usr/bin/env bash
# inner shell of the spawn test: an R session that starts a child R process every 10 s
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/review-nlfix-spawn.R > dev/review-nlfix-spawn.log 2>&1
