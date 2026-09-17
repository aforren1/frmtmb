#!/bin/sh
# Reviewer: every prior-touching test file, one R process each.
cd /c/Users/adf44/source/r/frmtmb-wt-priorform
while read pkg filt; do
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/priorform-rev-priortrace-run1.R "$pkg" "$filt" 2>&1 | grep "^PRIORTRACE" | sed "s|^|$pkg |" >> dev/priorform-rev-priortrace-log.txt || echo "$pkg $filt NOLINE" >> dev/priorform-rev-priortrace-log.txt
done < dev/priorform-rev-priortrace-files.txt
echo DONE >> dev/priorform-rev-priortrace-log.txt
