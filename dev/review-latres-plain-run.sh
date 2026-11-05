#!/bin/bash
cd /c/Users/adf44/source/r/frmtmb-wt-latres
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rm -f dev/review-latres-out/plain.txt
for arm in base lane base lane; do
  REV_ARM=$arm "$R" dev/review-latres-plain.R dev/review-latres-out/plain.txt > /dev/null 2>&1
done
