#!/bin/bash
cd /c/Users/adf44/source/r/frmtmb-wt-latres
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export REV_TK_RESTRICT=1
mkdir -p dev/review-latres-out/tk2
sed 's#review-latres-out/tk/#review-latres-out/tk2/#' dev/review-latres-out/tk.jobs > dev/review-latres-out/tk2.jobs
xargs -P 7 -L 1 "$R" dev/review-latres-tk.R < dev/review-latres-out/tk2.jobs > dev/review-latres-out/tk2.log 2>&1
