#!/usr/bin/env bash
# The review's off-grid designs on the final build, into dev/impsamp-out/offgrid/
# (its scripts, its designs, other output directories): bash dev/impsamp-offgrid.sh
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/impsamp-out/offgrid
mkdir -p impsamp-out/offgrid-log
sed "s#^out <- \"C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out/multid#out <- \"$O/multid#" review-impsamp-multid.R > impsamp-multid-run.R
DS=$(grep '^DESIGNS=' review-impsamp-accuracy.sh | sed 's/DESIGNS="//; s/"$//')
{
  for d in $DS; do for s in 1 2 3 4 5; do echo "review-impsamp-accuracy.R $d $s $O/accuracy"; done; done
  for d in zip hurdle betafam cens_gauss trunc_pois tiny_x_slope huge_group zero_wt scale1e5 pois1e5 year; do for s in 1 2 3; do echo "review-impsamp-accuracy2.R $d $s $O/accuracy2"; done; done
  while read -r a b; do echo "impsamp-multid-run.R $a $b"; done < review-impsamp-out/multid-jobs.txt
} > impsamp-out/offgrid-jobs.txt
wc -l impsamp-out/offgrid-jobs.txt
xargs -P "${1:-6}" -L 1 bash -c '"'"$R"'" "$0" "$1" "$2" "$3" > impsamp-out/offgrid-log/"$1"-"$2".log 2>&1' < impsamp-out/offgrid-jobs.txt
echo offgrid-done
