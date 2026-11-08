#!/usr/bin/env bash
# The release review's fixes, the rest, 1 R process at a time beside
# dev/rel080-chain4.sh: nlstart's 450 rows and the reviewer's fresh rows,
# the bigx grid and dev/bigx-prod.R, OpenBLAS on dev/rel080-log/ob2.jobs,
# R CMD check on core, the ledger, the RECORD=1 survey.
#   bash dev/rel080-chain5.sh
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel080/dev/rel080-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel080 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-rel080/dev/stan-cache
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/rel080-log
O=$L/rr-fix3
mkdir -p $O/nlstart-scan $O/nlfresh $O/nlf1 $O/nlf2
note() { echo "$(date '+%F %T') CHAIN5 $*" >> $L/chain.txt; }
run() { local n="$1"; shift; "$RS" dev/rel080-rerun.R "$@" > "$O/$n.txt" 2>&1; }
note "nl start"
for f in cumulative sratio cratio acat hurdle_cumulative bernoulli poisson gaussian gamma; do
  NLSTART_ARM=after run nlstart-scan/$f dev/nlstart-scan.R $f
done
run nlstart-scansum dev/nlstart-scansum.R dev/nlstart-out/scan-before $O/nlstart-scan
for fs in gaussian,poisson,bernoulli binomial gamma,cumulative; do
  REV_ARM=after run nlfresh/$fs dev/review-nlstart-fresh.R $fs ref
  cp "$O/nlfresh/$fs.txt" "$O/nlf1/$fs.log"
done
REV_BODIES=logis3,gomp,sig4 REV_ARM=after run nlfresh2 dev/review-nlstart-fresh2.R gaussian,poisson,bernoulli,binomial,gamma,cumulative ref
cp $O/nlfresh2.txt $O/nlf2/all.log
"$RS" dev/review-nlstart-freshsum.R $O/nlf1,$O/nlf2 dev/review-nlstart-out/fresh-before,dev/review-nlstart-out/fresh2-before > $O/nlstart-freshsum.txt 2>&1
note "nl end"
note "bigx start"
BIGX_ARM=after run bigx-grid dev/bigx-grid.R dev/bigx-out/grid-fix3.tsv
run bigx-gridsum dev/bigx-summary.R dev/bigx-out/grid-fix3.tsv
BIGX_ARM=after run bigx-prod dev/bigx-prod.R
note "bigx end"
note "ob start"
bash dev/ciharden-openblas.sh 0.3.26 lapack > $L/ob26-build3.txt 2>&1
bash dev/ciharden-openblas.sh 0.3.32 lapack > $L/ob32-build3.txt 2>&1
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r19
for c in ref ob0.3.26-lapack-t4 ob0.3.32-lapack-t4; do
  bash dev/ciharden-scan.sh r080-ob2-${c%%-lapack*} $c base 1 $L/ob2.jobs > /dev/null 2>&1
done
note "ob end"
# R CMD check runs two R processes: after chain4 has ended
until grep -q "CHAIN4 END" $L/chain.txt; do sleep 60; done
note "check start"
bash dev/rel080-check.sh frmtmb > $L/check-3.txt 2>&1
note "check end"
note "ledger start"
REC_JOBS=1 bash dev/rel080-ledger.sh > $L/ledger-run4.txt 2>&1
"$RS" dev/rel080-ledger-diff.R > $L/ledger-diff.txt 2>&1
"$RS" dev/rel080-bin2-diff.R > $L/bin2-diff.txt 2>&1
note "ledger end"
note "survey start"
EQUIV_JOBS=1 bash dev/rel080-survey.sh
"$RS" dev/rel080-equiv-delta.R > $L/equiv-delta.txt 2>&1
note "survey end"
note END
