#!/usr/bin/env bash
# Each lane's key table on the release build (rellib-r19), one R process at
# a time, through dev/rel080-rerun.R. Output: dev/rel080-log/rr/<name>.txt.
#   bash dev/rel080-keytables.sh <batch>
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel080/dev/rel080-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel080 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-rel080/dev/stan-cache
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/rel080-log/rr
mkdir -p $O
run() {  # run <name> <script> [args]
  local n="$1"; shift
  echo "$(date '+%T') start $n" >> $O/progress.txt
  "$RS" dev/rel080-rerun.R "$@" > "$O/$n.txt" 2>&1
  echo "$(date '+%T') end $n exit $?" >> $O/progress.txt
}
case "$1" in
  b1)
    run loglik-routes dev/loglik-routes.R after
    run latperf-agree dev/latperf-agree.R
    run latperf-prof-one dev/latperf-prof.R lane 8000 one
    run latperf-prof-two dev/latperf-prof.R lane 8000 two
    run latperf-pat-zi dev/review-latperf-pattern2.R 8000 zipzi
    run latperf-pat-hu dev/review-latperf-pattern2.R 8000 hurdle
    run latperf-pat-hu500 dev/review-latperf-pattern2.R 8000 hurdle 500
    run bigx-fx2 dev/review-bigx-fx2.R
    run nlstart-fix3-inv dev/review-nlstart-fix3-inv.R
    run readers-sc2 dev/review-readers-sc2.R
    run readers-sc3 dev/review-readers-sc3.R
    run readers-sc3a dev/review-readers-sc3a.R
    run readers-inf dev/review-readers-inf.R
    ;;
  b2)
    # the lanes' output folders are copied into this tree first (gitignored),
    # so that the summaries find their before arms and nothing is written
    # into a lane tree
    run osa-time dev/review-osa-time.R
    run arstep-cost dev/review-arstep-cost.R p600,p1000,a300 2
    run nestquad-nodes-bern dev/nestquad-nodes.R bern asym
    BIGX_ARM=after run bigx-grid dev/bigx-grid.R dev/bigx-out/grid-r19.tsv
    run bigx-gridsum dev/bigx-summary.R dev/bigx-out/grid-r19.tsv
    run fitedge-gpf-fresh dev/review-fitedge-gpfresh.R fresh lane
    run fitedge-gpf-lane96 dev/review-fitedge-gpfresh.R lane96 lane
    for p in G N P; do run arstep-fresh-$p dev/review-arstep-fresh.R lane $p 0; done
    run arstep-ng dev/arstep-ng.R lane 1 pois_low,pois_mid,bern,binom5,negbin,gamma,lognorm,student,gsigma,gar 0
    for f in gaussian poisson bernoulli lognormal; do
      MFX_ARM=after run mfxse-rate-$f dev/mfxse-rate.R $f dev/mfxse-rate-out/r19-$f.tsv
    done
    run mfxse-ratesum dev/mfxse-ratesum.R r19
    ;;
  b3)
    mkdir -p $O/nlstart-scan
    for f in cumulative sratio cratio acat hurdle_cumulative bernoulli poisson gaussian gamma; do
      NLSTART_ARM=after run nlstart-scan/$f dev/nlstart-scan.R $f
    done
    run nlstart-scansum dev/nlstart-scansum.R dev/nlstart-out/scan-before $O/nlstart-scan
    mkdir -p $O/nlfresh $O/nlfresh2
    for fs in gaussian,poisson,bernoulli binomial gamma,cumulative; do
      REV_ARM=after run nlfresh/$fs dev/review-nlstart-fresh.R $fs ref
    done
    REV_BODIES=logis3,gomp,sig4 REV_ARM=after run nlfresh2/all dev/review-nlstart-fresh2.R gaussian,poisson,bernoulli,binomial,gamma,cumulative ref
    ;;
  b4)
    # the final build (core 23:32:17, after the bigx x nlstart fix): every
    # table a nonlinear fit reaches, into rr-final/
    O=dev/rel080-log/rr-final
    mkdir -p $O/nlstart-scan $O/nlfresh $O/nlfresh2
    run nlstart-fix3-inv dev/review-nlstart-fix3-inv.R
    run fitedge-gpf-fresh dev/review-fitedge-gpfresh.R fresh lane
    run fitedge-gpf-lane96 dev/review-fitedge-gpfresh.R lane96 lane
    run readers-sc2 dev/review-readers-sc2.R
    run readers-sc3 dev/review-readers-sc3.R
    run readers-inf dev/review-readers-inf.R
    for f in cumulative sratio cratio acat hurdle_cumulative bernoulli poisson gaussian gamma; do
      NLSTART_ARM=after run nlstart-scan/$f dev/nlstart-scan.R $f
    done
    run nlstart-scansum dev/nlstart-scansum.R dev/nlstart-out/scan-before $O/nlstart-scan
    for fs in gaussian,poisson,bernoulli binomial gamma,cumulative; do
      REV_ARM=after run nlfresh/$fs dev/review-nlstart-fresh.R $fs ref
    done
    REV_BODIES=logis3,gomp,sig4 REV_ARM=after run nlfresh2/all dev/review-nlstart-fresh2.R gaussian,poisson,bernoulli,binomial,gamma,cumulative ref
    ;;
  *) echo "unknown batch"; exit 1 ;;
esac
echo "$(date '+%T') BATCH $1 END" >> $O/progress.txt
