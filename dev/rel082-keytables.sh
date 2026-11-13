#!/usr/bin/env bash
# Round of 2026-11-11, consolidation: the twelve lanes' key tables, rerun on
# the merged build (rellib-r21, frmtmb 0.82.0) so that the release record
# can set each number beside the lane's own.
#
#   bash dev/rel082-keytables.sh [table ...]
#
# Tables, in the order they run when none is named (cheap ones first):
#   quadseam boundfam-held ppdraw wcheck nlstop ordstart silent priortab
#   sampfront arspat boundfam-suite gpstop ordstart-sweep
# and, named only, gpstop-attrib (a probe of the gpstop rows that moved).
#
# The lane scripts are copies of the lanes' dev/ files and stay as the
# lanes wrote them. Where one hard-codes a lane library or sources a file
# by absolute path into a lane tree, a copy goes to
# dev/rel082-log/kt/scripts/ (ignored by git) with that path pointed at
# rellib-r21 or at the copy, and the copy runs from the tree root. The
# lane trees are only read: their reference logs and base traces.
# Output: dev/rel082-log/kt/<table>/. One R process at a time, because
# the consolidator's test tiers run beside this.
set -u
cd "$(dirname "$0")/.."
. dev/rel082-env.sh

SRC=C:/Users/adf44/source/r
REL=$SRC/frmtmb-wt-rel082
R21=$SRC/rellib-r21
R20=$SRC/rellib-r20
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
KT=dev/rel082-log/kt
SC=$KT/scripts
mkdir -p "$SC"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$REL/dev/stan-cache
# Several lane scripts print no frmtmb library line; this profile hook
# prints one at load, so every log can be checked for rellib-r21.
if [ ! -f "$SC/libhook.R" ]; then
  printf '%s\n' \
    'setHook(packageEvent("frmtmb", "onLoad"), function(pkgname, pkgpath) {' \
    '  cat("lib-hook: frmtmb",' \
    '      utils::packageDescription(pkgname, lib.loc = dirname(pkgpath))$Version,' \
    '      dirname(pkgpath), "\n")' \
    '})' > "$SC/libhook.R"
fi
export R_PROFILE_USER=$REL/$SC/libhook.R

# rr <log> <script> [args]: one R run, its exit code at the end of its log
rr() {
  local log=$1; shift
  mkdir -p "$(dirname "$log")"
  echo "run $(date +%H:%M:%S) $*" >&2
  "$RS" "$@" > "$log" 2>&1
  echo "EXIT $?" >> "$log"
}
# cp_lib <src> <dest> <lane>: copy with the lane library replaced by r21
cp_lib() {
  sed "s#$SRC/wt-$3-lib#$R21#g" "$1" > "$2"
}

t_quadseam() {
  local O=$KT/quadseam Q=$SRC/frmtmb-wt-quadseam/dev/reviews/quadseam-review-log
  mkdir -p "$O"
  cp_lib dev/reviews/quadseam-review-acc.R "$SC/quadseam-review-acc.R" quadseam
  sed "s#$SRC/frmtmb-wt-quadseam/dev/reviews/quadseam-review-acc.R#$SC/quadseam-review-acc.R#" \
    dev/quadseam-fr-sweep.R > "$SC/quadseam-fr-sweep.R"
  rr "$O/acc3-lane.txt" "$SC/quadseam-review-acc.R" lane \
    bern500,pois200,bern2sd3,poissd2 6:35
  awk -f dev/reviews/quadseam-short-short.awk "$Q/acc3-glmer.txt" \
    "$Q/acc3-base.txt" "$O/acc3-lane.txt" > "$O/acc3-short.txt"
  # worse or better than base: |exact - reported| per fit, lane against base
  awk -F'\t' '$1 == "ROW" && $5 != "ERROR" {
      k = $2 " " $3
      for (i = 5; i < NF; i += 2) if ($i == "dref") e = $(i + 1) + 0
      if (e < 0) e = -e
      v[$4, k] = e; if ($4 == "lane") K[k] = 1 }
    END { for (k in K) { n++; d = v["lane", k] - v["base", k]
        if (d > 1e-6) w++; if (d < -1e-6) b++
        if (v["lane", k] > ml) ml = v["lane", k]
        if (v["base", k] > mb) mb = v["base", k] }
      printf "fits %d lane worse than base by >1e-6 %d better %d max|err| lane %.3e base %.3e\n",
        n, w, b, ml, mb }' "$Q/acc3-base.txt" "$O/acc3-lane.txt" >> "$O/acc3-short.txt"
  rr "$O/fr-sweep-lane.txt" "$SC/quadseam-fr-sweep.R" "$R21" \
    bern5,bern500,bern20big,bern2sd3,bin500,bin50,binz,pois8,pois200,poissd2 1:5 1e-4,1e-5
  awk '$1 == "FIT" { t = $5; n[t]++
      for (i = 6; i < NF; i++) {
        if ($i == "err") { e = $(i + 1) + 0; if (e < 0) e = -e; if (e > m[t]) m[t] = e }
        if ($i == "rounds" && $(i + 1) + 0 > 0) r[t]++
        if ($i == "warn" && $(i + 1) != "-") w[t]++ } }
    END { for (t in n) printf "grad_tol %s fits %d warned %d rounds>0 %d max|err| %.3e\n",
      t, n[t], w[t], r[t], m[t] }' "$O/fr-sweep-lane.txt" > "$O/fr-sweep-sum.txt"
  rr "$O/runaway3-lane.txt" dev/quadseam-runaway3.R "$R21"
}

t_boundfam_held() {
  local O=$KT/boundfam
  mkdir -p "$O"
  rr "$O/held-r21.txt" dev/boundfam-held.R "$R21" "$O"
  rr "$O/held-r20.txt" dev/boundfam-held.R base "$O"
  # the first 200 values are fitted()'s Est.Error, the rest predict()'s
  paste "$O/held-lane.txt" "$O/held-base.txt" | awk '
    NR <= 200 { d = $1 - $2; if (d < 0) d = -d; if (d > m) m = d }
    END { printf "max |r21 - r20| fitted Est.Error over 200 rows %.3e\n", m }' \
    > "$O/held-diff.txt"
}

t_ppdraw() {
  local O=$KT/ppdraw P=$SRC/frmtmb-wt-ppdraw/dev
  mkdir -p "$O"
  rr "$O/fix-cases.txt" dev/ppdraw-fix-cases.R "$R21"
  rr "$O/fix-merged-b.txt" dev/ppdraw-fix-merged-b.R "$R21"
  rr "$O/sweep-lane.txt" dev/ppdraw-sweep.R "$R21" "$O/ppdraw-sweep-lane.tsv"
  sed -e "s#dev/ppdraw-sweep-base.tsv#$P/ppdraw-sweep-base.tsv#" \
      -e "s#dev/ppdraw-sweep-lane.tsv#$O/ppdraw-sweep-lane.tsv#" \
      dev/ppdraw-table.R > "$SC/ppdraw-table.R"
  rr "$O/table.txt" "$SC/ppdraw-table.R"
  diff "$P/ppdraw-table.txt" "$O/table.txt" > "$O/table.diff"
}

t_wcheck() {
  local O=$KT/wcheck-out W=$SRC/frmtmb-wt-wcheck/dev
  mkdir -p "$O"
  for p in mv cens prior; do
    rr "$O/cases-$p-lane.txt" dev/wcheck-cases.R "$R21" "$p"
    cp "$W/wcheck-out/cases-$p-base.txt" "$O/"
  done
  sed "s#dev/wcheck-out/#$O/#g" dev/wcheck-table.R > "$SC/wcheck-table.R"
  rr "$O/table.txt" "$SC/wcheck-table.R"
  diff "$W/wcheck-out/fr3/table.txt" "$O/table.txt" > "$O/table.diff"
  rr "$O/common-lane.txt" dev/reviews/wcheck-review-common.R "$R21"
  diff "$W/reviews/wcheck-review-common-base.txt" "$O/common-lane.txt" > "$O/common.diff"
  rr "$O/routes-lane.txt" dev/reviews/wcheck-short2-routes.R "$R21"
  diff "$W/wcheck-out/fr3/routes-lane.txt" "$O/routes-lane.txt" > "$O/routes.diff"
}

t_nlstop() {
  local O=$KT/nlstop N=$SRC/frmtmb-wt-nlstop/dev
  mkdir -p "$O"
  cp_lib dev/nlstop-common.R "$SC/nlstop-common.R" nlstop
  sed -i "s#$SRC/frmtmb-wt-nlstop/dev/stan-cache#$REL/dev/stan-cache#" "$SC/nlstop-common.R"
  for f in dev/nlstop-prior.R dev/nlstop-4850-ulp.R \
           dev/reviews/nlstop-review-intercept.R dev/reviews/nlstop-short-attack.R; do
    sed "s#source(\"dev/nlstop-common.R\")#source(\"$SC/nlstop-common.R\")#" \
      "$f" > "$SC/$(basename "$f")"
  done
  cp_lib "$N/nlstop-out/fix/zi-probe.R" "$SC/nlstop-zi-probe.R" nlstop
  export NLSTOP_ARM=after NOT_CRAN=true
  rr "$O/prior-after.txt" "$SC/nlstop-prior.R"
  # the last field of a ROW line is the gap
  awk '$1 == "ROW" { n++; g = $NF + 0; if (g < 0) g = -g; if (g > 1e-4) b++ }
    END { printf "ROW lines %d with |gap| > 1e-4 %d\n", n, b }' \
    "$O/prior-after.txt" > "$O/prior-sum.txt"
  rr "$O/4850-ulp-after.txt" "$SC/nlstop-4850-ulp.R"
  rr "$O/intercept-after.txt" "$SC/nlstop-review-intercept.R"
  rr "$O/attack-p1p2.txt" "$SC/nlstop-short-attack.R" p1 p2
  rr "$O/zi-probe.txt" "$SC/nlstop-zi-probe.R"
  unset NLSTOP_ARM NOT_CRAN
}

t_ordstart() {
  local O=$KT/ordstart OT=$SRC/frmtmb-wt-ordstart/dev
  mkdir -p "$O"
  cp_lib dev/ordstart-starts.R "$SC/ordstart-starts.R" ordstart
  for f in cum_cum cum_sratio; do
    rr "$O/lane-$f.log" "$SC/ordstart-starts.R" lane 1:20 "$O/lane-$f.tsv" $f 0
  done
  rr "$O/compare-40.txt" dev/ordstart-compare.R \
    "$OT/ordstart-out/base-cumcum.tsv" "$O/lane-cum_cum.tsv" \
    "$OT/ordstart-out/base-cumsr.tsv" "$O/lane-cum_sratio.tsv"
  diff "$OT/ordstart-q4-out/compare-40.txt" "$O/compare-40.txt" > "$O/compare-40.diff"
}

t_ordstart_sweep() {
  local O=$KT/ordstart/sweep OT=$SRC/frmtmb-wt-ordstart/dev/ordstart-q4-out/sweep
  mkdir -p "$O"
  cp_lib dev/reviews/ordstart-review-sweep.R "$SC/ordstart-review-sweep.R" ordstart
  for s in sweep wide steep; do
    rr "$O/lane-$s.log" "$SC/ordstart-review-sweep.R" lane $s "$O/lane-$s"
  done
  rr "$O/sweep-sum.txt" dev/ordstart-fix-sweep-sum.R \
    "$OT/base-sweep" "$O/lane-sweep" "$OT/base-wide" "$O/lane-wide" \
    "$OT/base-steep" "$O/lane-steep"
}

t_silent() {
  local O=$KT/silent
  mkdir -p "$O"
  cp_lib dev/silent-dharma-unif.R "$SC/silent-dharma-unif.R" silent
  cp_lib dev/reviews/silent-short-cox.R "$SC/silent-short-cox.R" silent
  cp_lib dev/reviews/silent-review-cox.R "$SC/silent-review-cox.R" silent
  rr "$O/dharma-unif-lane.txt" "$SC/silent-dharma-unif.R" lane
  rr "$O/short-cox.txt" "$SC/silent-short-cox.R"
  rr "$O/review-cox-lane.txt" "$SC/silent-review-cox.R" lane
}

t_priortab() {
  local O=$KT/priortab
  mkdir -p "$O"
  PRIORTAB_LIB=$R21 rr "$O/compare.txt" dev/reviews/priortab-review-compare.R
  diff "$SRC/frmtmb-wt-priortab/dev/reviews/priortab-review-compare.txt" \
    "$O/compare.txt" > "$O/compare.diff"
}

t_sampfront() {
  local O=$KT/sampfront S=$SRC/frmtmb-wt-sampfront/dev
  mkdir -p "$O"
  for f in dev/review-sampfront-probe.R dev/review-sampfront-probe-real.R \
           dev/review-sampfront-subset.R dev/review-sampfront-subset2.R \
           dev/reviews/sampfront-short-design.R dev/reviews/sampfront-short-place.R; do
    cp_lib "$f" "$SC/$(basename "$f")" sampfront
    sed -i "s#$SRC/frmtmb-wt-sampfront/dev/stan-cache#$REL/dev/stan-cache#" \
      "$SC/$(basename "$f")"
  done
  rr "$O/probe-lane.txt" "$SC/review-sampfront-probe.R"
  rr "$O/probe-real-lane.txt" "$SC/review-sampfront-probe-real.R"
  rr "$O/subset-lane.txt" "$SC/review-sampfront-subset.R"
  rr "$O/subset2-lane.txt" "$SC/review-sampfront-subset2.R"
  rr "$O/design-lane.txt" "$SC/sampfront-short-design.R"
  rr "$O/place-lane.txt" "$SC/sampfront-short-place.R"
  # the lane's logs, without the library lines, timings and exit codes
  local strip='/^(lib|lib-hook|arm|frmtmb.sample|EXIT)/d; s/[0-9.]+ ?s(ec)?\b//g'
  for x in probe probe-real subset subset2; do
    diff <(sed -E "$strip" "$S/sampfront-out/fix/$x-lane.txt") \
         <(sed -E "$strip" "$O/$x-lane.txt") > "$O/$x.diff"
  done
  for x in design place; do
    diff <(sed -E "$strip" "$S/reviews/sampfront-short-out/$x-lane.txt") \
         <(sed -E "$strip" "$O/$x-lane.txt") > "$O/$x.diff"
  done
}

t_arspat() {
  local O=$KT/arspat
  mkdir -p "$O"
  export ARSPAT_LIB=$R21
  rr "$O/l1b-all.txt" dev/arspat-l1b.R
  # l1c's case list holds l1b's as well; the lane's log has the sar cases only
  rr "$O/l1c-sar.txt" dev/arspat-l1c.R mi_sar_error mi_sar_lag
  rr "$O/l1d-student.txt" dev/arspat-l1d.R \
    mi_fcor_student mi_sar_error_student mi_sar_lag_student
  rr "$O/l1d-fix.txt" dev/arspat-l1d.R \
    mi_fcor_student mi_sar_error_student mi_sar_lag_student \
    skew_se_ar skew_se_ar_cond skew_se_cosy
  unset ARSPAT_LIB
  sed "s#dev/arspat-out/#$O/#g" dev/arspat-l1-summary.R > "$SC/arspat-l1-summary.R"
  rr "$O/l1-summary.md" "$SC/arspat-l1-summary.R"
}

t_boundfam_suite() {
  local O=$KT/boundfam B=$SRC/frmtmb-wt-boundfam/dev
  mkdir -p "$O"
  sed "s#$SRC/frmtmb-wt-boundfam/dev/release/run-tests.R#$REL/dev/release/run-tests.R#" \
    dev/boundfam-run-traced2.R > "$SC/boundfam-run-traced2.R"
  sed -e "s#$SRC/wt-boundfam-lib#$R21#" \
      -e "s#^ROOT=.*#ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel082#" \
      -e "s#^RUN=.*#RUN=\"\$ROOT/$SC/boundfam-run-traced2.R\"#" \
      -e "s#^OUT=.*#OUT=\"\$ROOT/$O/suite2-\$arm\"#" \
      dev/boundfam-par2.sh > "$SC/boundfam-par2.sh"
  # the merged tests, in the lane's job-list format
  for f in tests/testthat/test-*.R; do
    echo "frmtmb /c/Users/adf44/source/r/frmtmb-wt-rel082/$f"
  done > "$SC/boundfam-jobs-core.txt"
  NOT_CRAN=true bash "$SC/boundfam-par2.sh" lane "$SC/boundfam-jobs-core.txt" 1 \
    > "$O/par2.log" 2>&1
  rr "$O/cmp-base.txt" dev/boundfam-cmp.R "$B/boundfam-suite-base" "$O/suite2-lane"
  rr "$O/cmp-lanefix.txt" dev/boundfam-cmp.R "$B/boundfam-suite3-lane" "$O/suite2-lane"
  rr "$O/classify.txt" dev/boundfam-classify.R "$O/suite2-lane" "$B/boundfam-suite-baseall"
}

t_gpstop() {
  local O=$KT/gpstop-out G=$SRC/frmtmb-wt-gpstop
  mkdir -p "$O"
  sed -e "s#^RV_WT <- .*#RV_WT <- \"$REL\"#" -e "s#^RV_LIB <- .*#RV_LIB <- \"$R21\"#" \
      -e "s#^RV_OUT <- .*#RV_OUT <- \"$REL/$O\"#" \
      dev/reviews/gpstop-review-lib.R > "$SC/gpstop-review-lib.R"
  sed -e "s#^GS_WT <- .*#GS_WT <- \"$REL\"#" -e "s#^GS_LIB <- .*#GS_LIB <- \"$R21\"#" \
      -e "s#^GS_OUT <- .*#GS_OUT <- \"$REL/$O\"#" \
      dev/gpstop-lib.R > "$SC/gpstop-lib.R"
  sed "s#$G/dev/reviews/gpstop-review-lib.R#$REL/$SC/gpstop-review-lib.R#" \
    dev/gpstop-fx2-sweep.R > "$SC/gpstop-fx2-sweep.R"
  sed "s#$G/dev/gpstop-lib.R#$REL/$SC/gpstop-lib.R#" \
    dev/gpstop-cost-body.R > "$SC/gpstop-cost-body.R"
  rr "$O/fx2-sweep-gauss.txt" "$SC/gpstop-fx2-sweep.R" gaussian 40,120 1,2 approx
  rr "$O/fx2-sweep-pb.txt" "$SC/gpstop-fx2-sweep.R" poisson,bernoulli 40,120 1,2 none
  rr "$O/fx2-sweep-student.txt" "$SC/gpstop-fx2-sweep.R" student 40,120 1,2 big
  rr "$O/fx2-sweep-sum.txt" dev/reviews/gpstop-review-sum.R \
    "$O/fx2-sweep-gauss.txt" "$O/fx2-sweep-pb.txt" "$O/fx2-sweep-student.txt"
  rr "$O/bcost-fresh.txt" "$SC/gpstop-cost-body.R" fresh
  rr "$O/bcost-lane96.txt" "$SC/gpstop-cost-body.R" lane96
}

# Not a lane table: the student rows of the gpstop sweep that moved on the
# merged build, fitted plainly on rellib-r20 and on each lane library, to
# find the lane that moves them. Its script is kt/scripts/gpstop-attrib.R.
t_gpstop_attrib() {
  local O=$KT/gpstop-out
  for l in base "$R21" $SRC/wt-*-lib; do
    rr "$O/attrib-$(basename "$l").txt" "$SC/gpstop-attrib.R" "$l"
  done
  grep -h "^ATTRIB" "$O"/attrib-*.txt > "$O/attrib-sum.txt"
}

all="quadseam boundfam-held ppdraw wcheck nlstop ordstart silent priortab
sampfront arspat boundfam-suite gpstop ordstart-sweep"
for t in ${@:-$all}; do
  echo "table $t $(date +%H:%M:%S)" >&2
  "t_${t//-/_}"
done
echo "done $(date +%H:%M:%S)" >&2
