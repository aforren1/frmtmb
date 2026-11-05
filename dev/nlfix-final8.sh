#!/usr/bin/env bash
# Fix round 2, the runs on the final build, after the ord_fit_check() addition.
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RC="/c/Program Files/R/R-4.6.1/bin/R.exe"
LIB=C:/Users/adf44/source/r/wt-nlfix-lib2
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export NOT_CRAN=true
step() { echo "$(date +%H:%M:%S) $*" >> dev/nlfix-final8.progress; }
: > dev/nlfix-final8.progress
step install
bash dev/nlfix-install.sh > dev/nlfix-install-final.log 2>&1
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
for p in coupling eam latent learn ode sample spatial survival; do
  "$RC" CMD INSTALL --library="$LIB" --no-multiarch extensions/frmtmb.$p > dev/nlfix-install-$p.log 2>&1 || echo "INSTALL FAILED $p" >> dev/nlfix-final8.progress
done
unset R_LIBS
step "review reruns"
export NLFIX_LIB="$LIB;C:/Users/adf44/source/r/rellib-r16"
FRM_ARM=mutant REV_LIB=C:/Users/adf44/source/r/wt-nlfix-lib2 "$R" dev/nlfix-a8.R > dev/nlfix-a8.log 2>&1
cp dev/nlfix-a8-mutant.tsv dev/nlfix-a8-after.tsv
"$R" dev/nlfix-a8cmp.R > dev/nlfix-a8cmp.log 2>&1
FRM_ARM=mutant REV_LIB=C:/Users/adf44/source/r/wt-nlfix-lib2 "$R" dev/review-nlfix-a2.R > dev/nlfix-rev-a2.log 2>&1
FRM_ARM=mutant REV_LIB=C:/Users/adf44/source/r/wt-nlfix-lib2 "$R" dev/review-nlfix-a4c.R > dev/nlfix-rev-a4c.log 2>&1
"$R" dev/nlfix-cost2.R > dev/nlfix-cost2.log 2>&1
"$R" dev/nlfix-d1check.R > dev/nlfix-d1check.log 2>&1
"$R" dev/nlfix-startcheck.R > dev/nlfix-startcheck.log 2>&1
"$R" dev/nlfix-extfam.R > dev/nlfix-extfam.log 2>&1
"$R" dev/nlfix-asym.R > dev/nlfix-asym.log 2>&1
"$R" dev/nlfix-sar.R > dev/nlfix-sar.log 2>&1
"$R" dev/nlfix-arnl2.R > dev/nlfix-arnl2.log 2>&1
"$R" dev/nlfix-u955.R > dev/nlfix-u955.log 2>&1
"$R" dev/nlfix-r1.R > dev/nlfix-r1.log 2>&1
"$R" dev/nlfix-r2.R > dev/nlfix-r2.log 2>&1
"$R" dev/nlfix-acat.R > dev/nlfix-acat.log 2>&1
"$R" dev/nlfix-lgscan.R C:/Users/adf44/source/r/wt-nlfix-lib2 > dev/nlfix-lgscan-after.log 2>&1
"$R" dev/nlfix-cost3.R > dev/nlfix-cost3b.log 2>&1
"$R" dev/nlfix-cost4.R > dev/nlfix-cost4.log 2>&1
FRM_ARM=mutant REV_LIB=C:/Users/adf44/source/r/wt-nlfix-lib2 "$R" dev/review-nlfix-c4.R > dev/nlfix-rc-c4-after.log 2>&1
FRM_ARM=mutant REV_LIB=C:/Users/adf44/source/r/wt-nlfix-lib2 "$R" dev/review-nlfix-c1b.R > dev/nlfix-rc-c1b-after.log 2>&1
FRM_ARM=mutant REV_LIB=C:/Users/adf44/source/r/wt-nlfix-lib2 "$R" dev/review-nlfix-c3b.R > dev/nlfix-rc-c3b-after.log 2>&1
unset NLFIX_LIB
step "saddle table"
mv dev/nlfix-saddle-out dev/nlfix-saddle-r6-out 2>/dev/null
mv dev/nlfix-saddle.done dev/nlfix-saddle-r6.done 2>/dev/null
bash dev/nlfix-saddle-run.sh > /dev/null 2>&1
step "L1 default-start specs"
bash dev/nlfix-l1-run.sh dev/nlfix-l3-out sad_cumulative sad_sratio sad_cratio sad_acat sad_hurdle_cumulative sad_probit sad_sat sad_gaussian sad_poisson powzero powzero_fun > dev/nlfix-l4.log 2>&1
step "core ungated"
mv dev/nlfix-suite-core dev/nlfix-suite-core-r6 2>/dev/null
bash dev/nlfix-suite.sh core 0 8 frmtmb > dev/nlfix-suite-core.log 2>&1
step "extensions ungated"
mv dev/nlfix-suite-ext dev/nlfix-suite-ext-r6 2>/dev/null
bash dev/nlfix-suite.sh ext 0 8 frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spatial frmtmb.survival > dev/nlfix-suite-ext.log 2>&1
step "gated"
gated=$(grep -l "SKIP=[1-9]" dev/nlfix-suite-core/*.txt | sed -E 's/.*frmtmb--(.*)[.]txt/\1/' | paste -sd'|')
mv dev/nlfix-suite-gated dev/nlfix-suite-gated-r6 2>/dev/null
ONLY="($gated)[.]R" bash dev/nlfix-suite.sh gated 1 8 frmtmb > dev/nlfix-suite-gated.log 2>&1
gateds=$(grep -l "SKIP=[1-9]" dev/nlfix-suite-ext/frmtmb.sample--*.txt | sed -E 's/.*frmtmb.sample--(.*)[.]txt/\1/' | paste -sd'|')
mv dev/nlfix-suite-gatedsample dev/nlfix-suite-gatedsample-r6 2>/dev/null
ONLY="($gateds)[.]R" bash dev/nlfix-suite.sh gatedsample 1 8 frmtmb.sample > dev/nlfix-suite-gatedsample.log 2>&1
step "sweep and equiv"
mv dev/nlfix-equiv-out dev/nlfix-equiv-r6-out 2>/dev/null
mkdir -p dev/nlfix-equiv-out
FRMTMB_EQUIV_LIB="$LIB;C:/Users/adf44/source/r/rellib-r16" "$R" dev/equiv-nl-sweep.R dev/nlfix-equiv-out > dev/nlfix-sweep.log 2>&1
base=/c/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out
cp "$base"/p*.tsv dev/nlfix-equiv-out/
cp "$base"/p*.rds dev/nlfix-equiv-out/ 2>/dev/null
for id in $(tr -d '\r' < dev/nlfix-equiv-ids.txt); do
  rm -f "dev/nlfix-equiv-out/$id.tsv" "dev/nlfix-equiv-out/$id.rds"
done
R_PARALLEL_PORT=random FRMTMB_EQUIV_LIB="$LIB;C:/Users/adf44/source/r/rellib-r16" bash dev/equiv-run.sh 8 dev/nlfix-equiv-out > dev/nlfix-equiv-run.log 2>&1
"$R" dev/nlfix-equiv-compare.R > dev/nlfix-equiv-compare.log 2>&1
step sessions
export FRMTMB_RELLIB="$LIB"
bash dev/release/nlfix-run-sessions.sh nlfixsess3 8 > dev/nlfix-sessions3.out 2>&1
cp dev/release/nlfixsess3.log dev/nlfix-sessions3.log
unset FRMTMB_RELLIB
step openblas
mkdir -p dev/nlfix-log2
for ver in 0.3.26 0.3.32; do
  RB="dev/ciharden-out/Rob$ver-lapack/bin/x64/Rscript.exe"
  for f in tests/testthat/test-nl-family.R tests/testthat/test-nl-stan-funs.R tests/testthat/test-beta-binomial-density.R tests/testthat/test-cdfs.R tests/testthat/test-cens-trunc.R tests/testthat/test-numerical-robustness.R; do
    R_PARALLEL_PORT=12050 OPENBLAS_NUM_THREADS=4 "$RB" dev/nlfix-testrun.R frmtmb "$f" dev/nlfix-log2/ob$ver-$(basename "$f" .R).txt > /dev/null 2>&1
  done
  FRMTMB_BRMS_FIT_TESTS=true R_PARALLEL_PORT=12050 OPENBLAS_NUM_THREADS=4 "$RB" dev/nlfix-testrun.R frmtmb tests/testthat/test-nl-brms.R dev/nlfix-log2/ob$ver-test-nl-brms.txt > /dev/null 2>&1
  R_PARALLEL_PORT=12050 OPENBLAS_NUM_THREADS=4 "$RB" dev/nlfix-acat.R > dev/nlfix-acat-ob$ver.log 2>&1
done
step "R CMD check"
for p in extensions/frmtmb.sample extensions/frmtmb.eam extensions/frmtmb.learn extensions/frmtmb.survival; do
  n=$(basename "$p")
  R_PARALLEL_PORT=$((12010 + ${#n})) bash dev/nlfix-check2.sh "$p" > dev/nlfix-check2-$n.out 2>&1 &
done
wait
# core alone: its R CMD build copies the whole worktree, dev/nlfix-check included
R_PARALLEL_PORT=12030 bash dev/nlfix-check2.sh . > dev/nlfix-check2-frmtmb.out 2>&1
step done
echo FINAL8 DONE > dev/nlfix-final8.done
