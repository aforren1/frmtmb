#!/usr/bin/env bash
# Each mutant against the lane's new test files, one R process per file, 6 at once.
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
out=dev/review-nlfix-mut-out
mkdir -p "$out"
T=tests/testthat
jobs=()
for m in none is_primary guard_off guard_obj_off guard_design_off domain_off d_Phi d_cloglog d_lse d_lde d_hypot d_invlogit plain_l1pe step_zero kept_off thres_rows ac_drop ac_refit mm_first shadow_off; do
  jobs+=("$m $T/test-nl-family.R frmtmb false")
  jobs+=("$m $T/test-nl-stan-funs.R frmtmb false")
done
for m in none is_primary thres_rows ac_drop mm_first; do
  jobs+=("$m $T/test-nl-brms.R frmtmb true")
done
for m in none d_Phi d_cloglog d_lse d_lde; do
  jobs+=("$m $T/test-nl-stan-funs.R frmtmb true")
done
for m in none is_primary ac_drop; do
  jobs+=("$m extensions/frmtmb.sample/tests/testthat/test-nl-body-draws.R frmtmb.sample false")
done
printf '%s\n' "${jobs[@]}" | xargs -P 6 -L 1 bash -c '
  m=$0; f=$1; p=$2; g=$3
  tag=$m-$(basename $f .R)-$g
  if [ "$g" = true ]; then export FRMTMB_BRMS_FIT_TESTS=true; fi
  timeout 1500 "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/review-nlfix-mutrun.R $m $f $p > dev/review-nlfix-mut-out/$tag.log 2>&1 || echo "EXIT nonzero $tag" >> dev/review-nlfix-mut-out/$tag.log
'
cat "$out"/*.log | grep -E '^MUT|BAD|EXIT' > dev/review-nlfix-mut.log
