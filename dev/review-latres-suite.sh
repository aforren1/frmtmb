#!/usr/bin/env bash
# The whole gated run on the final build, which the lane did not repeat. One R
# process per test file, N at a time, on the lane's library, logs under
# dev/review-latres-out/gated. Extensions from the lane's tree.
#   bash dev/review-latres-suite.sh <N>
source "$(dirname "$0")/latres-env.sh"
cd "$LATRES_ROOT" || exit 1
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/stan-cache
export FRMTMB_BRMS_FIT_TESTS=true
N="$1"
out=dev/review-latres-out/gated
mkdir -p "$out"
for p in core frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spatial frmtmb.survival; do
  if [ "$p" = core ]; then
    dir=tests/testthat; pk=frmtmb
  else
    dir="extensions/$p/tests/testthat"; pk="$p"
  fi
  for f in "$dir"/test-*.R; do
    b="$(basename "$f" .R)"
    [ -s "$out/$pk--$b.txt" ] && grep -q "^RESULT" "$out/$pk--$b.txt" && continue
    while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
    ( "$RSCRIPT" dev/latres-test1.R "$f" "$pk" > "$out/$pk--$b.txt" 2>&1 ) &
  done
done
wait
grep -h "^RESULT" "$out"/*.txt > dev/review-latres-out/gated.summary
echo "ran $(wc -l < dev/review-latres-out/gated.summary) of $(ls "$out"/*.txt | wc -l) files"
awk '{split($2,pk," "); for(i=4;i<=NF;i++){split($i,kv,"="); s[$2,kv[1]]+=kv[2]; k[$2]=1}}
     END{for(p in k){printf "%s tests=%d pass=%d fail=%d err=%d skip=%d warn=%d\n",p,s[p,"tests"],s[p,"pass"],s[p,"fail"],s[p,"err"],s[p,"skip"],s[p,"warn"]}}' \
  dev/review-latres-out/gated.summary | sort
awk '{for(i=4;i<=NF;i++){split($i,kv,"="); if ((kv[1]=="fail"||kv[1]=="err") && kv[2]>0) {print; next}}}' dev/review-latres-out/gated.summary
grep -L "^lib: C:/Users/adf44/source/r/wt-latres-lib" "$out"/*.txt | head
