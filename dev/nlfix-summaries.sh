#!/usr/bin/env bash
# The count blocks of dev/nlfix-findings.md, from the logs. Run after both final
# chains, then Rscript dev/nlfix-report.R.
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
{
  echo "Before, the new test files on rellib-r16 (0.78.0):"
  for f in family stanfuns nlbrms sample; do
    sed -n 3p dev/nlfix-suite-r16/$f.txt
  done
  echo "Behavioral failures among them, first lines:"
  grep -hE "reaches the likelihood of cumulative|inv_logit|sampled\]" \
    dev/nlfix-suite-r16/family.txt dev/nlfix-suite-r16/stanfuns.txt \
    dev/nlfix-suite-r16/sample.txt | cut -c1-150 | head -6
  grep -A2 "reaches the likelihood of cumulative" dev/nlfix-suite-r16/family.txt | cut -c1-150 | head -4
} > dev/nlfix-summary-r16.txt
{
  echo "Core, ungated, one file per process, lane library:"
  tail -3 dev/nlfix-suite-core.log
  grep -h "^FILE" dev/nlfix-suite-core/frmtmb--test-nl*.txt
} > dev/nlfix-summary-core.txt
{
  echo "Extensions, ungated, one file per process, lane core, base extensions:"
  tail -3 dev/nlfix-suite-ext.log
  for p in coupling eam latent learn ode sample spatial survival; do
    grep -h "^FILE" dev/nlfix-suite-ext/frmtmb.$p--*.txt |
      awk -v p=$p '{ n++; split($3,a,"="); split($4,b,"="); split($5,c,"="); split($6,d,"="); pa+=a[2]; fa+=b[2]; er+=c[2]; sk+=d[2] } END { printf "frmtmb.%s files=%d PASS=%d FAIL=%d ERROR=%d SKIP=%d\n", p, n, pa, fa, er, sk }'
  done
  grep -h "^FILE" dev/nlfix-suite-ext/frmtmb.sample--test-nl-body-draws.txt
} > dev/nlfix-summary-ext.txt
{
  echo "Gated (FRMTMB_BRMS_FIT_TESTS=true), the files that skipped ungated:"
  tail -3 dev/nlfix-suite-gated.log
  grep -h "^FILE" dev/nlfix-suite-gated/*.txt | grep -E "SKIP=[1-9]"
  grep -h "^FILE" dev/nlfix-suite-gated/frmtmb--test-nl-brms.txt
  grep -h "^FILE" dev/nlfix-suite-gated/frmtmb--test-nl-stan-funs.txt
  if [ -f dev/nlfix-suite-gatedsample.log ]; then
    echo "frmtmb.sample, gated:"
    tail -3 dev/nlfix-suite-gatedsample.log
  fi
  echo "OpenBLAS 0.3.26 and 0.3.32 with LAPACK, new test files (gated for test-nl-brms):"
  grep -H "^FILE" dev/nlfix-log2/ob0*.txt 2>/dev/null | sed 's|dev/nlfix-log2/||'
} > dev/nlfix-summary-gated.txt
{
  echo "One session, forward, reverse and preload:"
  grep -aE "^SESSION|^LEAK|SESS|ran [0-9]+ of" dev/nlfix-sessions3.log | cut -c1-170
} > dev/nlfix-summary-sess.txt
{
  for p in frmtmb frmtmb.sample frmtmb.eam frmtmb.learn frmtmb.survival; do
    echo "$p: $(grep -a '^Status' dev/nlfix-check2/$p/check.log)"
    grep -aE "^\* checking .*(WARNING|NOTE|ERROR)$" dev/nlfix-check2/$p/check.log | cut -c1-120
  done
} > dev/nlfix-summary-check.txt
