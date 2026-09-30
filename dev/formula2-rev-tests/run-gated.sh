#!/bin/bash
# Reviewer driver: the gated files, both arms, counts; then the brms
# suite files once more in record mode for each row's message.
cd /c/Users/adf44/source/r/frmtmb-wt-formula2
O=dev/formula2-rev-tests
{
for arm in before after; do
  for f in test-brms-suite-brmsformula.R test-brms-suite-standata.R \
           test-brms-suite-priors.R test-brms-suite-methods.R \
           test-brms-likelihood.R test-dpar-equate.R test-family-list.R \
           test-formula-cmc.R test-formula-dot.R; do
    echo "$arm frmtmb $f $O/g-$arm-$f.log"
  done
  echo "$arm frmtmb.sample test-brms-suite-brmsformula.R $O/g-$arm-sample-test-brms-suite-brmsformula.R.log"
done
} | xargs -P 20 -n 4 sh -c \
  '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/formula2-rev-runtest-gated.R "$0" "$1" "$2" TRUE > "$3" 2>&1'
for arm in before after; do
  for f in test-brms-suite-brmsformula.R test-brms-suite-standata.R \
           test-brms-suite-priors.R test-brms-suite-methods.R; do
    echo "$arm $f"
  done
done | xargs -P 8 -n 2 sh -c \
  'rm -f "dev/formula2-rev-tests/rec-$0-$1.tsv"; "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/formula2-rev-runtest-gated.R "$0" frmtmb "$1" TRUE "C:/Users/adf44/source/r/frmtmb-wt-formula2/dev/formula2-rev-tests/rec-$0-$1.tsv" > "dev/formula2-rev-tests/rec-$0-$1.log" 2>&1'
echo ALLDONE
