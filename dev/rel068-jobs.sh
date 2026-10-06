#!/usr/bin/env bash
# Job lists for the 0.68.0 tiers (dev/rel068-log/*.jobs), from the
# release tree as it stands:
# every test file of all eight packages (ungated), run-gated.ps1's list
# (gated), and the scale files.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
pk="frmtmb frmtmb.eam frmtmb.sample frmtmb.coupling frmtmb.learn frmtmb.latent frmtmb.ode frmtmb.spline"
dir() { if [ "$1" = frmtmb ]; then echo "$ROOT/tests/testthat"; else echo "$ROOT/extensions/$1/tests/testthat"; fi; }
: > $ROOT/dev/rel068-log/suite.jobs; : > $ROOT/dev/rel068-log/scale.jobs
for p in $pk; do
  for f in $(ls "$(dir $p)"/test-*.R | sort); do echo "$p $f" >> $ROOT/dev/rel068-log/suite.jobs; done
  [ -f "$(dir $p)/test-scale.R" ] && echo "$p $(dir $p)/test-scale.R" >> $ROOT/dev/rel068-log/scale.jobs
done
c=$(dir frmtmb); sd=$(dir frmtmb.sample)
{
  for f in $(ls $c/test-bcm-*.R | sort); do echo "frmtmb $f"; done
  for nm in test-brms-agreement.R test-brms-likelihood.R test-brms-methods.R \
            test-brms-priors.R test-brms-port.R test-rl-example.R \
            test-drmtmb-agreement.R test-fuzz.R; do echo "frmtmb $c/$nm"; done
  for f in $(ls $c/test-brms-suite-*.R | sort); do echo "frmtmb $f"; done
  for f in $(ls $sd/test-brms-suite-*.R | sort); do echo "frmtmb.sample $f"; done
  echo "frmtmb.learn $(dir frmtmb.learn)/test-stan-identity.R"
  echo "frmtmb.sample $sd/test-loo.R"
  echo "frmtmb.sample $sd/test-sampling-ported.R"
} > $ROOT/dev/rel068-log/gated.jobs
while read -r p f; do [ -f "$f" ] || echo "MISSING $f"; done < $ROOT/dev/rel068-log/gated.jobs
wc -l $ROOT/dev/rel068-log/suite.jobs $ROOT/dev/rel068-log/gated.jobs $ROOT/dev/rel068-log/scale.jobs
