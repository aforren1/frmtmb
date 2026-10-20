source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
cat(brms::stancode(y | thres(gr = gg) ~ x, data = D$ord, family = brms::cumulative()))
print(brms::default_prior(y | thres(gr = gg) ~ x, data = D$ord, family = brms::cumulative()))
cat(brms::stancode(y ~ x, data = D$ord, family = brms::cumulative()))
