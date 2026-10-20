source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
set.seed(1)
n <- 40
d <- data.frame(x = rnorm(n, 3, 1), g = factor(rep(1:5, length.out = n)))
d$y <- 20 + 5 * d$x + rnorm(n, 0, 10)
od <- cap_obj(y ~ x, d, gaussian())
of <- cap_obj(y ~ x, d, gaussian(), prior = "flat")
print(od$env$par)
print(of$env$par)
print(od$env$random)
p <- od$env$par
cat("lp at par:", logprior_at(od, of, p), "\n")
cat(brms::stancode(y ~ x, data = d, family = gaussian()), "\n")
print(brms::default_prior(y ~ x, data = d, family = gaussian()))
