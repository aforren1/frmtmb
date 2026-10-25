# Review: two priors on one ar() parameter (class ar and theta/thetaac_1)
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(x = rnorm(60), t = rep(1:20, 3), s = rep(1:3, each = 20))
d$y <- d$x + rnorm(60)
f <- frm(y ~ x + ar(t, s), data = d, prior = c(
  set_prior("normal(0, 0.5)", class = "ar"),
  set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1")))
cat(format(packageVersion("frmtmb")), "\n"); print(prior_summary(f))
