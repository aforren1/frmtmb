# Reviewer, lane wt-priorform: where does frm_sample()'s own default stack
# put one slot twice? Resolves priors as frm_sample() does, without
# sampling. Lane build, seed 1.
#   Rscript dev/priorform-rev-sampledup.R
.libPaths(c("C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + dd$x + rnorm(10)[dd$g])
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
r <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, NULL))
cat("defaults only:\n"); print(as.data.frame(r$effective)[, 1:6])
pl <- set_prior("normal(0, 5)", class = "b")
fit2 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, prior = pl)
r2 <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit2, pl, base = fit2$prior))
cat("MAP fit prior passed again on the call:\n"); print(as.data.frame(r2$effective)[, 1:6])
