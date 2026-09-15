# Item 2.5d, the half that needs brms: does frmtmb's RESPONSE scale
# agree with brms's, and by how much does the DEFAULT of predict()
# disagree?
#
# The tolerance is never an absolute number.  The comparison is a
# ratio, and the yardstick is the run's own posterior uncertainty in
# the same quantity, so a wider posterior loosens the bound by itself.
#
#   Rscript dev/generics-scale-brms.R
LIB <- "C:/Users/adf44/source/r/generics-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("StanHeaders", as.character(packageVersion("StanHeaders")),
    " rstan", as.character(packageVersion("rstan")), "\n")
suppressMessages(library(brms))
suppressMessages(library(frmtmb))

set.seed(2026)
n <- 400
dd <- data.frame(x = rnorm(n), g = factor(rep(1:20, 20)))
eta <- 8 + 0.4 * dd$x + rnorm(20, 0, 0.3)[dd$g]
dd$y <- exp(rnorm(n, eta, 0.4))

fit <- frm(bf(y ~ x + (1 | g)) + lognormal(), data = dd)

cache <- Sys.getenv("FRMTMB_STAN_CACHE", "dev/stan-cache")
bf_file <- file.path(cache, "generics-scale-brmsfit.rds")
if (file.exists(bf_file)) {
  bfit <- readRDS(bf_file)
  cat("brms fit read from cache\n")
} else {
  bfit <- brm(y ~ x + (1 | g), family = brms::lognormal(), data = dd,
              chains = 2, iter = 2000, refresh = 0, seed = 2026,
              backend = "rstan")
  saveRDS(bfit, bf_file)
  cat("brms fit sampled and cached\n")
}

be <- fitted(bfit)                 # posterior mean of E[Y|x], response
bp <- predict(bfit)                # posterior predictive mean, response
fe <- fitted(fit)                  # frmtmb, response scale
fl <- predict(fit)                 # frmtmb DEFAULT, link scale

cat("\n== row 1, the four numbers ==\n")
cat(sprintf("frmtmb predict()   default (link)   %.6f\n", fl[1]))
cat(sprintf("frmtmb fitted()    (response)       %.3f\n", fe[1]))
cat(sprintf("brms   fitted()    Estimate         %.3f\n", be[1, 1]))
cat(sprintf("brms   predict()   Estimate         %.3f\n", bp[1, 1]))
cat(sprintf("brms predict() / exp(frmtmb link)   %.6f\n",
            bp[1, 1] / exp(fl[1])))
cat(sprintf("exp(sigma^2/2), sigma = %.6f       %.6f\n",
            sigma(fit), exp(sigma(fit)^2 / 2)))

# agreement on the response scale, as a RATIO to the posterior
# uncertainty the design itself reports
rel <- (fe - be[, 1]) / be[, 1]
mcse <- be[, 2] / be[, 1]           # posterior SD of E[Y|x], relative
cat("\n== agreement on the response scale, n = 400 rows ==\n")
cat(sprintf("max |frmtmb/brms - 1|                      %.6f\n",
            max(abs(rel))))
cat(sprintf("median relative posterior SD of brms epred %.6f\n",
            median(mcse)))
cat(sprintf("max |difference| / its own posterior SD    %.4f\n",
            max(abs(fe - be[, 1]) / be[, 2])))
cat(sprintf("median |difference| / its own posterior SD %.4f\n",
            median(abs(fe - be[, 1]) / be[, 2])))

cat("\n== the DEFAULT of predict() is a different scale ==\n")
cat(sprintf("max |frmtmb predict() / brms predict() - 1| %.6f\n",
            max(abs(fl / bp[, 1] - 1))))
cat(sprintf("ratio brms predict() to frmtmb predict()   %.1f\n",
            median(bp[, 1] / fl)))

cat("\n== sigma: link scale in summary(), response in sigma() ==\n")
cat(sprintf("frmtmb summary() sigma Estimate  %.6f (log scale)\n",
            summary(fit)$coefficients$sigma[1, 1]))
cat(sprintf("frmtmb sigma()                   %.6f\n", sigma(fit)))
bs <- posterior_summary(bfit)
cat(sprintf("brms   sigma posterior mean      %.6f\n",
            bs[grep("^sigma$", rownames(bs)), 1]))
cat("DONE\n")
