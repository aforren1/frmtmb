# dev/quadseam-acc.R bern seed 5: the lane is 9.0e-7 from glmer(nAGQ =
# 25) where base is 4.2e-8. At each build's estimate: the fit's
# log-likelihood, and an independent reference by stats::integrate()
# per group (plain R, no tape), so the rule's own error is read apart
# from where the optimizer stopped.
#   Rscript dev/quadseam-bern5.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(5)
d <- data.frame(g = factor(rep(1:100, each = 4)), x = rnorm(400))
d$y <- rbinom(400, 1, plogis(-0.5 + 0.7 * d$x + rnorm(100, 0, 1.2)[d$g]))
f <- suppressMessages(frm(y ~ x + (1 | g), family = bernoulli(), data = d,
                          quadrature = TRUE))
ref_ll <- function(b0, b1, s) {
  sum(vapply(split(d, d$g), function(dg) {
    h <- function(u) vapply(u, function(ui) {
      exp(sum(dbinom(dg$y, 1, plogis(b0 + b1 * dg$x + ui), log = TRUE)) +
            dnorm(ui, 0, s, log = TRUE))
    }, 0)
    log(integrate(h, -Inf, Inf, rel.tol = 1e-12)$value)
  }, 0))
}
b <- f$estimates$beta
s <- exp(f$estimates$theta[1])
l25 <- glmer(y ~ x + (1 | g), d, family = binomial, nAGQ = 25)
cat(sprintf("fit ll %.10f  ref at fit est %.10f  fit-ref %.3e\n",
            as.numeric(logLik(f)), ref_ll(b[1], b[2], s),
            as.numeric(logLik(f)) - ref_ll(b[1], b[2], s)))
gb <- fixef(l25)
gs <- attr(VarCorr(l25)$g, "stddev")
cat(sprintf("glmer25 ll %.10f  ref at glmer est %.10f  glmer-ref %.3e\n",
            as.numeric(logLik(l25)), ref_ll(gb[1], gb[2], gs),
            as.numeric(logLik(l25)) - ref_ll(gb[1], gb[2], gs)))
cat(sprintf("est fit %.6f %.6f %.6f  glmer %.6f %.6f %.6f\n", b[1], b[2], s,
            gb[1], gb[2], gs))
