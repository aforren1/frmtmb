# log_lik:462 after lane catmulti's cumulative fix: brms's
# log_lik_cumulative() and frmtmb's ord_cumulative_logpmf() on brms's
# block (tests.log_lik.R:446, three thresholds all 0, logit), against
# the exact values: log F(disc * (0 - mu)) for y = 1,
# log(1 - F(disc * (0 - mu))) for y = 4, -Inf for the two zero-width
# categories. Seed 462 (the block draws its own; any seed shows it).
#   Rscript dev/rel072-cumtie-probe.R <lib>
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(brms)})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
set.seed(462)
ns <- 50; nobs <- 8; nthres <- 3; ncat <- 4
prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
prep$dpars <- list(mu = array(rnorm(ns * nobs), dim = c(ns, nobs)),
                   disc = rexp(ns), hu = rbeta(ns, 1, 1))
prep$thres$thres <- array(0, dim = c(ns, nthres))
prep$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
prep$family$link <- "logit"
lb <- sapply(1:nobs, brms:::log_lik_cumulative, prep = prep)
lk <- frmtmb:::ord_link("logit", "cumulative")
lf <- sapply(1:nobs, function(i) {
  vapply(1:ns, function(s) {
    frmtmb:::ord_cumulative_logpmf(prep$data$Y[i], prep$dpars$mu[s, i],
                                   prep$thres$thres[s, ], lk,
                                   prep$dpars$disc[s])
  }, 0)
})
ex <- sapply(1:nobs, function(i) {
  y <- prep$data$Y[i]
  e <- prep$dpars$disc * (0 - prep$dpars$mu[, i])
  if (y == 1) plogis(e, log.p = TRUE) else
    if (y == 4) plogis(e, lower.tail = FALSE, log.p = TRUE) else rep(-Inf, ns)
})
cat("brms: NaN in", sum(is.nan(lb)), "of", length(lb), "values\n")
cat("frmtmb: NaN in", sum(is.nan(lf)), "| -Inf in", sum(lf == -Inf),
    "(the zero-width categories:", sum(prep$data$Y %in% 2:3) * ns, ")\n")
fin <- is.finite(ex)
cat("frmtmb against the exact values, finite cells: max abs diff",
    format(max(abs(lf[fin] - ex[fin])), digits = 3), "| -Inf cells equal:",
    identical(lf[!fin], ex[!fin]), "\n")
