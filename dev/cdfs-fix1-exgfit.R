# Fix round 1, finding 5: how far RTMBdist::dexgauss()'s approximation
# from sigma / beta of about 20 moves frmtmb's own exgaussian
# likelihood. On the given library: the objective at fixed parameters
# beside the closed form (brms::dexgaussian()), and fits whose maximum
# lies at a large sigma / beta, re-scored with the closed form at their
# own estimates.
#   Rscript dev/cdfs-fix1-exgfit.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", find.package("frmtmb"), "\n")
# gaussian data, the ex-gaussian's beta -> 0 limit
set.seed(20261022)
d <- data.frame(x = stats::rnorm(400))
d$y <- 1 + 0.5 * d$x + stats::rnorm(400)
ob <- frm(y ~ x, family = exgaussian(), data = d, dry_run = "objective")$obj
cat("== the objective at mu = 1 + 0.5 x, sigma = 1, beta = 1 / u\n")
for (u in c(10, 20, 50, 100, 1000)) {
  p <- c(1, 0.5, log(1), log(1 / u))
  mu <- 1 + 0.5 * d$x
  closed <- sum(brms::dexgaussian(d$y, mu, 1, 1 / u, log = TRUE))
  cat(sprintf("u %5g: -objective %.6f closed form %.6f diff %.4g\n", u,
              -ob$fn(p), closed, -ob$fn(p) - closed))
}
cat("== fits, re-scored with the closed form\n")
for (s in 1:4) {
  set.seed(20261022 + s)
  d <- data.frame(x = stats::rnorm(400))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(400) + stats::rexp(400, 30)
  fit <- tryCatch(frm(y ~ x, family = exgaussian(), data = d),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("seed %d: ERROR %s\n", s, substr(conditionMessage(fit), 1, 80)))
    next
  }
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  closed <- sum(brms::dexgaussian(d$y, dp$mu, dp$sigma, dp$beta, log = TRUE))
  cat(sprintf(paste("seed %d: logLik %.6f closed form %.6f diff %.3g |",
                    "u %.4g | b %s | conv %d\n"),
              s, as.numeric(logLik(fit)), closed,
              as.numeric(logLik(fit)) - closed, mean(dp$sigma / dp$beta),
              paste(sprintf("%.5f", fixef_by_dpar(fit)$mu), collapse = " "),
              fit$opt$convergence))
}
cat("END\n")
