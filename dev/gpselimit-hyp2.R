# Lane gpselimit: is the product hypothesis() reports for a white field the
# one the fit determines? The two-dimensional seed-1 fields of
# gpselimit-hyp.R under exp_quad and matern32: the common prior sd of the
# basis coefficients from the fit (hsgp_sds() at the estimate), the
# hypothesis estimate, and 2 log sdgp + D log lscale from theta and from
# the brms values.
# Usage: Rscript gpselimit-hyp2.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-lib.R")
a <- commandArgs(TRUE)
gs_libs(a[1])
ns <- asNamespace("frmtmb")
set.seed(5150)
n <- 150
d2 <- data.frame(x1 = runif(n), x2 = runif(n), y = rnorm(n))
for (cv in c("exp_quad", "matern32")) {
  form <- bf(y ~ gp(x1, x2, k = 4, iso = TRUE, cov = cv))
  f0 <- suppressWarnings(suppressMessages(frm(form, data = d2,
                                              family = gaussian())))
  bk <- Filter(function(b) b$covstruct == "hsgp", f0$frame$re_blocks)[[1]]
  lp <- f0$frame$linpreds[[ns$gp_block_lp(bk, f0$frame)]]
  Z <- as.matrix(ns$lp_Z(f0, lp))[, bk$c_idx, drop = FALSE]
  set.seed(5161)
  dd <- d2
  dd$y <- as.vector(Z %*% rnorm(ncol(Z))) + 0.1 * rnorm(n)
  fit <- suppressMessages(frm(form, data = dd, family = gaussian()))
  bk <- Filter(function(b) b$covstruct == "hsgp", fit$frame$re_blocks)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  sds <- ns$hsgp_sds(th, bk$aux_omega, bk$gp_iso, ns$gp_cov_of(bk))
  bv <- ns$gp_brms_values(fit, fit$estimates$theta)
  h <- hypothesis(fit, "2 * log(sdgp_gpx1x2) + 2 * log(lscale_gpx1x2) = 0",
                  class = "")$hypothesis
  cat(sprintf(paste0("HYP2 %s | coef prior sd %s | 2th1+2th2 %.6f | ",
                     "brms values %s | hyp %.6f se %.6f | ll %.6f\n"),
              cv, paste(format(range(sds), digits = 6), collapse = ".."),
              2 * th[1] + 2 * th[2],
              paste(names(bv), format(bv, digits = 6), collapse = " "),
              h$Estimate, h$Est.Error, as.numeric(logLik(fit))))
}
