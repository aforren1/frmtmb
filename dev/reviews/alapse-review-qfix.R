# Review of lane alapse, check 1: confirm with frm() itself that the
# higher maxima of alapse-review-qglob.R / -multistart.R exist: refit
# with the quantile fixed at the better value and compare logLik.
#   Rscript dev/reviews/alapse-review-qfix.R lane
lib_arg <- "lane"
source("dev/reviews/alapse-review-common.R")
mkB <- function(s) {
  set.seed(2000 + s)
  d <- data.frame(x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(100)
  d
}
for (cs in list(c(109, 0.7247), c(22, 0.1487), c(140, 0.3565))) {
  d <- mkB(cs[1])
  r0 <- cap(frm(y ~ x, family = asym_laplace(), data = d))
  r1 <- cap(frm(bf(y ~ x, quantile = cs[2]), family = asym_laplace(),
                data = d))
  q0 <- ns$eval_dpars(r0$v)[[1]]$quantile[1]
  cat(sprintf("B %d: estimated quantile %.4f logLik %.6f code %d [%s] | fixed at %.4f logLik %.6f | slope %.4f vs %.4f\n",
              cs[1], q0, as.numeric(logLik(r0$v)), r0$v$opt$convergence,
              paste(r0$wc, collapse = ","), cs[2],
              as.numeric(logLik(r1$v)),
              ns$fixef_by_dpar(r0$v)$mu[2], ns$fixef_by_dpar(r1$v)$mu[2]))
  # and from a start near the better quantile
  r2 <- cap(frm(y ~ x, family = asym_laplace(), data = d,
                init_dpars = list(quantile = cs[2])))
  if (!inherits(r2$v, "error")) {
    cat(sprintf("   started at quantile %.4f: estimated %.4f logLik %.6f code %d cert %s [%s]\n",
                cs[2], ns$eval_dpars(r2$v)[[1]]$quantile[1],
                as.numeric(logLik(r2$v)), r2$v$opt$convergence,
                isTRUE(r2$v$opt$ald_cert$ok), paste(r2$wc, collapse = ",")))
  } else cat("   init_dpars error:", conditionMessage(r2$v), "\n")
}
cat("DONE\n")
