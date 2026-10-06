# The simultaneous critical value over an exact gp(), inside, at the edge
# of and past the observed positions: the construction of the R9 table
# of dev/reviews/2026-09-08-diffcurve.md (y ~ fac + gp(x), 60
# observations on [0, 6], noise 0.2, nsim = 20000, seed 1), rebuilt here
# with its data generator written down, since that review did not keep
# one. Run on each arm: Rscript dev/gpby-crit.R base|lane
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spline"), "\n")
set.seed(1)
n <- 60
d <- data.frame(x = stats::runif(n, 0, 6),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(d$x) + ifelse(d$fac == "B", 0.5, 0) + stats::rnorm(n, 0, 0.2)
fit <- frm(bf(y ~ fac + gp(x)), data = d)
# the max-deviation quantile of a draw from S standardized by `div`
crit_of <- function(S, div, nsim = 20000, seed = 1) {
  set.seed(seed)
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)
  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow(S))
  z <- matrix(stats::rnorm(nrow(S) * nsim), nrow(S), nsim)
  unname(stats::quantile(apply(abs((L %*% z) / div), 2, max), 0.95,
                         type = 8))
}
for (g in list(c(1, 5), c(5.5, 6.5), c(7, 12))) {
  nd <- data.frame(x = seq(g[1], g[2], length.out = 25),
                   fac = factor("A", levels = levels(d$fac)))
  cv <- frm_curve(fit, newdata = nd, nsim = 20000, seed = 1)
  lb <- frm_lp_basis(fit, newdata = nd, re_formula = NA)
  A <- as.matrix(lb$A)
  S0 <- A %*% lb$V %*% t(A)
  r <- lb$extra_var / (diag(S0) + lb$extra_var)
  # the reviewer's bound: the same draw standardized by its own scale
  rescaled <- crit_of(S0, sqrt(diag(S0)))
  cat(sprintf(paste0("GRID [%.1f, %.1f] r in [%.5f, %.5f] | frm_curve ",
                     "crit %.5f | A V A' self-standardized %.5f | ratio ",
                     "%.5f\n"),
              g[1], g[2], min(r), max(r), cv$.crit_sim[1], rescaled,
              cv$.crit_sim[1] / rescaled))
}
