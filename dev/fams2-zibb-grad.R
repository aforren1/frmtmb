# frmtmb's own gradient at the zibb row's optimum, by seed
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
for (s in 16:24) {
  set.seed(s)
  n <- 300
  db <- data.frame(x = rnorm(n), tr = sample(5:15, n, TRUE))
  mu <- plogis(-0.3 + 0.5 * db$x)
  yb <- stats::rbinom(n, db$tr, stats::rbeta(n, mu * 4, (1 - mu) * 4))
  db$y <- ifelse(runif(n) < plogis(-1 + 0.4 * db$x), 0L, yb)
  fit <- frm(bf(y | trials(tr) ~ x, zi ~ x), family = zero_inflated_beta_binomial(), data = db)
  g <- fit$obj$gr(fit$obj$env$last.par.best)
  cat("seed", s, "max|gr|", format(max(abs(g)), digits = 3), " conv", fit$opt$convergence %||% NA, "\n")
}
