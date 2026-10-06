# Reviewer: fit the n = 140, seed 11 data set of dev/cifix-scan2.R
# several times in one process, with and without a gc() before each
# fit, and print theta_1 to 17 digits.
# Usage: Rscript dev/cifixrev-determ4.R <lib or "base"> <gc|nogc|torture>
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
n <- 140L
set.seed(11)
d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
  stats::rnorm(n, 0, 0.3)
for (k in 1:6) {
  if (args[2] == "gc") invisible(gc())
  fit <- suppressWarnings(
    frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)), data = d))
  cat("FIT", k, sprintf("%.17g", fit$opt$par[6]), fit$opt$evals,
      fit$opt$iterations, "\n")
}
