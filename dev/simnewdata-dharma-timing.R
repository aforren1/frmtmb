# The examples-timing NOTE of the punch-round check names
# dharma_residuals (12.25 s elapsed). Is it the lane or the load? The
# example's simulate() call (the fitted rows, re_formula = NULL, 250
# replicates), base against lane in two processes each, interleaved,
# with a fixed arithmetic CONTROL in each process, minimum of 5.
#   SIMNEWDATA_LIB=base Rscript dev/simnewdata-dharma-timing.R base
#   Rscript dev/simnewdata-dharma-timing.R lane
lab <- commandArgs(trailingOnly = TRUE)[1]
source("dev/simnewdata-prelude.R")
if (identical(lab, "base")) .libPaths(.libPaths()[-1L])
suppressMessages(library(frmtmb))
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
tm <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
r <- replicate(5, c(
  sim = tm(for (i in 1:8) simulate(fit, nsim = 250, seed = 1)),
  control = tm(for (i in 1:40) sum(sort(runif(1e5))))))
cat(sprintf("%s: min of 5, 8 x simulate(nsim = 250) %.2f s, control %.2f s\n",
            lab, min(r[1, ]), min(r[2, ])))
