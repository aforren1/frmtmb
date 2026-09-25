# Punch round 1, item 4: the cost of simulate(newdata = ) per replicate,
# against simulate() on the fitted rows as the CONTROL, interleaved in
# one process, minimum of 5 rounds, each block 2000 replicates (over
# 1.2 s). Reviewer's design (dev/simnewdata-review/rv-timing.R).
#   Rscript dev/simnewdata-timing.R <label>
lab <- commandArgs(trailingOnly = TRUE)[1]
source("dev/simnewdata-prelude.R")
suppressMessages(library(frmtmb))
set.seed(7)
n <- 500
d <- data.frame(x = rnorm(n), g = factor(rep(1:25, 20)))
d$y <- 0.5 + 0.8 * d$x + rnorm(25)[d$g] + rnorm(n, 0, 0.5)
fit <- frm(bf(y ~ x + (1 | g), sigma ~ x), data = d)
tm <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
r <- replicate(5, c(
  fitted_rows = tm(simulate(fit, nsim = 2000, seed = 1)),
  newdata = tm(simulate(fit, nsim = 2000, seed = 1, newdata = d)),
  newdata_NA = tm(simulate(fit, nsim = 2000, seed = 1, newdata = d,
                           re_formula = NA))))
print(r)
m <- apply(r, 1, min)
cat(sprintf("%s: per-replicate ms, min of 5: fitted rows %.3f, newdata %.3f, newdata NA %.3f; newdata / fitted %.2f\n",
            lab, 1000 * m[1] / 2000, 1000 * m[2] / 2000, 1000 * m[3] / 2000,
            m[2] / m[1]))
