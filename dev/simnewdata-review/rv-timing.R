# Cost of simulate(newdata = <the fitted rows>) against simulate(),
# interleaved, minimum of 3 rounds, blocks over 1.2 s.
#   Rscript dev/simnewdata-review/rv-timing.R
source("dev/simnewdata-review/rv-prelude.R")
set.seed(7)
n <- 500
d <- data.frame(x = rnorm(n), g = factor(rep(1:25, 20)))
d$y <- 0.5 + 0.8 * d$x + rnorm(25)[d$g] + rnorm(n, 0, 0.5)
fit <- frm(bf(y ~ x + (1 | g), sigma ~ x), data = d)
tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; expr
  proc.time()[["elapsed"]] - t0 }
r <- replicate(3, c(
  fitted_rows = tm(simulate(fit, nsim = 2000, seed = 1)),
  newdata = tm(simulate(fit, nsim = 2000, seed = 1, newdata = d)),
  control = tm(simulate(fit, nsim = 2000, seed = 1))))
print(r)
cat("per-replicate ms (min of 3): fitted", 1000 * min(r[1, ]) / 2000,
    " newdata", 1000 * min(r[2, ]) / 2000, " control",
    1000 * min(r[3, ]) / 2000, "\n")
cat("pp_check default ndraws:", format(formals(frmtmb:::pp_check.frmtmb_fit)$ndraws), "\n")
