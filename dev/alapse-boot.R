# Lane alapse: lane squant's dev/squant-boot-ald.R on this lane's build
# (frm_bootstrap() on asym_laplace(), gaussian data, 400 rows, quantile
# 0.9, 200 draws, seeds 20261117 to 20261119). squant's Monte Carlo sd
# of rq()'s slope on this design is 0.2977.
#   Rscript dev/alapse-boot.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/alapse-common.R")
for (s in 1:3) {
  set.seed(20261116 + s)
  n <- 400
  dd <- data.frame(x = runif(n))
  dd$y <- 1 + 2 * dd$x + rnorm(n)
  fit <- suppressWarnings(frm(bf(y ~ x, quantile = 0.9),
                              family = asym_laplace(), data = dd))
  b <- suppressWarnings(frm_bootstrap(fit, nsim = 200, seed = s))
  r <- quantreg::rq(y ~ x, tau = 0.9, data = dd)
  cat(sprintf(paste0("seed %d: slope bootstrap SE %.4f (converged %d of",
                     " 200) | rq nid SE %.4f\n"), s, sd(b$t[, 2]),
              sum(b$converged), summary(r, se = "nid")$coefficients[2, 2]))
}
