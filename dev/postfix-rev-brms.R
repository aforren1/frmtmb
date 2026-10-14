# Reviewer of lane postfix: brms 2.23.0's conditional_effects(method =
# "posterior_predict") on a poisson and a bernoulli model, with and
# without `transform`, for comparison with frmtmb.sample's draws method
# (dev/postfix-rev-sample.R reads the rds this writes).
# Usage: Rscript dev/postfix-rev-brms.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
source("dev/postfix-rev-data.R")
ic <- list(x = c(-1, -0.5, 0, 0.5, 1))
q <- function(e) suppressWarnings(suppressMessages(e))
bp <- q(brm(y ~ x, family = poisson(), data = dp, chains = 4, iter = 2000,
            seed = 1, refresh = 0))
bb <- q(brm(y ~ x, family = bernoulli(), data = db, chains = 4, iter = 2000,
            seed = 1, refresh = 0))
set.seed(2)
out <- list(
  pois = conditional_effects(bp, "x", method = "posterior_predict",
                             int_conditions = ic)$x,
  pois_sqrt = conditional_effects(bp, "x", method = "posterior_predict",
                                  int_conditions = ic, transform = sqrt)$x,
  pois_epred = conditional_effects(bp, "x", int_conditions = ic)$x,
  bern = conditional_effects(bb, "x", method = "posterior_predict",
                             int_conditions = ic)$x,
  fixef_pois = fixef(bp), fixef_bern = fixef(bb))
for (k in names(out)[1:4]) {
  cat("--", k, "\n")
  print(out[[k]][, c("x", "estimate__", "se__", "lower__", "upper__")])
}
print(out$fixef_pois); print(out$fixef_bern)
## the convention, read off the draws
pp <- posterior_predict(bp, newdata = data.frame(x = ic$x))
cat("median", apply(pp, 2, median), "\n")
cat("mad   ", apply(pp, 2, mad), "\n")
cat("q2.5  ", apply(pp, 2, quantile, 0.025), "\n")
saveRDS(out, "dev/postfix-rev-log/brms-ce.rds")
