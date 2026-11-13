# Review of lane sampfront: does the GEV formula route sample the same
# posterior as the fit route (backlog 3215)?
#   Rscript dev/review-sampfront-gev.R
# Data seed 4 (as the lane's repro), sampler seed 1, 2 chains, 1000 iter.
.libPaths(c("C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb.sample", find.package("frmtmb.sample"), "\n")
set.seed(4)
rgev <- function(n, mu, s, xi) mu + s * ((-log(runif(n)))^(-xi) - 1) / xi
n <- 150
dv <- data.frame(x = rnorm(n), g = factor(rep(1:5, length.out = n)))
dv$y <- rgev(n, 0.5 * dv$x + c(-0.4, 0, 0.3, 0.1, -0.2)[dv$g], 1, 0.15)
smp <- function(...) {
  suppressWarnings(suppressMessages(
    frm_sample(..., chains = 2, iter = 1000, refresh = 0, seed = 1)))
}
summ <- function(s, vars) {
  m <- as.matrix(s)[, vars, drop = FALSE]
  rbind(mean = colMeans(m), sd = apply(m, 2, sd))
}
for (fm in list(y ~ 1, y ~ x, y ~ x + (1 | g))) {
  cat("\n==", deparse(fm), "\n")
  t1 <- system.time(a <- tryCatch(smp(fm, family = gen_extreme_value(),
                                      data = dv), error = function(e) e))
  fit <- suppressWarnings(frm(fm, family = gen_extreme_value(), data = dv))
  t2 <- system.time(b <- tryCatch(smp(fit), error = function(e) e))
  if (inherits(a, "error")) {
    cat("formula route ERROR:", conditionMessage(a), "\n"); next
  }
  if (inherits(b, "error")) {
    cat("fit route ERROR:", conditionMessage(b), "\n"); next
  }
  vars <- intersect(c("b_Intercept", "b_x", "sigma", "xi", "sd_g__Intercept"),
                    colnames(as.matrix(a)))
  sa <- summ(a, vars); sb <- summ(b, vars)
  cat("formula route:\n"); print(round(sa, 3))
  cat("fit route:\n"); print(round(sb, 3))
  # difference of means in units of the posterior sd's Monte Carlo error
  z <- (sa["mean", ] - sb["mean", ]) /
    sqrt(sa["sd", ]^2 / 500 + sb["sd", ]^2 / 500)
  cat("z of mean difference (rough, ignores autocorrelation):",
      paste(names(z), round(z, 2)), "\n")
  cat("elapsed formula", round(t1[3], 1), "fit", round(t2[3], 1), "\n")
}
cat("\nDONE\n")
