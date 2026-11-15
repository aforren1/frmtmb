# Short check of parity2's fix round: the width test of
# test-ce-new-levels-snl.R (min(boot 200 / Wald) > 0.8) over 20 ulp
# perturbations of y and over 30 boot seeds at the test's data.
#   Rscript dev/reviews/parity2-short-width.R > <out> 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-parity2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
cat("BLAS:", extSoftVersion()[["BLAS"]], "\n")
set.seed(100)
d <- data.frame(x = rnorm(60), g = factor(rep(1:30, each = 2)))
d$y <- 1 + 0.5 * d$x + rnorm(30)[d$g] + rnorm(60, sd = 2)
ic <- list(x = c(-1, 0, 1))
wid <- function(ce) ce$x$upper__ - ce$x$lower__
ratio <- function(dd, seed = 2) {
  w <- character()
  fit <- withCallingHandlers(
    frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  wb <- wid(suppressMessages(conditional_effects(
    fit, "x", re_formula = NULL, band = "boot", boot = 200, seed = seed,
    int_conditions = ic)))
  ww <- wid(conditional_effects(fit, "x", re_formula = NULL,
                                int_conditions = ic))
  c(min = min(wb / ww), nwarn = length(w))
}
r0 <- ratio(d)
cat("test data, seed 2: min ratio", format(r0[1], digits = 4),
    "warnings", r0[2], "\n")
t0 <- proc.time()[["elapsed"]]
ulp <- t(vapply(1:20, function(k) {
  dk <- d
  dk$y <- d$y * (1 + 2e-15 * k)
  ratio(dk)
}, c(min = 0, nwarn = 0)))
cat("ulp k=1..20 min ratio:", format(ulp[, 1], digits = 4), "\n")
cat("ulp range", format(range(ulp[, 1]), digits = 4), "warnings",
    sum(ulp[, 2]), "\n")
seeds <- t(vapply(1:30, function(s) ratio(d, s), c(min = 0, nwarn = 0)))
cat("boot seeds 1..30 min ratio:", format(seeds[, 1], digits = 4), "\n")
cat("seeds range", format(range(seeds[, 1]), digits = 4), "below 0.8:",
    sum(seeds[, 1] <= 0.8), "of 30; mean", format(mean(seeds[, 1]),
                                                 digits = 4),
    "sd", format(sd(seeds[, 1]), digits = 3), "\n")
cat("elapsed", round(proc.time()[["elapsed"]] - t0, 1), "s\ndone\n")
