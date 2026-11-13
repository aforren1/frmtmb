# Short check: the draw_edge_c cliff on a continuous family of data sets.
# The lane's "slope sd 0.05 seed 5" data (changed, r 0.495), with the
# group-slope component scaled by lambda, so the slope log sd and its
# upper side's fall move smoothly. Per lambda: the decision statistic
# (fall of the log-likelihood at 1 conditional sd up; steep when >= 2),
# whether the law changed the draws, and a new level's predictive spread
# at x = 2 (4000 draws).
# Usage: Rscript ppdraw-short-cliff.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
lane <- exists("draw_edge_csd", envir = asNamespace("frmtmb"))
set.seed(5)
dd <- data.frame(g = factor(rep(1:20, each = 10)), x = stats::rnorm(200))
a <- stats::rnorm(20, 0, 0.7)[dd$g]
b <- stats::rnorm(20, 0, 0.05)[dd$g] * dd$x
e <- stats::rnorm(200)
nd <- data.frame(x = 2, g = factor("new"))
cat(sprintf("%6s %9s %8s %8s %9s %8s %8s %8s %8s\n", "lambda", "theta2",
            "sd", "csd", "fall@csd", "changed", "pred_sd", "w90", "w98"))
for (lam in seq(0.5, 2, by = 0.1)) {
  dd$y <- 1 + 0.5 * dd$x + a + lam * b + e
  fit <- suppressWarnings(suppressMessages(frm(bf(y ~ x + (1 + x | g)),
                                               data = dd)))
  ds <- frmtmb:::fit_draw_space(fit)
  V <- ds$V
  j <- which(rownames(V) == "theta_2")
  p <- fit$opt$par
  off <- nrow(V) - length(p)
  sdv <- sqrt(diag(V))
  fall <- csdj <- NA
  if (lane) {
    csd <- frmtmb:::draw_edge_csd(fit, V, sdv, off, FALSE)
    csdj <- csd[j]
    f0 <- fit$obj$fn(p)
    qq <- p
    qq[j - off] <- p[j - off] + min(csd[j], sdv[j])
    fall <- fit$obj$fn(qq) - f0
    fit$obj$fn(p)
  }
  said <- FALSE
  set.seed(1)
  pr <- withCallingHandlers(
    suppressWarnings(predict(fit, newdata = nd, ndraws = 4000,
                             summary = FALSE, allow_new_levels = TRUE)),
    frmtmb_draw_held = function(m) {
      said <<- TRUE
      invokeRestart("muffleMessage")
    })
  qs <- stats::quantile(pr, c(0.01, 0.05, 0.95, 0.99))
  cat(sprintf("%6.2f %9.3f %8.3g %8.3g %9.4f %8s %8.4f %8.4f %8.4f\n", lam,
              p[j - off], sdv[j], csdj, fall, said, stats::sd(pr),
              qs[3] - qs[2], qs[4] - qs[1]))
}
