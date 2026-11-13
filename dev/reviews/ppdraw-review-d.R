# Review D: a zero inflation at the end of its range on LARGE data. The
# flat side's change is about the optimizer's last gradient along the
# parameter, not the sample size; does the optimizer leave it below 1e-4
# when the objective is large? Also a hurdle on large data.
# Usage: Rscript ppdraw-review-d.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (length(args) && args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
source("C:/Users/adf44/source/r/frmtmb-wt-ppdraw/dev/reviews/ppdraw-review-common.R")
one <- function(n, seed, fam = "zip") {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  if (fam == "zip") {
    # mean 0.7: many zeros, all from the poisson
    d$y <- stats::rpois(n, exp(-0.35 + 0.3 * d$x))
    t0 <- proc.time()[["elapsed"]]
    fit <- q(frm(bf(y ~ x), family = zero_inflated_poisson(), data = d))
  } else {
    d$y <- stats::rnbinom(n, mu = exp(0.5 + 0.3 * d$x), size = 2)
    t0 <- proc.time()[["elapsed"]]
    fit <- q(frm(bf(y ~ x), family = zero_inflated_negbinomial(), data = d))
  }
  tf <- proc.time()[["elapsed"]] - t0
  fs <- flat_sides(fit)
  j <- grep("^zi_", rownames(fs))
  g <- fit$obj$gr(fit$opt$par)
  r <- pz(fit, ndraws = 200L, seed = 1)
  cat(sprintf(paste("%-5s n %7d seed %2d fit %.1fs objective %.6g conv %s",
                    "zi est %.3f sd %.3g minus2sd %.3g plus2sd %.3g",
                    "grad_zi %.3g z %.2f held %s share_na %.3f\n"),
              fam, n, seed, tf, fit$opt$objective,
              fit$opt$convergence, fs[j, "est"], fs[j, "sd"],
              fs[j, "minus2sd"], fs[j, "plus2sd"], g[j], r$z,
              length(r$held) > 0, r$na))
  if (length(r$warns)) cat("  warnings:", substr(r$warns, 1, 150), sep = "\n   ")
}
for (n in c(2e4, 2e5, 1e6)) for (seed in 1:2) one(as.integer(n), seed)
for (n in c(2e5, 1e6)) one(as.integer(n), 1, "zinb")
