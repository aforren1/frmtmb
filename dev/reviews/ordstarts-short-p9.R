# Reviewer short check of lane ordstarts (build fix round 1), item 2:
# constant weights 1/1000, 1 and 10 on cum_sratio without random effects,
# n 1000, seeds 1 to 20. Same maximum? Per-unit logLik (logLik / w) and
# the largest parameter difference from w = 1, for the default (D) and
# for starts = 0 (A0, the optimizer alone).
#   Rscript ordstarts-short-p9.R > ordstarts-short-p9.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
mk <- function(seed, n) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
q <- function(e) suppressWarnings(suppressMessages(e))
fam <- mixture(cumulative(), sratio())
nd <- 0
na0 <- 0
for (seed in 1:20) {
  d <- mk(seed, 1000)
  res <- list()
  for (w in c(1, 1 / 1000, 10)) {
    d$w <- w
    for (arm in c("A0", "D")) {
      ctl <- if (arm == "A0") frmtmb_control(starts = 0) else
        frmtmb_control()
      f <- q(frm(bf(y | weights(w) ~ x + z), family = fam, data = d,
                 control = ctl))
      ms <- f$opt$mixture_starts
      res[[paste(arm, w)]] <- list(ll = as.numeric(logLik(f)) / w,
                                   par = f$opt$par, code = f$opt$convergence,
                                   sw = if (is.null(ms)) NA else
                                     ms[["swapped"]])
    }
  }
  for (arm in c("A0", "D")) {
    r1 <- res[[paste(arm, 1)]]
    for (w in c(1 / 1000, 10)) {
      r <- res[[paste(arm, w)]]
      dll <- r$ll - r1$ll
      dp <- max(abs(r$par - r1$par))
      bad <- abs(dll) > 1e-4 || dp > 1e-2
      if (bad && arm == "D") nd <- nd + 1
      if (bad && arm == "A0") na0 <- na0 + 1
      cat(sprintf("seed %2d %-2s w %-6g | per-unit ll - w=1's %+.3g | max|dpar| %.3g | code %d (w=1: %d) | swapped %s (w=1: %s)%s\n",
                  seed, arm, w, dll, dp, r$code, r1$code, r$sw, r1$sw,
                  if (bad) "  DIFFERENT" else ""))
    }
  }
}
cat("D differs from w = 1 on", nd, "of 40 (seed, w) pairs; A0 on", na0,
    "of 40\n")
cat("DONE p9\n")
