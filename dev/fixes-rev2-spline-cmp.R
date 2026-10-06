# Reviewer of lane fixes, re-check: compare dev/fixes-rev2-spline.R's
# two outputs (lane against base).
a <- readRDS("dev/fixes-rev2-log/spline-wt-fixes-lib.rds")
b <- readRDS("dev/fixes-rev2-log/spline-rellib-r5.rds")
r <- function(u, v) max(abs(u - v)) / max(abs(v))
for (k in names(a)) {
  cat(sprintf(paste0("%-4s dlogLik %.3g; curve est %.3g se %.3g ",
                     "lower_sim %.3g; deriv est %.3g se %.3g; feature %.3g\n"),
              k, a[[k]]$ll - b[[k]]$ll, r(a[[k]]$est, b[[k]]$est),
              r(a[[k]]$se, b[[k]]$se), r(a[[k]]$lsim, b[[k]]$lsim),
              r(a[[k]]$dest, b[[k]]$dest), r(a[[k]]$dse, b[[k]]$dse),
              r(a[[k]]$feat, b[[k]]$feat)))
  cat("   fixef lane:", format(a[[k]]$fe, digits = 5), "\n   fixef base:",
      format(b[[k]]$fe, digits = 5), "\n")
}
