# Review of lane squant: relative error of the tape's d lF / d mu of
# quantile_elf()'s cdf_log against the analytic -f/F, over t, for the
# settings of squant-review-tape4.txt (p 0.1) and the default.
# Run (worktree root):
#   Rscript dev/reviews/squant-review-tape5.R > dev/reviews/squant-review-tape5.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
dp <- function(mu, sigma, p) list(mu = mu, sigma = sigma, quantile = p)
for (p in c(0.1, 0.5, 0.9)) {
  for (lam in c(0.25, 1)) {
    fam <- quantile_elf(smoothing = lam)
    out <- character()
    for (t in c(-10, -20, -40, -60, -80, -100, -125, -150, -200, -300, -500,
                -689)) {
      q <- t * lam
      tp <- RTMB::MakeTape(function(m) {
        fam$cdf_log(q, dp(m, 1, p), list())[["lF"]]
      }, 0)
      g <- as.numeric(tp$jacobian(0))
      v <- fam$cdf_log(q, dp(0, 1, p), list())[["lF"]]
      an <- -exp(fam$lpdf(q, dp(0, 1, p), list()) - v)
      out <- c(out, sprintf("%g:%.1f:%.2g", t, v, abs(g / an - 1)))
    }
    cat(sprintf("p %.1f lam %.2f (t:lF:rel err) %s\n", p, lam,
                paste(out, collapse = " ")))
  }
}
cat("DONE\n")
