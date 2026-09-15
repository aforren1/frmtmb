# genrev round 2: the truncation error of the naive lognormal mean.
# The lane's construction (dev/generics-scale2.R) verbatim, then the
# statistic the page NAMES ("relative, at the worst row") beside the one
# the script COMPUTES (max |diff| / max fitted), then the closed form.
LIB <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(2026)
n <- 400L
dd <- data.frame(x = rnorm(n), g = factor(rep(1:20, 20)))
eta <- 8 + 0.4 * dd$x + rnorm(20, 0, 0.3)[dd$g]
dd$y <- exp(rnorm(n, eta, 0.4))
say <- function(...) cat(sprintf(...))
for (lb in c(2000, 500, 5000)) {
  d3 <- dd[dd$y > lb, ]
  d3$lbv <- lb
  f3 <- frm(bf(y | trunc(lb = lbv) ~ x + (1 | g)) + lognormal(),
            data = d3)
  mu <- predict(f3, type = "link"); s <- sigma(f3); ft <- fitted(f3)
  naive <- exp(mu + s^2 / 2)
  a <- (log(lb) - mu) / s
  # truncated-below lognormal mean over the untruncated one
  closed <- 1 - stats::pnorm(-a) / stats::pnorm(s - a)
  say("lb = %5d  rows %3d  bound position a: %.2f .. %.2f sigma\n",
      lb, nrow(d3), min(a), max(a))
  say("   lane's statistic max|diff|/max(fitted)   %.4f\n",
      max(abs(ft - naive)) / max(abs(ft)))
  say("   relative error at the worst row          %.4f\n",
      max(abs(ft - naive) / ft))
  say("   median row                               %.4f\n",
      stats::median(abs(ft - naive) / ft))
  say("   closed form 1 - Phi(-a)/Phi(s-a), max    %.4f  max |diff| vs row %.2e\n",
      max(closed), max(abs(closed - abs(ft - naive) / ft)))
}
cat("GENREVDONE\n")
