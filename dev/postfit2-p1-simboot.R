# Punch round 1, the frm_bootstrap() decision of 2026-09-29: the
# dev/simnewdata-boot.R setup (y ~ s(x), 200 rows, seed 21, 40 refits,
# seed 1), the bootstrap SD of the intercept and of the smooth's SD
# parameter at re_formula = NA and NULL, on the base build and this lane.
#   Rscript dev/postfit2-p1-simboot.R <base|lane>
arm <- commandArgs(trailingOnly = TRUE)[1]
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "base")) libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(21)
d <- data.frame(x = stats::runif(200), g = factor(rep(1:10, 20)))
d$y <- 2 * sin(2 * pi * d$x) + stats::rnorm(10, 0, 0.5)[d$g] +
  stats::rnorm(200, 0, 0.3)
fs <- frm(bf(y ~ s(x)), data = d)
sm_sd <- function(f) exp(f$estimates$theta[[1]])
FUN <- function(f) c(fixef(f, flatten = TRUE)[1], sd_s = sm_sd(f))
for (rf in list(NA, NULL)) {
  b <- frm_bootstrap(fs, FUN = FUN, nsim = 40, seed = 1, re_formula = rf)
  cat(sprintf("%s: re_formula %s, estimate %.4f / %.4f; bootstrap sd %.4f / %.4f; converged %d of 40\n",
              arm, if (is.null(rf)) "NULL" else "NA", b$t0[1], b$t0[2],
              stats::sd(b$t[, 1], na.rm = TRUE),
              stats::sd(b$t[, 2], na.rm = TRUE), sum(b$converged)))
}
cat(sprintf("%s: intercept Wald se %.4f\n", arm, sqrt(vcov(fs)[1, 1])))
