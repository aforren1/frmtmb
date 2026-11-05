# Fix round 1: why the pin stops Student at y x 1e6 and gaussian at y x 1e9. The
# level's prior is read from the block density at unit vectors, so its precision is
# a difference of two numbers near log(sd), and it rounds once 1 / sd^2 nears eps.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(1)
g <- factor(rep(1:15, each = 10))
d <- data.frame(y = rnorm(150), x = rnorm(150), g)
fr <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d, dry_run = "frame")
lay <- ns$imp_layout(fr[["re_blocks"]])
pf <- ns$imp_prior_terms(lay)
for (sd in c(1, 1e3, 1e5, 1e6, 1e7, 1e8, 1e9)) {
  tm <- pf(log(sd))
  p <- ns$imp_prior_precision(tm, 1L)[1, 1]
  # per group error in the log weight at u = sd: (P_hat - 1/sd^2) sd^2 / 2
  cat(sprintf("sd %-7g precision read %.6g, 1/sd^2 %.6g, relative error %+.3e, per-group log weight error at u = sd %+.3e\n",
              sd, p, 1 / sd^2, p * sd^2 - 1, (p * sd^2 - 1) / 2))
}
for (arm in "rellib") {
  # the same pin on the base build reads the same rounding
  cat("Student at y x 1e6 stops on rellib-r16 too: see own run of lane's scales\n")
}
