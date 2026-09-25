# Punch round 1, item 7 as the user decided it (2026-09-24):
# frm_bootstrap() keeps its 0.62.0 whole-model bootstrap, which redraws
# the group effects AND the penalized smooth coefficients in every
# replicate. Saves frm_bootstrap() at a fixed seed on a smooth fit and a
# mixed fit, for a bitwise comparison of two builds, and measures the
# default against re_formula = NULL (condition on the fitted random
# effects and smooths) on the smooth fit.
#   SIMNEWDATA_LIB=base Rscript dev/simnewdata-boot.R base
#   Rscript dev/simnewdata-boot.R lane
#   Rscript dev/simnewdata-boot.R compare
arm <- commandArgs(trailingOnly = TRUE)[1]
OUT <- "dev/simnewdata-log"
if (identical(arm, "compare")) {
  a <- readRDS(file.path(OUT, "boot-base.rds"))
  b <- readRDS(file.path(OUT, "boot-lane.rds"))
  for (k in names(a)) {
    cat(sprintf("%-28s %s\n", k, if (identical(a[[k]], b[[k]]))
      "identical" else "DIFFERS"))
  }
  quit(save = "no")
}
source("dev/simnewdata-prelude.R")
suppressMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(21)
d <- data.frame(x = stats::runif(200), g = factor(rep(1:10, 20)))
d$y <- 2 * sin(2 * pi * d$x) + stats::rnorm(10, 0, 0.5)[d$g] +
  stats::rnorm(200, 0, 0.3)
fs <- frm(bf(y ~ s(x)), data = d)
fm <- frm(bf(y ~ x + (1 | g)), data = d)
# the smooth's standard deviation is its block's exp(theta); one fixed
# effect is the intercept
sm_sd <- function(f) exp(f$estimates$theta[[1]])
FUN <- function(f) c(fixef(f, flatten = TRUE)[1], sd_s = sm_sd(f))
res <- list()
res[["smooth default"]] <- frm_bootstrap(fs, FUN = FUN, nsim = 40, seed = 1)
res[["mixed default"]] <- frm_bootstrap(fm, nsim = 40, seed = 1)
res[["smooth NULL"]] <- frm_bootstrap(fs, FUN = FUN, nsim = 40, seed = 1,
                                      re_formula = NULL)
res[["mixed NULL"]] <- frm_bootstrap(fm, nsim = 40, seed = 1,
                                     re_formula = NULL)
for (k in c("smooth default", "smooth NULL")) {
  t <- res[[k]]$t
  cat(sprintf("%s: %s, estimate %.4f / %.4f; bootstrap sd %.4f / %.4f; converged %d of 40\n",
              arm, k, res[[k]]$t0[1], res[[k]]$t0[2],
              stats::sd(t[, 1], na.rm = TRUE), stats::sd(t[, 2], na.rm = TRUE),
              sum(res[[k]]$converged)))
}
cat(sprintf("%s: intercept Wald se %.4f\n", arm, sqrt(vcov(fs)[1, 1])))
saveRDS(res, file.path(OUT, paste0("boot-", arm, ".rds")))
