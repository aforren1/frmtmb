# Review of lane parity2, item 1497, part B again on a NON-singular fit
# with strong shrinkage: 30 groups of 2 rows, sd_g 1, sigma 2. The data
# seed is the first of 100.. whose ML sd_g is at least 0.7; the data go
# to snl2-data.rds for dev/reviews/parity2-review-snl2-brms.R.
#   Rscript dev/reviews/parity2-review-snl2.R
.libPaths(c("C:/Users/adf44/source/r/wt-parity2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
G <- 30
for (s in 100:200) {
  set.seed(s)
  d <- data.frame(x = rnorm(2 * G), g = factor(rep(1:G, each = 2)))
  d$y <- 1 + 0.5 * d$x + rnorm(G)[d$g] + rnorm(2 * G, sd = 2)
  fb <- suppressMessages(suppressWarnings(
    frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)))
  sdg <- VarCorr(fb)$g$sd[1, "Estimate"]
  if (sdg >= 0.7) break
}
cat("data seed", s, " ML sd_g", sdg, " sigma",
    VarCorr(fb)$residual__$sd[1, "Estimate"], "\n")
saveRDS(d, "C:/Users/adf44/source/r/frmtmb-wt-parity2/dev/reviews/parity2-review-out/snl2-data.rds")
bk <- fb$frame$re_blocks[[1]]
cat("sd of the 30 modes:", sd(fb$estimates$b[bk$b_idx]), "\n")
ic <- list(x = c(-1, 0, 1))
wid <- function(ce) ce$x$upper__ - ce$x$lower__
for (snl in c("uncertainty", "gaussian", "old_levels")) {
  ce <- suppressMessages(conditional_effects(
    fb, "x", int_conditions = ic, re_formula = NULL, band = "boot",
    boot = 400, seed = 2, sample_new_levels = snl))
  cat("fit boot 400 seed 2", snl, "width:", format(wid(ce), digits = 4), "\n")
}
cw <- conditional_effects(fb, "x", int_conditions = ic, re_formula = NULL)
cat("fit wald width:", format(wid(cw), digits = 4), "\n")
cwo <- conditional_effects(fb, "x", int_conditions = ic, re_formula = NULL,
                           sample_new_levels = "old_levels")
cat("fit wald with old_levels identical to default:", identical(cw, cwo), "\n")
smp <- suppressMessages(suppressWarnings(
  frm_sample(fb, chains = 2, iter = 2000, seed = 1, refresh = 0)))
for (snl in c("uncertainty", "gaussian")) {
  ce <- conditional_effects(smp, "x", int_conditions = ic,
                            re_formula = NULL, seed = 2,
                            sample_new_levels = snl)
  cat("NUTS draws", snl, "width:", format(wid(ce), digits = 4), "\n")
}
r <- smp$draws[, grep("^r_g\\[", colnames(smp$draws))]
cat("sd over draws and levels of r_g draws:", sd(as.vector(r)), "\n")
