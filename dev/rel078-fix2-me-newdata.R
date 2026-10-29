# The release review's finding 8, the me() newdata item, measured: on
# new data posterior_epred() plugs the observed x in for the latent value
# and ignores sdx, where brms 2.23.0 draws a new latent per draw from
# normal(x, sdx) (dev/brms-suite/brms/R/prepare_predictions.R, the
# "sample new values of latent variables" branch). The model and data
# are test-rescor-me.R's (data seed 29), sampled with seed 5, 2 chains of
# 800. brms's analogue is computed on frmtmb's own draws: b0 + b (x +
# sdx z), z ~ N(0, 1) per draw (seed 1).
#   Rscript dev/rel078-fix2-me-newdata.R > dev/rel078-log/fix2-me-newdata.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")),
    format(packageVersion("frmtmb.sample")), "\n")
ex <- parse("C:/Users/adf44/source/r/frmtmb-wt-rel078/dev/rel078-fix-rescorme-test.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% c("rm_data", "rm_form", "rm_prior")) eval(e)
}
d <- rm_data()
ds <- suppressWarnings(suppressMessages(frm_sample(
  rm_form(), data = d, family = gaussian(), prior = rm_prior(),
  chains = 2, iter = 800, seed = 5, refresh = 0)))
m <- as.matrix(ds)
nd <- d[1:5, ]
ep_in <- posterior_epred(ds)[, 1, "y"]
ep_new <- posterior_epred(ds, newdata = nd)[, 1, "y"]
nd0 <- nd; nd0$sx <- 0
ep_new0 <- posterior_epred(ds, newdata = nd0)[, 1, "y"]
b0 <- m[, "b_y_Intercept"]; b1 <- m[, "bsp_y_mexosx"]
set.seed(1)
brms_like <- b0 + b1 * (nd$xo[1] + nd$sx[1] * rnorm(nrow(m)))
plug <- b0 + b1 * nd$xo[1]
cat(sprintf("row 1, response y, sd over %d draws of the expected value:\n",
            nrow(m)))
cat(sprintf("  in-sample (latent draws):           %.4f\n", sd(ep_in)))
cat(sprintf("  newdata, frmtmb:                    %.4f\n", sd(ep_new)))
cat(sprintf("  newdata with sx = 0, frmtmb:        %.4f\n", sd(ep_new0)))
cat(sprintf("  plug-in b0 + b x on the draws:      %.4f (max |diff| to frmtmb %.2g)\n",
            sd(plug), max(abs(plug - ep_new))))
cat(sprintf("  brms's draw, b0 + b (x + sdx z):    %.4f\n", sd(brms_like)))
cat(sprintf("  ratio frmtmb / brms analogue:       %.3f\n",
            sd(ep_new) / sd(brms_like)))
q <- function(v) unname(quantile(v, c(0.025, 0.975)))
cat(sprintf("  95%% interval width: frmtmb %.4f, brms analogue %.4f\n",
            diff(q(ep_new)), diff(q(brms_like))))
