# Object size of a draws object with and without saved warmup, on the
# reviewer's 400-group case (dev/review-priordraws-probe.R, P7): data
# seed 1, sampler seed 2, 2 chains of 1000. Also times inc_warmup.
#   Rscript dev/priordraws-size.R <label>
lab <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(1)
big <- data.frame(x = rnorm(2000), g = factor(rep(1:400, 5)))
big$y <- 1 + 0.5 * big$x + rnorm(400)[big$g] + rnorm(2000)
pr <- set_prior("normal(0, 1)", class = "b")
mb <- function(o) as.numeric(object.size(o)) / 2^20
for (sw in c(TRUE, FALSE)) {
  bd <- q(frm_sample(bf(y ~ x + (1 | g)), data = big, family = gaussian(),
                     prior = pr, chains = 2, iter = 1000, seed = 2,
                     refresh = 0, save_warmup = sw))
  f <- tempfile()
  saveRDS(bd, f)
  t <- system.time(a <- as_draws_array(bd, inc_warmup = TRUE))[["elapsed"]]
  cat(sprintf("%s save_warmup=%s: draws %.1f MB, warmup %.1f MB, stanfit %.1f MB, total %.1f MB, saveRDS %.1f MB; inc_warmup array %d x %d x %d in %.2f s\n",
              lab, sw, mb(bd$draws), mb(bd$warmup), mb(bd$stanfit), mb(bd),
              file.size(f) / 2^20, dim(a)[1], dim(a)[2], dim(a)[3], t))
  if (sw) saveRDS(unclass(a), sprintf(
    "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/priordraws-log/size-%s.rds",
    lab))
}
