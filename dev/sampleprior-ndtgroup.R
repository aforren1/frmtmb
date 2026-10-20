# Lane sampleprior, fix round 2, finding 9: frm_sample()'s default on a
# frmtmb.eam ndt FRACTION under ndt_group(), against the ML fraction,
# for wiener(), lba() and rdm() on three data seeds each. Sampler seed 4,
# 2 chains x 400.
#   Rscript dev/sampleprior-ndtgroup.R <lane|r9> <wiener|lba|rdm>
a <- commandArgs(TRUE)
arm <- a[1]
which <- a[2]
lib <- if (identical(arm, "r9")) character(0) else
  "C:/Users/adf44/source/r/wt-sampleprior-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
  library(frmtmb.eam)
})
cat("arm", arm, "lib:", find.package("frmtmb.sample"), "\n")
for (seed in 20261121:20261123) {
  set.seed(seed)
  if (which == "wiener") {
    d <- ddm_simulate(120, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
    form <- bf(rt | dec(upper) + ndt_group(subj) ~ 1)
    fam <- wiener()
  } else if (which == "lba") {
    d <- lba_simulate(200, v = c(2.4, 1.4), A = 0.5, k = 0.4, ndt = 0.2)
    form <- bf(rt | vint(choice) + ndt_group(subj) ~ 1)
    fam <- lba(2)
  } else {
    d <- rdm_simulate(200, v = c(2.4, 1.4), A = 0.5, k = 0.5, ndt = 0.2)
    form <- bf(rt | vint(choice) + ndt_group(subj) ~ 1)
    fam <- rdm(2)
  }
  d$subj <- factor(rep(1:4, length.out = nrow(d)))
  fit <- suppressWarnings(frm(form, family = fam, data = d))
  s <- Filter(function(z) identical(z$dpar, "ndt"),
              unclass(frmtmb.sample:::default_priors_for(fit)))[[1L]]
  ml <- frm_linpred(fit, dpar = "ndt", type = "response")[1L]
  t0 <- proc.time()[["elapsed"]]
  ds <- suppressWarnings(suppressMessages(frm_sample(
    fit, chains = 2, iter = 400, seed = 4, refresh = 0)))
  v <- as.matrix(ds)[, "ndt"]
  cat(sprintf(paste0("%s seed %d: default %s ub %s | ML fraction %.4f | ",
                     "posterior q01/50/99 %s, max %.4f | rhat %.3f | ",
                     "%.0f s\n"), which, seed, s$prior, format(s$ub), ml,
              paste(sprintf("%.4f", stats::quantile(v, c(.01, .5, .99))),
                    collapse = " "), max(v),
              posterior::rhat(matrix(v, ncol = 2)),
              proc.time()[["elapsed"]] - t0))
}
