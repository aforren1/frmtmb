# Short check, item 2: laplace = TRUE's wall clock, base (centered, 1e-6) against
# now (factor, 1e-12), at ordinary noise, by n. Seed 6701; 2 chains of 300.
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1], ";")[[1]], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
for (n in c(40L, 100L, 200L)) {
  set.seed(6701)
  d <- data.frame(x = round(runif(n, 0, 6), 3))
  d$y <- sin(2 * d$x) + 0.3 * rnorm(n)
  t0 <- proc.time()[["elapsed"]]
  fl <- suppressMessages(suppressWarnings(frm_sample(bf(y ~ gp(x)), data = d, family = gaussian(), laplace = TRUE,
                                    chains = 2, iter = 300, warmup = 150, seed = 8, refresh = 0, cores = 1)))
  dd <- posterior::as_draws_df(fl)
  cat(sprintf("LCOST %s n %d: %.0f s, sigma median %.4f sdgp median %.3f\n", basename(a[1]), n,
              proc.time()[["elapsed"]] - t0, median(dd$sigma), median(dd$sdgp_gpx)))
}
