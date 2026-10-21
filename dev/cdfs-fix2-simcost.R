# Fix round 2, N3: what inversion costs against rejection in
# sim_response(), so the fallback can be placed: rejection with the
# inverse-CDF fallback after 20 passes (the default), against inversion
# for every row rejection's first pass leaves outside (inv_after = 0).
#   Rscript dev/cdfs-fix2-simcost.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), "\n")
tm <- function(lab, fam, dp, av, n, inv_after) {
  set.seed(1)
  t0 <- proc.time()[["elapsed"]]
  for (r in 1:20) y <- ns$sim_response(ns$as_frmtmb_family(fam), dp, av, n,
                                       inv_after = inv_after)
  cat(sprintf("%-46s inv_after %2d: %7.1f ms per draw of %d rows\n", lab,
              inv_after, 1000 * (proc.time()[["elapsed"]] - t0) / 20, n))
}
n <- 10000
# most of the mass in the window: rejection's ordinary case
for (ia in c(20L, 0L)) {
  tm("gaussian, trunc(lb = -1), keeps 84 percent", gaussian(),
     list(mu = rep(0, n), sigma = rep(1, n)), list(trunc_lb = rep(-1, n)),
     n, ia)
  tm("poisson mu 3, trunc(lb = 1), keeps 95 percent", poisson(),
     list(mu = rep(3, n)), list(trunc_lb = rep(1, n)), n, ia)
  tm("weibull, trunc(lb = 0.5), keeps 10 percent", weibull(),
     list(mu = rep(0.2, n), shape = rep(1.5, n)),
     list(trunc_lb = rep(0.5, n)), n, ia)
  tm("Gamma shape 3, trunc(lb = 6), keeps 6e-3", Gamma(link = "log"),
     list(mu = rep(1, n), shape = rep(3, n)), list(trunc_lb = rep(6, n)),
     n, ia)
}
cat("END\n")
