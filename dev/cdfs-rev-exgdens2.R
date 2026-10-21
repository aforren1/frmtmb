# Reviewer: exgaussian()'s density (RTMBdist::dexgauss, unchanged by
# this lane) against brms::dexgaussian() as sigma / beta grows; on the
# base build. ARM=base or lane.
arm <- Sys.getenv("ARM", "base")
lp <- c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") lp <- c("C:/Users/adf44/source/r/wt-cdfs-rev-lib", lp)
.libPaths(lp)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), " RTMBdist", format(packageVersion("RTMBdist")), "\n")
f <- frmtmb:::as_frmtmb_family(exgaussian())
y <- seq(-3, 3, by = 0.5)
for (u in c(1, 2, 5, 8, 10, 12, 15, 20, 40, 100, 1e3)) {
  bt <- 1 / u
  a <- f$lpdf(y, list(mu = 0, sigma = 1, beta = bt), list())
  b <- brms::dexgaussian(y, 0, 1, bt, log = TRUE)
  cat(sprintf("sigma / beta %-6g max |lpdf - brms| %.3e  at y = %g\n", u, max(abs(a - b)),
              y[which.max(abs(a - b))]))
}
cat("RTMBdist::dexgauss directly, mu 0 sigma 1 lambda 20, y = 0:",
    RTMBdist::dexgauss(0, 0, 1, 20, log = TRUE), " closed form:",
    log(20) + 20 * (0 - 0) + 20^2 / 2 + pnorm(0 - 20, log.p = TRUE), "\n")
