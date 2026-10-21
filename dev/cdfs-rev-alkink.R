# Reviewer: asym_laplace AD vs FD on the truncated objective is the
# density's kink, not the CDF: the same check untruncated, and central
# differences at 1e-7 / 1e-9.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(77)
n <- 200
x <- rnorm(n)
y <- rnorm(n)
keep <- y >= -1.5 & y <= 2
d <- data.frame(x = x[keep], y = y[keep], lbv = -1.5, ubv = 2)
for (fo in list(bf(y ~ x, sigma ~ 1, quantile ~ 1),
                bf(y | trunc(lb = lbv, ub = ubv) ~ x, sigma ~ 1, quantile ~ 1))) {
  ob <- frm(fo, family = asym_laplace(), data = d, dry_run = "objective")$obj
  p <- ob$par
  g <- as.numeric(ob$gr(p))
  fdr <- numDeriv::grad(ob$fn, p, method = "Richardson")
  fd7 <- vapply(seq_along(p), function(j) { e <- replace(numeric(length(p)), j, 1e-7)
    (ob$fn(p + e) - ob$fn(p - e)) / 2e-7 }, 0)
  cat(deparse(fo$formula), "\n  AD", signif(g, 7), "\n  FD Richardson", signif(fdr, 7),
      "\n  FD 1e-7", signif(fd7, 7), "\n")
}
