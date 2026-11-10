# Release review of 0.81.0: the fitted objective's GEV barrier just
# outside the support, against the exact reader, at t = 1 + xi z from
# +1e-6 down to -1e-2. The record says it falls by 1e8 per unit of t.
#   Rscript dev/review-rel081-barrier.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), format(packageVersion("frmtmb")), "\n")
rl <- get("row_lpdf", asNamespace("frmtmb"))
fam <- gen_extreme_value()
tt <- c(1e-6, 1e-9, 1e-12, 0, -1e-12, -1e-9, -1e-6, -1e-4, -1e-2)
for (xi in c(-0.9, -0.3, 0.3, 10)) {
  mu <- 0; sg <- 1
  # t = 1 + xi (y - mu) / sigma, so y = mu + sigma (t - 1) / xi
  y <- mu + sg * (tt - 1) / xi
  p <- list(mu = rep(mu, length(y)), sigma = rep(sg, length(y)),
            xi = rep(xi, length(y)))
  a <- as.numeric(rl(fam, y, y, p, list(), NULL))
  b <- tryCatch(as.numeric(rl(fam, y, y, p, list(), NULL, exact = "TRUE" == "TRUE")),
                error = function(e) rep(NA, length(y)))
  cat(sprintf("xi %5.1f  t %9.0e  objective %12.5g  exact %12.5g\n",
              xi, tt, a, b), sep = "")
}
cat("BARRIER DONE\n")
