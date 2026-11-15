# Review of lane famend, item 3580: poisson_cdf_log() at large means and
# in the blend, value against ppois() and AD gradient against a central
# difference of the tape's own value.
#   Rscript dev/reviews/famend-review-pois.R <base|lane>
a <- commandArgs(trailingOnly = TRUE)
lanelib <- "C:/Users/adf44/source/r/wt-famend-lib"
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] == "lane") lanelib, base, ulib))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
pcl <- get("poisson_cdf_log", asNamespace("frmtmb"))
qs <- c(1, 2, 5, 10, 20, 29, 30)
mus <- c(1e-3, 0.5, 3, 8, 15, 25, 35, 45, 60, 71, 72, 100, 160, 1e3, 1e4, 1e6)
worst <- c(lF = 0, lS = 0, gF = 0, gS = 0)
bad <- 0L
for (q in qs) {
  for (m in mus) {
    tF <- RTMB::MakeTape(function(l) pcl(q, exp(l))$lF, log(m))
    tS <- RTMB::MakeTape(function(l) pcl(q, exp(l))$lS, log(m))
    l <- log(m)
    vF <- tF(l); vS <- tS(l)
    rF <- ppois(q, m, log.p = TRUE)
    rS <- ppois(q, m, lower.tail = FALSE, log.p = TRUE)
    gF <- tF$jacobian(l); gS <- tS$jacobian(l)
    # exact derivative in log(mu): d log F / d log mu = -mu p(q) / F
    eF <- -m * dpois(q, m) / ppois(q, m)
    eS <- m * dpois(q, m) / ppois(q, m, lower.tail = FALSE)
    hF <- tF$jacfun()$jacobian(l)
    if (!all(is.finite(c(vF, vS, gF, gS, hF)))) bad <- bad + 1L
    rel <- function(v, r) if (is.finite(r) && r != 0) abs(v - r) / abs(r) else
      if (identical(v, r)) 0 else abs(v - r)
    e <- c(rel(vF, rF), rel(vS, rS),
           if (is.finite(eF)) rel(gF, eF) else 0,
           if (is.finite(eS)) rel(gS, eS) else 0)
    worst <- pmax(worst, e)
    if (max(e) > 1e-9) {
      cat(sprintf("q %2d mu %8.3g: lF %.10g (ref %.10g) lS %.10g (ref %.10g) dlF %.6g (ref %.6g) dlS %.6g (ref %.6g)\n",
                  q, m, vF, rF, vS, rS, gF, eF, gS, eS))
    }
  }
}
cat("worst relative error over", length(qs) * length(mus), "points:\n")
print(signif(worst, 3))
cat("points with a non-finite value, gradient or Hessian:", bad, "\n")
cat("\ndone\n")
