# Review lane osa: the upper end of a continuous residual under "cdf".
# Reference: -qnorm(S(y)) from pnorm()/pt() at the fit's estimates.
#   Rscript dev/review-osa-edge2.R [lib]
args <- commandArgs(TRUE)
lib1 <- if (length(args) >= 1) args[1] else "C:/Users/adf44/source/r/wt-osa-lib"
.libPaths(c(if (lib1 != "none") lib1, "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-edge2\n")
set.seed(410)
n <- 2000
d <- data.frame(x = rnorm(n))
d$y <- 1 + 0.5 * d$x + rnorm(n)
for (k in c(7, 8, 8.5, 9, 12, -9, -12)) {
  dk <- d
  dk$y[1] <- 1 + 0.5 * dk$x[1] + k
  f1 <- frm(y ~ x, family = gaussian(), data = dk)
  f2 <- frm(y | trunc(lb = -1e3) ~ x, family = gaussian(), data = dk)
  dp <- lapply(frmtmb:::eval_dpars(f2)[[1]], as.numeric)
  z <- (dk$y[1] - dp$mu[1]) / dp$sigma[1]
  w <- character(0)
  r2 <- withCallingHandlers(residuals(f2, type = "osa")[1, "Estimate"],
                            warning = function(c) {
                              w <<- c(w, conditionMessage(c))
                              invokeRestart("muffleWarning")
                            })
  r1 <- residuals(f1, type = "osa")[1, "Estimate"]
  cat(sprintf(paste("outlier %+5.1f: fullGaussian %8.4f | trunc() via cdf",
                    "%8s | (y - mu) / sigma %8.4f%s\n"), k, r1,
              format(r2, digits = 6), z,
              if (length(w)) paste("  WARN:", substr(w[1], 1, 90)) else ""))
}
cat("END review-osa-edge2\n")
