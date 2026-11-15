# Lane slopescale: is the gradient warning of slopescale-cost.R's gau1e5
# design (y ~ x + z + (1 + x + z | g), 1000 groups of 100) the base's
# too? Prints the conditions and the optimum.
# Usage: Rscript slopescale-gau1e5.R before|after
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", find.package("frmtmb"), "\n")
set.seed(2)
g <- factor(rep(1:1000, each = 100))
x <- rnorm(1e5)
z <- rnorm(1e5)
y <- 1 + rnorm(1000, 0, .5)[g] + rnorm(1000, 0, .3)[g] * x + 0.2 * z +
  rnorm(1e5)
f <- withCallingHandlers(
  frm(y ~ x + z + (1 + x + z | g), data = data.frame(y, g, x, z)),
  warning = function(w) {
    cat("WARN", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  })
print(f$opt$par, digits = 12)
cat("logLik", format(as.numeric(logLik(f)), digits = 15), "autoscaled",
    isTRUE(f$autoscaled), "\n")
