# Lane slopetrig: the year fit of the MEDIUM backlog item (seed 2, both
# forms), lane default: correlation parameter, correlation, its SE and
# the warnings. Usage: Rscript slopetrig-probe2.R
.libPaths(c("C:/Users/adf44/source/r/wt-slopetrig-lib", "C:/Users/adf44/source/r/rellib-r22", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(2); g <- factor(rep(1:20, each = 15)); z <- rnorm(300)
eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * z + 0.3 * z
d <- data.frame(g = g, x = 2000 + 3 * z, y = eta + rnorm(300))
for (fo in list(y ~ 1 + (1 + x | g), y ~ x + (1 + x | g))) {
  w <- character(0)
  f <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  vc <- suppressWarnings(VarCorr(f))$g$cor[2, , 1]
  cat(deparse(fo), "t", f$opt$par[names(f$opt$par) == "theta"][3],
      "cor", format(vc[["Estimate"]], digits = 10), "se", vc[["Est.Error"]],
      "\n  warnings:", substr(w, 1, 100), sep = " ", "\n")
}
