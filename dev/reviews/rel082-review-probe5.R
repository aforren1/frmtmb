# Release review of 0.82.0, probe 5: the escaped "NA/NaN function
# evaluation" warnings beside code 1 on a bound where the objective is
# Inf (NEWS, prior tables). usage: probe5.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")), "\n")
set.seed(11)
d <- data.frame(x = rnorm(120))
d$yn <- 1 - 2 * d$x + rnorm(120)
ws <- list()
f <- withCallingHandlers(
  frm(yn ~ x, data = d,
      prior = set_prior("lognormal(0, 1)", class = "b", coef = "x", lb = 0)),
  warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
cat("code", f$opt$convergence, "objective", f$opt$objective, "\n")
for (w in ws) {
  cat("W[", class(w)[1], "] call:", deparse1(conditionCall(w))[1], "|",
      substr(conditionMessage(w), 1, 80), "\n")
}
