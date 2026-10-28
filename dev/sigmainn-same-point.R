# The same model, read at the same point: the lane build's objective at
# rellib-r14's ML estimates carried to the innovation sd (log sigma
# intercept minus 0.5 log gamma(0), slopes unchanged) against
# rellib-r14's -logLik. An identity up to rounding.
# Needs dev/sigmainn-out/route-before.rds (dev/sigmainn-route.R before).
# Usage: Rscript dev/sigmainn-same-point.R
# Log: dev/sigmainn-out/same-point.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
.libPaths(c("C:/Users/adf44/source/r/wt-sigmainn-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
src <- readLines(file.path(wt, "dev", "sigmainn-route.R"))
# the data and formulas of the route script, without its arms
i1 <- grep("^if \\(arm %in% c\\(\"before\", \"after\"\\)\\)", src) - 1
arm <- "none"
eval(parse(text = src[grep("^src <- readLines", src):i1]))
B <- readRDS(file.path(wt, "dev", "sigmainn-out", "route-before.rds"))
for (m in names(dat)) {
  d <- dat[[m]]
  fit <- if (m == "sigx") {
    frm(bf(y ~ x + ar(time = t, gr = g, cov = TRUE), sigma ~ x2), data = d)
  } else frm(fo[[m]], data = d)
  b <- B[[m]]
  ac <- fit$frame$autocor[[1]]
  bd <- b$betad
  bd[["sigma_(Intercept)"]] <- bd[["sigma_(Intercept)"]] -
    frmtmb:::autocor_log_marg(b$thetaac, ac)
  par <- c(b$beta, bd, b$thetaac)
  stopifnot(length(par) == length(fit$obj$par))
  nll <- fit$obj$fn(unname(par))
  cat(sprintf("%-5s r14 logLik %.17g | lane objective at r14's point %.17g | rel diff %.2e\n",
              m, b$logLik, -nll, abs(-nll - b$logLik) / abs(b$logLik)))
}
