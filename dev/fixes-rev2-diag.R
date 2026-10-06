# Reviewer of lane fixes, re-check: what diagnose() says about a smooth
# fit now that smooth_fx_units() gives it par_units, without autoscale.
#   Rscript dev/fixes-rev2-diag.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
set.seed(3)
d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, verbose = FALSE))
d$w <- runif(200, 0, 1e4)
fit <- suppressWarnings(frm(bf(y ~ s(x1) + s(x2) + w), data = d))
cat("par_units NULL:", is.null(fit$par_units), " control autoscale:",
    format(fit$control$autoscale), "\n")
diagnose(fit)
