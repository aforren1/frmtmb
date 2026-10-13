# Lane optrap, punch round 1: test-edgecases.R's constant response
# (y = 2, y ~ x + (1 | g)): its unit and where the fit stops.
#   Rscript dev/optrap-p1-flat.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(159)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
dd$y <- 2
w <- character()
f <- withCallingHandlers(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
u <- tryCatch(frmtmb:::predictor_units(f$frame), error = function(e) "none")
cat("lib", dirname(find.package("frmtmb")), "| units", deparse(u),
    "| objective", f$opt$objective, "code", f$opt$convergence,
    f$opt$message, "| sigma", stats::sigma(f), "| warnings", length(w),
    "\n")
