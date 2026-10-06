# Reviewer: why does `a * b` (intercepts) get fixes' flat warning on the
# release build and the SE warning on the trial merge?
#   Rscript dev/nanse-rev-ab.R merge|release
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  release = "C:/Users/adf44/source/r/rellib-r6",
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("LIB", find.package("frmtmb"), "\n")
suppressMessages(trace("nl_flat_message", where = ns, print = FALSE,
  exit = quote(cat("  [nl_flat_message] returned",
                   if (is.null(returnValue())) "NULL" else "a message",
                   " opt$par", signif(opt$par, 6), "\n"))))
if (exists("se_check", ns)) {
  suppressMessages(trace("se_check", where = ns, print = FALSE,
    tracer = quote(cat("  [se_check] explained",
                       format(fit$cache$se_explained), " par",
                       signif(fit$opt$par, 6), "\n"))))
}
set.seed(955)
n <- 60
dd <- data.frame(x = rnorm(n), z = rnorm(n))
dd$y <- 3 + 0.5 * dd$x + 0.3 * dd$z + rnorm(n, 0, 0.4)
w <- character()
fit <- withCallingHandlers(
  frm(bf(y ~ a * b, a ~ 1, b ~ 1, nl = TRUE), data = dd,
      start = list(beta = c(1, 1))),
  warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
cat("code", fit$opt$convergence, "par", fit$opt$par, "\n")
for (x in w) cat("warn:", substr(x, 1, 300), "\n")
