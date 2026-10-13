# Lane ci691: what the SE check reads on the xbeta precise-response fits
# of test-xbeta-zibb-hurdle-cum.R ("xbeta fits a precise response"),
# which lost every standard error on the ubuntu and macOS runners.
#   Rscript dev/ci691-xbeta-probe.R <lib>
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r7"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r7",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
cat("BLAS:", extSoftVersion()[["BLAS"]], La_library(), "\n")
options(digits = 17)
for (phi in c(1e4, 2e4)) {
  set.seed(4200)
  z <- stats::rbeta(1000, 0.047 * phi, 0.953 * phi)
  d <- data.frame(y = pmin(pmax(1.1 * z - 0.05, 0), 1))
  ws <- character(0)
  fit <- withCallingHandlers(frm(y ~ 1, family = xbeta(), data = d),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  cat("\n==== phi", phi, "\n")
  cat("warnings:", length(ws), "\n")
  for (w in ws) cat("  ", substr(w, 1, 160), "\n")
  o <- fit$opt
  cat("convergence", o$convergence, o$message, "iter", o$iterations,
      "evals", paste(o$evaluations, collapse = " "), "\n")
  cat("par", format(o$par, digits = 17), "\n")
  cat("objective", format(o$objective, digits = 17), "\n")
  g <- fit$obj$gr(o$par)
  cat("gradient", format(g, digits = 6), "\n")
  hx <- fit$obj$he(o$par)
  cat("he(par):\n"); print(hx)
  cat("eigen he:", format(eigen(hx, symmetric = TRUE)$values, digits = 6),
      "\n")
  hh <- frmtmb:::fit_outer_hessian(fit)
  cat("fd H:\n"); print(hh$H)
  cat("fd E max", format(max(hh$E), digits = 6), "\n")
  cat("par_units", format(fit$par_units %||% NA), "\n")
  cat("fd H finite", all(is.finite(hh$H)), "\n")
  cat("vcov se", format(sqrt(diag(vcov(fit, full = TRUE))), digits = 8), "\n")
  cat("se_explained:\n"); str(fit$cache$se_explained)
  ch <- frmtmb:::cov_from_hessian(fit, hh$H, hh$E)
  cat("tier", ch$tier, "lost", paste(names(ch$lost), ch$lost), "\n")
  if (!is.null(ch$ev)) cat("ev", format(ch$ev, digits = 6), "\n")
}
