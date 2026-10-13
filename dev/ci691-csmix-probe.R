# Lane ci691: test-optima.R:163 asserts that the cs() ordinal mixture's
# reported objective is the objective at the reported estimates, bit for
# bit. On the ubuntu runners the two differed by 1.2e-13
# (374.16747089090501 against 374.16747089090489). This prints, for the
# test's seeds 3 and 12, where the reported objective came from.
#   Rscript dev/ci691-csmix-probe.R <lib>
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r7"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r7",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
cat("LAPACK:", La_library(), "\n")
optima_csmix_data <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
for (s in c(3, 12)) {
  ws <- character(0)
  fit <- withCallingHandlers(
    frm(bf(y ~ x + cs(z)), family = mixture(cumulative(), sratio()),
        data = optima_csmix_data(s)),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  o <- fit$opt
  v <- fit$obj$fn(o$par)
  cat(sprintf("\nseed %d code %d (%s) iter %s\n", s, o$convergence,
              o$message, format(o$iterations)))
  cat(sprintf("objective %.17g\nfn(par)   %.17g\nidentical %s  diff %.3g\n",
              o$objective, v, identical(v, o$objective), v - o$objective))
  cat("rejected_last", format(o$rejected_last), " nonfinite_trials",
      format(o$nonfinite_trials), " par_units",
      if (is.null(fit$par_units)) "NULL" else
        paste(format(fit$par_units, digits = 4), collapse = " "), "\n")
  cat("stationary_escape", format(o$stationary_escape), " mo_search",
      format(o$mo_search), "\n")
  cat("names(opt):", paste(names(o), collapse = " "), "\n")
  # repeat evaluations: is fn deterministic at this point?
  vv <- vapply(1:5, function(i) fit$obj$fn(o$par), 0)
  cat("repeat fn identical:", length(unique(vv)) == 1L, "\n")
  cat("warnings:", length(ws), "\n")
  for (w in ws) cat("  ", substr(w, 1, 120), "\n")
}
