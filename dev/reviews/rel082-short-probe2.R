# Short check of 0.82.0 after the fix pass, probe 2: boundfam's (b), the
# frmtmb_limit_fit subclass in nl_fit_off_saddle(). A nonlinear student()
# fit on gaussian data, where nu's limit is real and the restart cannot
# help: does it now restart (cost), and is the limit message said once?
# usage: rel082-short-probe2.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")), "\n")
run <- function(label, expr) {
  cls <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(w) {
      cls <<- c(cls, paste0("W:", class(w)[1]))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      cls <<- c(cls, paste0("M:", class(m)[1]))
      invokeRestart("muffleMessage")
    })
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(f, "error")) {
    cat(label, "ERROR", conditionMessage(f), "\n")
    return(invisible())
  }
  cat(sprintf("%-30s ll %10.4f secs %5.2f | %s\n", label,
              as.numeric(logLik(f)), el, paste(cls, collapse = " ")))
}
for (s in 1:6) {
  set.seed(700 + s)
  n <- 150
  x <- stats::runif(n, 0, 3)
  y <- 2 * exp(0.4 * x) + stats::rnorm(n, 0, 0.5)
  d <- data.frame(x = x, y = y)
  run(sprintf("student exp, gaussian s%d", s),
      frm(bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), data = d,
          family = student()))
  d$k <- stats::rpois(n, 2 * exp(0.4 * x))
  run(sprintf("negbinomial exp, poisson s%d", s),
      frm(bf(k ~ a * exp(b * x), a + b ~ 1, nl = TRUE), data = d,
          family = negbinomial(link = "identity")))
}
