# Lane gpstop: the new test's assertions on 20 roundings of the data, y *
# (1 + 2e-15 k), k = 1 to 20 (lane-env: a test must not depend on one
# rounding). Rows and starts as in test-gp-start.R, "a multiplicative body
# is linearized again at the fit's optimum". Args: arm.
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-lib.R")
gs_libs(a[1])
source(file.path(GS_WT, "dev/review-fitedge-designs.R"))
ll <- function(...) {
  q <- fe_fit(...)
  if (inherits(q$fit, "error")) NA_real_ else as.numeric(logLik(q$fit))
}
cases <- list(
  E = list(d = rv_mk("E", 2, 120L, 2), f = rv_spell$E$nl,
           th = c(log(1.5), log(0.3))),
  L = list(d = rv_mk("L", 2, 120L, 2), f = rv_spell$L$nl,
           th = c(log(0.3), log(0.1))))
fails <- 0L
for (nm in names(cases)) for (k in 1:20) {
  cs <- cases[[nm]]
  d <- cs$d
  d$y <- d$y * (1 + 2e-15 * k)
  ref <- ll(cs$f, data = d, start = list(theta = cs$th))
  def <- ll(cs$f, data = d)
  ok <- isTRUE(def > ref - 1e-5 * abs(ref))
  fails <- fails + !ok
  cat(sprintf("ULP %s k %2d default %.6f ref %.6f diff %.3e %s\n", nm, k, def,
              ref, def - ref, if (ok) "ok" else "FAIL"))
  flush.console()
}
cat(sprintf("ULPSUM rows %d fail %d\n", 40L, fails))
