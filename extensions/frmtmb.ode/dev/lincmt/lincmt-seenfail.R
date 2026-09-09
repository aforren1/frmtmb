# The gradient test, run against the UNFIXED primitive, to record the
# failure it pins. The only change is lincmt_phi(): the offset-only
# spelling this file started with.
.libPaths(c("C:/Users/adf44/source/r/lincmt-lib","C:/Users/adf44/source/r/reflib-r2","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB);library(testthat);library(frmtmb);library(frmtmb.ode)})
ns <- asNamespace("frmtmb.ode")
tt <- c(0, 0.5, 1, 2, 4, 8, 12, 18, 24, 30, 36, 48)
ev <- data.frame(time = 0, state = "depot", value = 100, ii = 12, addl = 3L)
f <- function(th) sum(frm_lincmt(parms = list(ka = exp(th[1]), ke = exp(th[2]), V = exp(th[3])),
                                 times = tt, ncmt = 1, depot = TRUE, events = ev))
th <- c(log(0.2), log(0.2), log(10))
near <- c(log(0.2 * (1 + 1e-4)), log(0.2), log(10))
report <- function(tag) {
  g <- as.numeric(MakeTape(f, th)$jacobian(th))
  gn <- as.numeric(MakeTape(f, near)$jacobian(near))
  cat(tag, "\n  value ", format(f(th), digits = 12), "\n")
  cat("  grad at ka == ke  ", paste(format(g, digits = 6), collapse = "  "), "\n")
  cat("  grad just off it  ", paste(format(gn, digits = 6), collapse = "  "), "\n")
  ok <- all(is.finite(g)) && max(abs(g)) < 2 * max(abs(gn))
  cat("  assertion max|g| < 2 * max|g_near|:", if (ok) "PASS" else "FAIL",
      " ratio", format(max(abs(g)) / max(abs(gn))), "\n")
}
report("FIXED (shipped lincmt_phi)")
old_phi <- function(x) { z <- x + 1e-300; -expm1(-z) / z }
environment(old_phi) <- ns
unlockBinding("lincmt_phi", ns)
assign("lincmt_phi", old_phi, envir = ns)
report("UNFIXED (offset only, no series blend)")
