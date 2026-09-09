# The load-independent count that IS comparable between the two paths:
# how many numerical ODE solves frm_ode() performs per group per
# objective evaluation, against how many arithmetic tape nodes the
# closed form uses for the same group.
.libPaths(c("C:/Users/adf44/source/r/lincmt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(RTMB); library(frmtmb); library(frmtmb.ode)
})
ode_dyn <- function(t, y, p) {
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}
doses <- data.frame(time = c(0, 12), state = "depot", value = 100,
                    ii = c(12, 12), addl = c(0L, 12L),
                    ss = c(TRUE, FALSE))
tt <- 144 + c(0.5, 1, 2, 4, 6, 8, 10, 12)

count_solves <- function(n_ss) {
  n <- 0L
  f <- RTMBode::ode
  tmp <- new.env()
  assign("ode", function(...) { n <<- n + 1L; f(...) }, envir = tmp)
  # frm_ode() calls RTMBode::ode() by its fully qualified name, so the
  # counter is installed in the RTMBode namespace itself
  ns <- asNamespace("RTMBode")
  unlockBinding("ode", ns)
  assign("ode", get("ode", envir = tmp), envir = ns)
  on.exit({ assign("ode", f, envir = ns); lockBinding("ode", ns) })
  invisible(frm_ode(ode_dyn, init = list(0, 0), times = tt,
                    parms = list(1, 0.15, 20),
                    states = c("depot", "central"), output = "central",
                    events = doses, n_ss = n_ss))
  n
}
for (n_ss in c(20L, 1L)) {
  cat("frm_ode(n_ss =", n_ss, "): ", count_solves(n_ss),
      "RTMBode::ode() calls for ONE group's 8 observations\n")
}
ndf <- function(f) nrow(MakeTape(f, numeric(3))$data.frame())
for (nss in list(20L, Inf)) {
  n <- ndf(function(p) sum(frm_lincmt(
    parms = list(ka = exp(p[1]), ke = exp(p[2]), V = exp(p[3])),
    times = tt, ncmt = 1, depot = TRUE, events = doses,
    n_ss = nss)))
  cat("frm_lincmt(n_ss =", format(nss), "):", n,
      "arithmetic tape nodes for the same group\n")
}
