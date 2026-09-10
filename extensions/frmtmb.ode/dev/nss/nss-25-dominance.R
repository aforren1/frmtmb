# NIT 2: the durable ground for rejecting `n_ss = "auto"` is DOMINANCE,
# not "worse than the status quo", which a floor at 20 repairs.
#
# Both arms on the same design at the same tolerances: what a floored
# "auto" would buy for what it costs, against what the correction buys
# for what it costs. Counted, not timed.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-25-dominance.R
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
nss_report_env()
two <- function(t, y, p) list(c(-p[4L] * y[1L],
                                p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] +
                                  p[3L] * y[3L],
                                p[2L] * y[2L] - p[3L] * y[3L]))
P <- list(0.15, 0.3, 0.02, 1.0)
tt <- seq(0, 24, length.out = 25)
ev <- data.frame(time = 0, state = 1L, value = 100, ii = 24, ss = TRUE)
exact <- as.numeric(frm_lincmt(
  parms = list(ke = 0.15, k12 = 0.3, k21 = 0.02, ka = 1.0, V = 10),
  times = tt, ncmt = 2, depot = TRUE,
  events = data.frame(time = 0, state = "depot", value = 100, ii = 24,
                      addl = 0L, ss = TRUE)))
sc <- max(abs(exact))
err <- function(n, ext) {
  v <- suppressWarnings(frm_ode(two, init = list(0, 0, 0), times = tt,
                                parms = P, events = ev, output = 2L,
                                n_ss = n, ss_tol = Inf,
                                ss_extrapolate = ext, atol = 1e-8,
                                rtol = 1e-8))
  max(abs(as.numeric(v) / 10 - exact)) / sc
}
count <- function(n, ext) {
  cnt <- 0L
  ns <- asNamespace("deSolve")
  orig <- get("lsoda", envir = ns)
  unlockBinding("lsoda", ns)
  assign("lsoda", function(...) { cnt <<- cnt + 1L; orig(...) },
         envir = ns)
  on.exit({ assign("lsoda", orig, envir = ns)
            lockBinding("lsoda", ns) }, add = TRUE)
  invisible(err(n, ext))
  cnt
}
cat("\n2 cmt oral, terminal half-life 107.1 h, ii 24, atol=rtol=1e-8\n")
cat(sprintf("%-34s %12s %10s\n", "arm", "error", "lsoda"))
for (z in list(list("today, n_ss = 20, truncated", 20L, FALSE),
               list("floored auto, n_ss = 134, truncated", 134L, FALSE),
               list("floored auto, n_ss = 266, truncated", 266L, FALSE),
               list("shipped, n_ss = 20, tail summed", 20L, TRUE))) {
  cat(sprintf("%-34s %12.3e %10d\n", z[[1]], err(z[[2]], z[[3]]),
              count(z[[2]], z[[3]])))
}
