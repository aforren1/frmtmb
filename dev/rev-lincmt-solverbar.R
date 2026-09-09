# Is the 1e-11-ish disagreement on a multi-dose schedule the closed
# form's error or the solver's? Run frm_ode() at three tolerances and
# see which of the two moves.
#
# If the closed form is right, tightening the solver moves frm_ode()
# TOWARD it and the residual disagreement shrinks with the tolerance.
# If the closed form were wrong, the disagreement would stall at the
# closed form's own error.
#
# Script path: dev/rev-lincmt-solverbar.R.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")

dyn2 <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[[4]] * y[1],
         p[[4]] * y[1] - (p[[1]] + p[[2]]) * y[2] + p[[3]] * y[3],
         p[[2]] * y[2] - p[[3]] * y[3]))
}
P <- list(ke = 0.2, k12 = 0.4, k21 = 0.1, ka = 1.1, V = 10)
cases <- list(
  "two doses, conc, long tail" = list(
    times = c(0, 1, 4, 12, 30, 60, 90),
    ev = data.frame(time = c(0, 14), state = c("depot", "central"),
                    value = c(120, 60), method = "add",
                    duration = c(0, 1.5), ii = 0, addl = 0L,
                    ss = FALSE)),
  "ss plus addl, observed late" = list(
    times = c(0, 4, 8, 40, 80),
    ev = data.frame(time = c(0, 8), state = "depot",
                    value = c(100, 100), method = "add",
                    duration = 0, ii = c(8, 8), addl = c(0L, 5L),
                    ss = c(TRUE, FALSE))),
  "30-dose addl block" = list(
    times = c(0, 120, 232, 240),
    ev = data.frame(time = 0, state = "depot", value = 100,
                    method = "add", duration = 0, ii = 8, addl = 29L,
                    ss = FALSE)))

cat(sprintf("\n%-30s %12s %12s %12s %12s\n", "case", "ode 1e-9",
            "ode 1e-12", "ode 1e-14", "1e-12 vs 14"))
for (nm in names(cases)) {
  z <- cases[[nm]]
  a <- frm_lincmt(parms = P, times = z$times, ncmt = 2, depot = TRUE,
                  events = z$ev, n_ss = 20L)
  o <- lapply(c(1e-9, 1e-12, 1e-14), function(tol)
    frm_ode(dyn2, init = list(0, 0, 0), times = z$times,
            parms = list(P$ke, P$k12, P$k21, P$ka),
            states = c("depot", "central", "peripheral1"),
            output = "central", events = z$ev, n_ss = 20L,
            atol = tol, rtol = tol) / P$V)
  sc <- max(abs(o[[3L]]))
  d <- vapply(o, function(b) max(abs(a - b)) / sc, 0)
  spread <- max(abs(o[[2L]] - o[[3L]])) / sc
  cat(sprintf("%-30s %12.3e %12.3e %12.3e %12.3e\n", nm, d[[1L]],
              d[[2L]], d[[3L]], spread))
}
cat("\nA disagreement that falls with the solver's tolerance belongs",
    "to\nthe solver. One that stalls belongs to the closed form.\n")
