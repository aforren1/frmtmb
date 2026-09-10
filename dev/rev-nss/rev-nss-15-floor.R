# REVIEW of lane nss round 1, attack 6: is 2.618e-09 really the
# integrator's floor?
#
# `dev/nss/nss-25-dominance.R` argues that floored "auto" at
# `n_ss` = 134 reaches 2.854e-09 for 135 solves while the correction
# reaches 2.623e-09 for 21, and calls 2.618e-09 at `n_ss` = 266 the
# integrator's floor. That floor is what makes the comparison fair: if
# the error were still falling with `n_ss`, "auto" would only need a
# bigger `n_ss` and the dominance claim would be about the constant
# rather than about the design.
#
# A floor is only a floor if it MOVES WITH THE INSTRUMENT. So the same
# grid is run at four tolerances. If 2.6e-09 is the integrator's own
# error it falls when `atol` and `rtol` fall; if it is the run-in's, it
# does not.
#
# Script path: dev/rev-nss/rev-nss-15-floor.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
rev_env()

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
err <- function(n, ext, tol) {
  v <- suppressWarnings(frm_ode(two, init = list(0, 0, 0), times = tt,
                                parms = P, events = ev, output = 2L,
                                n_ss = n, ss_tol = Inf,
                                ss_extrapolate = ext, atol = tol,
                                rtol = tol))
  max(abs(as.numeric(v) / 10 - exact)) / sc
}

cat("\n== the error against n_ss, at four tolerances,",
    "ss_extrapolate = FALSE ==\n")
cat("   solves per group is n_ss + 1\n\n")
cat(sprintf("%6s %8s", "n_ss", "solves"))
tols <- c(1e-6, 1e-8, 1e-10, 1e-12)
for (tl in tols) cat(sprintf(" %12s", format(tl)))
cat("\n")
for (n in c(20L, 60L, 134L, 200L, 266L, 400L, 600L)) {
  cat(sprintf("%6d %8d", n, n + 1L))
  for (tl in tols) cat(sprintf(" %12.3e", err(n, FALSE, tl)))
  cat("\n")
}
cat("\n== the same for the correction, ss_extrapolate = TRUE ==\n")
cat(sprintf("%6s %8s", "n_ss", "solves"))
for (tl in tols) cat(sprintf(" %12s", format(tl)))
cat("\n")
for (n in c(20L, 60L, 134L, 266L)) {
  cat(sprintf("%6d %8d", n, n + 1L))
  for (tl in tols) cat(sprintf(" %12.3e", err(n, TRUE, tl)))
  cat("\n")
}
cat("\nIf the plateau at large n_ss falls when the tolerance falls, it",
    "is\nthe integrator's floor and the lane's comparison is fair. If",
    "it does\nnot move, it is not a floor and the dominance claim",
    "rests on a\nconstant that could be raised.\n")
cat("\ndone\n")
