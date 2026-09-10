# Choosing the guard on the geometric-tail correction.
#
# nss-04 found two constructions where an unguarded correction is much
# WORSE than truncating: a state with no steady state at all (an AUC
# compartment), and a system that has already converged into the
# integrator's noise, where the ratio is formed from noise. This sweeps
# the candidate guards over every case, good and bad, and reports the
# worst row of each, because a rule is only as good as its worst row.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-05-guard.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")

ATOL <- 1e-8
RTOL <- 1e-8
NSS <- 20L

run_in <- function(dyn, p, n_state, dose_state, amt, ii, n) {
  y <- numeric(n_state)
  ys <- matrix(0, n + 1L, n_state)
  for (k in seq_len(n)) {
    y[dose_state] <- y[dose_state] + amt
    s <- deSolve::ode(y, c(0, ii), dyn, p, method = "lsoda",
                      atol = ATOL, rtol = RTOL)
    y <- as.numeric(s[2L, 1L + seq_len(n_state)])
    ys[k + 1L, ] <- y
  }
  ys
}

# The rules. Each takes the last three iterates and returns a state.
# `damp` is Tikhonov damping by the integrator's own noise floor: a
# difference no bigger than the tolerance carries no information about
# the ratio, so it must not be allowed to set one.
mk_rule <- function(mode = c("each", "joint"), cap = 0.999, damp = 0) {
  mode <- match.arg(mode)
  function(y2, y1, y0) {
    d1 <- y1 - y0
    d2 <- y2 - y1
    del <- (damp * (ATOL + RTOL * abs(y2)))^2
    r <- if (mode == "each") (d1 * d2) / (d1 * d1 + del) else
      sum(d1 * d2) / (sum(d1 * d1) + sum(del))
    r <- pmin(pmax(r, 0), cap)
    y2 + d2 * r / (1 - r)
  }
}
rules <- list(
  plain            = function(y2, y1, y0) y2,
  joint            = mk_rule("joint", 0.999, 0),
  each             = mk_rule("each",  0.999, 0),
  joint_damp       = mk_rule("joint", 0.999, 4),
  each_damp        = mk_rule("each",  0.999, 4),
  each_damp_c99    = mk_rule("each",  0.99,  4),
  each_damp_c95    = mk_rule("each",  0.95,  4),
  each_damp16      = mk_rule("each",  0.999, 16))

cases <- list()
add <- function(tag, dyn, p, n_state, amt, ii, nref, keep = NULL) {
  ys <- run_in(dyn, p, n_state, 1L, amt, ii, nref)
  cases[[tag]] <<- list(ys = ys, ref = ys[nref + 1L, ],
                        keep = if (is.null(keep)) seq_len(n_state) else keep)
}

two_oral <- function(t, y, p) {
  list(c(-p$ka * y[1L],
         p$ka * y[1L] - (p$ke + p$k12) * y[2L] + p$k21 * y[3L],
         p$k12 * y[2L] - p$k21 * y[3L]))
}
one_oral <- function(t, y, p) list(c(-p$ka * y[1L],
                                     p$ka * y[1L] - p$ke * y[2L]))
mm <- function(t, y, p) {
  cc <- y[2L] / p$V
  list(c(-p$ka * y[1L], p$ka * y[1L] - p$Vmax * cc / (p$Km + cc)))
}
auc <- function(t, y, p) list(c(-p$ka * y[1L],
                                p$ka * y[1L] - p$ke * y[2L],
                                y[2L] / p$V))
osc <- function(t, y, p) list(c(y[2L], -p$w2 * y[1L] - p$z * y[2L]))
twoexp <- function(t, y, p) list(c(-p$k1 * y[1L],
                                   -p$k2 * y[2L] + y[1L]))

# the compartment grid nss-03 measured, where the correction must win
add("2cmt t1/2 23 ii 8",   two_oral, list(ka=1.1,ke=0.2,k12=0.4,k21=0.1),   3L, 100, 8,   4000L)
add("2cmt t1/2 107 ii 24", two_oral, list(ka=1.0,ke=0.15,k12=0.3,k21=0.02), 3L, 100, 24,  4000L)
add("2cmt t1/2 107 ii 12", two_oral, list(ka=1.0,ke=0.15,k12=0.3,k21=0.02), 3L, 100, 12,  6000L)
add("2cmt t1/2 265 ii 24", two_oral, list(ka=1.0,ke=0.1,k12=0.2,k21=0.008), 3L, 100, 24,  6000L)
add("2cmt t1/2 670 ii 24", two_oral, list(ka=1.0,ke=0.08,k12=0.15,k21=0.003),3L,100, 24, 12000L)
# the easy cases, where the correction must not do harm
add("1cmt ke*ii 2.4",  one_oral, list(ka=1.0,ke=0.1),  2L, 100, 24, 600L)
add("1cmt ke*ii 12",   one_oral, list(ka=3,  ke=0.5),  2L, 100, 24, 400L)
add("1cmt ke*ii 24",   one_oral, list(ka=3,  ke=1),    2L, 100, 24, 400L)
add("1cmt ke*ii 48",   one_oral, list(ka=3,  ke=2),    2L, 100, 24, 400L)
# nonlinear
add("MM Km 0.5", mm, list(ka=1,Vmax=8,Km=0.5,V=10), 2L, 100, 24, 4000L)
add("MM Km 5",   mm, list(ka=1,Vmax=8,Km=5,  V=10), 2L, 100, 24, 4000L)
add("MM Km 50",  mm, list(ka=1,Vmax=8,Km=50, V=10), 2L, 100, 24, 4000L)
# a state that never settles; only the two PK states are scored
add("AUC state", auc, list(ka=1,ke=0.1,V=10), 3L, 100, 24, 400L, keep = 1:2)
# oscillatory
add("osc ii 1", osc, list(w2=1,z=0.1), 2L, 1, 1, 8000L)
add("osc ii 3", osc, list(w2=1,z=0.1), 2L, 1, 3, 8000L)
add("osc ii 6", osc, list(w2=1,z=0.1), 2L, 1, 6, 8000L)
# nearly equal modes
add("2 modes 1.5x",   twoexp, list(k1=0.02,k2=0.03),   2L, 100, 24, 12000L)
add("2 modes 1.05x",  twoexp, list(k1=0.02,k2=0.021),  2L, 100, 24, 12000L)
add("2 modes 1.005x", twoexp, list(k1=0.02,k2=0.0201), 2L, 100, 24, 12000L)

score <- function(rule, cs) {
  ys <- cs$ys
  y <- rule(ys[NSS + 1L, ], ys[NSS, ], ys[NSS - 1L, ])
  k <- cs$keep
  max(abs(y[k] - cs$ref[k])) / max(abs(cs$ref[k]))
}

nm <- names(rules)
cat(sprintf("%-22s", "case"))
for (r in nm) cat(sprintf("%12s", r))
cat("\n")
tab <- matrix(NA_real_, length(cases), length(rules),
              dimnames = list(names(cases), nm))
for (i in seq_along(cases)) {
  cat(sprintf("%-22s", names(cases)[i]))
  for (j in seq_along(rules)) {
    v <- score(rules[[j]], cases[[i]])
    tab[i, j] <- v
    cat(sprintf("%12.3e", v))
  }
  cat("\n")
}
cat(sprintf("%-22s", "WORST ROW"))
for (j in seq_along(rules)) cat(sprintf("%12.3e", max(tab[, j])))
cat("\n")
cat(sprintf("%-22s", "worst vs plain"))
for (j in seq_along(rules))
  cat(sprintf("%12.3g", max(tab[, j] / tab[, "plain"])))
cat("\n")
cat(sprintf("%-22s", "best vs plain"))
for (j in seq_along(rules))
  cat(sprintf("%12.3g", min(tab[, j] / tab[, "plain"])))
cat("\n")
saveRDS(tab, "C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/nss-05-guard.rds")
