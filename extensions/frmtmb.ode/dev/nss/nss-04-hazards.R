# Where a geometric-tail extrapolation of the run-in can go wrong.
#
# Every row here is a system frm_ode() accepts today. The question is
# not whether extrapolation helps on a compartment model (nss-03 says
# it does) but whether it turns a merely-short answer into a wild one
# on a system that is nonlinear, that has a state with no steady state
# at all, that oscillates, or that has already converged into the
# integrator's noise.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-04-hazards.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")

ATOL <- 1e-8
RTOL <- 1e-8

# One run-in, exactly as ode_run_in() performs it for a bolus.
run_in <- function(dyn, p, n_state, dose_state, amt, ii, n,
                   atol = ATOL, rtol = RTOL) {
  y <- numeric(n_state)
  ys <- matrix(0, n + 1L, n_state)
  for (k in seq_len(n)) {
    y[dose_state] <- y[dose_state] + amt
    s <- deSolve::ode(y, c(0, ii), dyn, p, method = "lsoda",
                      atol = atol, rtol = rtol)
    y <- as.numeric(s[2L, 1L + seq_len(n_state)])
    ys[k + 1L, ] <- y
  }
  ys
}

r_dot <- function(d1, d2) sum(d1 * d2) / sum(d1 * d1)
# the same estimate formed one state at a time; eps keeps a state whose
# difference has already vanished from dividing by zero
r_each <- function(d1, d2) (d1 * d2) / (d1 * d1 + 1e-300)

clamp <- function(r, lo = 0, hi = 0.999) pmin(pmax(r, lo), hi)

ext_joint <- function(ys, n, cap = TRUE) {
  d <- diff(ys)
  r <- r_dot(d[n - 1L, ], d[n, ])
  if (cap) r <- clamp(r)
  ys[n + 1L, ] + d[n, ] * r / (1 - r)
}
ext_each <- function(ys, n, cap = TRUE) {
  d <- diff(ys)
  r <- r_each(d[n - 1L, ], d[n, ])
  if (cap) r <- clamp(r)
  ys[n + 1L, ] + d[n, ] * r / (1 - r)
}
rel <- function(y, ref) max(abs(y - ref)) / max(abs(ref))

report <- function(tag, ys, ref, n = 20L) {
  cat(sprintf("%-34s plain %10.3e | joint %10.3e | each %10.3e |",
              tag, rel(ys[n + 1L, ], ref),
              rel(ext_joint(ys, n), ref), rel(ext_each(ys, n), ref)))
  cat(sprintf(" uncapped joint %10.3e\n",
              rel(ext_joint(ys, n, cap = FALSE), ref)))
}

cat("=== H1. Michaelis-Menten elimination, a nonlinear system ===\n")
cat("Reference is the same run-in taken to 4000 cycles.\n")
mm <- function(t, y, p) {
  cc <- y[2L] / p$V
  list(c(-p$ka * y[1L],
         p$ka * y[1L] - p$Vmax * cc / (p$Km + cc)))
}
for (km in c(0.5, 5, 50)) {
  p <- list(ka = 1.0, Vmax = 8, Km = km, V = 10)
  ys <- run_in(mm, p, 2L, 1L, 100, 24, 4000L)
  ref <- ys[4001L, ]
  report(sprintf("MM  Km = %-5g", km), ys, ref)
}

cat("\n=== H2. a state with no steady state: an AUC compartment ===\n")
cat("dA3/dt = A2 / V. It grows without bound, so no run-in length is\n",
    "right for it. The question is what it does to the other two.\n",
    sep = "")
auc <- function(t, y, p) {
  list(c(-p$ka * y[1L],
         p$ka * y[1L] - p$ke * y[2L],
         y[2L] / p$V))
}
p <- list(ka = 1.0, ke = 0.1, V = 10)
ys <- run_in(auc, p, 3L, 1L, 100, 24, 400L)
ref3 <- ys[401L, ]
# the PK states alone have a limit; read the damage on those two only
refpk <- ys[401L, 1:2]
pk <- function(y) y[1:2]
cat(sprintf("%-34s plain %10.3e | joint %10.3e | each %10.3e |",
            "PK states, AUC state present",
            rel(pk(ys[21L, ]), refpk),
            rel(pk(ext_joint(ys, 20L)), refpk),
            rel(pk(ext_each(ys, 20L)), refpk)))
cat(sprintf(" uncapped joint %10.3e\n",
            rel(pk(ext_joint(ys, 20L, cap = FALSE)), refpk)))
ys2 <- run_in(function(t, y, p) list(c(-p$ka * y[1L],
                                       p$ka * y[1L] - p$ke * y[2L])),
              p, 2L, 1L, 100, 24, 400L)
cat(sprintf("%-34s plain %10.3e | joint %10.3e | each %10.3e\n",
            "the same two states, AUC absent",
            rel(ys2[21L, ], ys2[401L, ]),
            rel(ext_joint(ys2, 20L), ys2[401L, ]),
            rel(ext_each(ys2, 20L), ys2[401L, ])))

cat("\n=== H3. an oscillating system, complex cycle-map eigenvalues ===\n")
cat("A damped oscillator dosed every ii. Not pharmacokinetics, but\n",
    "frm_ode() takes an arbitrary dyn and will be handed one.\n",
    sep = "")
osc <- function(t, y, p) list(c(y[2L], -p$w2 * y[1L] - p$z * y[2L]))
for (ii in c(1.0, 3.0, 6.0)) {
  p <- list(w2 = 1, z = 0.1)
  ys <- run_in(osc, p, 2L, 1L, 1, ii, 4000L)
  ref <- ys[4001L, ]
  d <- diff(ys)
  cat(sprintf("  ii = %-4g r_dot = %8.4f  ", ii,
              r_dot(d[19L, ], d[20L, ])))
  report("", ys, ref)
}

cat("\n=== H4. already converged: the difference is solver noise ===\n")
cat("ke * ii large, so after three cycles y_n - y_(n-1) is below the\n",
    "integrator's own tolerance and the ratio is formed from noise.\n",
    sep = "")
one <- function(t, y, p) list(c(-p$ka * y[1L],
                                p$ka * y[1L] - p$ke * y[2L]))
for (ke in c(0.5, 1, 2)) {
  p <- list(ka = 3, ke = ke)
  ys <- run_in(one, p, 2L, 1L, 100, 24, 400L)
  ref <- ys[401L, ]
  d <- diff(ys)
  cat(sprintf("  ke*ii = %-5g r_dot(19,20) = %11.4e  ", ke * 24,
              r_dot(d[19L, ], d[20L, ])))
  report("", ys, ref)
}

cat("\n=== H5. two nearly equal modes ===\n")
cat("The extrapolation removes ONE mode. With two of the same size it\n",
    "removes their combination, and what is left is the other one.\n",
    sep = "")
two <- function(t, y, p) list(c(-p$k1 * y[1L], -p$k2 * y[2L] + y[1L]))
for (dd in c(0.5, 0.05, 0.005)) {
  p <- list(k1 = 0.02, k2 = 0.02 * (1 + dd))
  ys <- run_in(two, p, 2L, 1L, 100, 24, 8000L)
  ref <- ys[8001L, ]
  report(sprintf("k2/k1 - 1 = %-8g", dd), ys, ref)
}
