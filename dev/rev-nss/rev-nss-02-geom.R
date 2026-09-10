# REVIEW of lane nss, attack 2: where does `d * r / (1 - r)` stop being
# the tail?
#
# The lane discloses a constructed oscillator where extrapolation is
# 1.39x worse than truncation and says the field does not write it. The
# question the brief asks is whether an ORDINARY population PK schedule
# can put it there. This script answers it on EXACT cycle maps, so
# nothing is confounded with solver error.
#
# The mechanism to hunt. Write the cycle-start error as a sum of modes,
#   y_inf - y_n = sum_i C_i a_i^n,   a_i = exp(mu_i * ii),
# and the successive difference as
#   d_n = y_n - y_{n-1} = sum_i C_i a_i^{n-2} (1 - a_i) * a_i^0 ... ,
# so the ratio the run-in reads is a weighted mean of the a_i with
# weights u_i = C_i a_i^{n-2} (1 - a_i). When every u_i has the SAME
# sign that mean lies between the slowest and the fastest a_i, so the
# ratio UNDERSTATES the dominant mode and the correction undershoots:
# it can only improve. When two u_i have OPPOSITE signs the mean is an
# extrapolation outside the range and the ratio OVERSTATES it, without
# bound as the two weights approach cancellation.
#
# Opposite signs are not exotic in pharmacokinetics: the central
# compartment of any oral model carries the absorption mode with the
# sign opposite to the disposition modes, which is what makes a
# concentration rise before it falls. So the overshoot regime is
# ka near lambda_z, and ka near lambda_z is flip-flop kinetics, which
# extended-release and depot formulations are written to produce.
#
# The sweep below is ordinary therapeutics: one and two compartments,
# oral, half-lives 1.4 h to 693 h, ka 0.02 to 3 /h, ii 8, 12, 24, 48 h,
# at the shipped n_ss = 20 and the shipped atol = rtol = 1e-8.
#
# Script path: dev/rev-nss/rev-nss-02-geom.R
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages(library(frmtmb.ode))
ext <- frmtmb.ode:::ode_ss_extrapolate

# Exact cycle map of a linear mammillary model with a depot.
# States: depot, central, peripheral (2 cmt only).
# The run-in the package performs records the TROUGH: the cycle-start
# state BEFORE that cycle's dose, so y_k = A (y_{k-1} + e).
build <- function(ke, ka, k12 = 0, k21 = 0, ncmt = 1, ii = 24,
                  dose = 100) {
  M <- if (ncmt == 1)
    matrix(c(-ka, 0, ka, -ke), 2, 2)
  else
    matrix(c(-ka, 0, 0,
             ka, -(ke + k12), k21,
             0, k12, -k21), 3, 3, byrow = TRUE)
  eg <- eigen(M)
  A <- Re(eg$vectors %*% diag(exp(eg$values * ii)) %*%
            solve(eg$vectors))
  e <- c(dose, rep(0, nrow(M) - 1L))
  list(A = A, Ae = as.numeric(A %*% e), e = e,
       yinf = as.numeric(solve(diag(nrow(A)) - A, A %*% e)),
       lam = -max(Re(eg$values[Re(eg$values) < -1e-14])))
}

# The run-in exactly as ode_run_in() performs it.
runin <- function(m, n_ss) {
  y <- numeric(nrow(m$A))
  keep <- vector("list", n_ss + 1L)
  keep[[1L]] <- y
  for (k in seq_len(n_ss)) {
    y <- as.numeric(m$A %*% (y + m$e))
    keep[[k + 1L]] <- y
  }
  # ode_run_in() hands ode_ss_extrapolate() the last three cycle-start
  # states, y_{n-2}, y_{n-1}, y_n, with y_0 = 0 the first
  keep[seq(n_ss - 1L, n_ss + 1L)]
}

one <- function(ke, ka, k12 = 0, k21 = 0, ncmt = 1, ii = 24,
                n_ss = 20L, atol = 1e-8, rtol = 1e-8) {
  m <- build(ke, ka, k12, k21, ncmt, ii)
  kp <- runin(m, n_ss)
  y <- kp[[3L]]
  ex <- ext(kp[[1L]], kp[[2L]], y, atol, rtol)
  sc <- max(1e-12, max(abs(m$yinf)))
  list(trunc = max(abs(y - m$yinf)) / sc,
       extrap = max(abs(as.numeric(ex$y) - m$yinf)) / sc,
       et = abs(y - m$yinf) / sc,
       ee = abs(as.numeric(ex$y) - m$yinf) / sc,
       r = as.numeric(ex$r), lam = m$lam,
       a_dom = exp(-m$lam * ii))
}

## ---- A. the mechanism, isolated: 1 cmt oral, ka swept past ke ------
cat("\n== A. one-compartment oral, ke = 0.03 (t1/2 23.1 h), ii = 24,",
    "n_ss = 20 ==\n")
cat("   ratio > 1 means the shipped default is WORSE than truncation\n\n")
cat(sprintf("%9s %9s %11s %11s %8s %9s\n", "ka", "ka/ke", "truncated",
            "extrapolated", "ratio", "r(central)"))
ke <- 0.03
worst <- list(ratio = 0)
for (ka in ke * c(0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 1.001, 1.01, 1.02,
                  1.05, 1.1, 1.2, 1.5, 2, 5, 20)) {
  o <- one(ke, ka, ii = 24)
  rt <- o$extrap / o$trunc
  cat(sprintf("%9.5f %9.3f %11.3e %11.3e %8.2f %9.5f\n", ka, ka / ke,
              o$trunc, o$extrap, rt, o$r[2L]))
  if (is.finite(rt) && rt > worst$ratio)
    worst <- list(ratio = rt, ka = ka, o = o)
}
cat(sprintf("\nworst ratio in A: %.2f at ka = %.5f\n", worst$ratio,
            worst$ka %||% NA))

## ---- B. an ordinary-therapeutics sweep -----------------------------
cat("\n== B. sweep of ordinary schedules, n_ss = 20 ==\n")
grid <- expand.grid(
  thalf = c(1.4, 3.5, 6, 12, 23, 48, 107, 265, 693),
  karat = c(0.2, 0.5, 0.8, 0.95, 1.0, 1.05, 1.2, 2, 5, 20, 60),
  ii = c(8, 12, 24, 48), ncmt = c(1, 2))
rows <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  ke <- log(2) / g$thalf
  if (g$ncmt == 1) {
    k12 <- 0; k21 <- 0; lz <- ke
  } else {
    # a peripheral compartment that actually slows the terminal phase
    k12 <- 2 * ke; k21 <- ke / 3
    b <- ke + k12 + k21
    lz <- (b - sqrt(b * b - 4 * ke * k21)) / 2
  }
  ka <- g$karat * lz
  if (ka > 5 || ka < 1e-4) next
  o <- tryCatch(one(ke, ka, k12, k21, g$ncmt, g$ii),
                error = function(e) NULL)
  if (is.null(o)) next
  rows[[length(rows) + 1L]] <- data.frame(
    ncmt = g$ncmt, thalf = g$thalf, ii = g$ii, karat = g$karat,
    ka = ka, lz = lz, trunc = o$trunc, extrap = o$extrap,
    ratio = o$extrap / o$trunc)
}
d <- do.call(rbind, rows)
d <- d[is.finite(d$ratio), ]
cat("rows:", nrow(d), "\n")
cat("rows where extrapolation is WORSE than truncation:",
    sum(d$ratio > 1), "\n")
cat("rows where it is worse AND the truncated error is above 1e-6:",
    sum(d$ratio > 1 & d$trunc > 1e-6), "\n\n")
o <- d[order(-d$ratio), ]
cat("worst 15 by ratio:\n")
print(format(head(o, 15), digits = 3), row.names = FALSE)
cat("\nworst 10 by ABSOLUTE extrapolated error:\n")
print(format(head(d[order(-d$extrap), ], 10), digits = 3),
      row.names = FALSE)
cat("\nquantiles of ratio:\n")
print(signif(quantile(d$ratio, c(0, .5, .9, .99, 1)), 3))
cat("\ndone\n")
