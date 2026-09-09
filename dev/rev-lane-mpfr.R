# frm_lincmt() against a 240-bit reference.
#
# The reference never forms an eigenvalue and never divides by a
# difference of rate constants: it exponentiates the rate matrix by
# scaling and squaring in Rmpfr, and integrates it with the augmented
# matrix [[M, b], [0, 0]] for an infusion. So it agrees with the closed
# form only if the closed form's algebra is right, and it keeps every
# digit where the closed form's textbook spelling would keep none.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src2.R")
suppressPackageStartupMessages(library(Rmpfr))
PREC <- 240L

mp <- function(x) mpfr(x, PREC)

expm_mp <- function(A) {
  n <- nrow(A)
  nrm <- max(as.numeric(apply(abs(A), 1, sum)))
  s <- max(0L, as.integer(ceiling(log2(max(nrm, 1e-300)))) + 6L)
  B <- A / mp(2)^s
  I <- mpfrArray(0, PREC, c(n, n))
  for (i in 1:n) I[i, i] <- mp(1)
  T <- I
  P <- I
  for (k in 1:45) {
    P <- P %*% B / mp(k)
    T <- T + P
  }
  for (k in seq_len(s)) T <- T %*% T
  T
}

# rate matrix, depot first when present
mk_mat <- function(ncmt, depot, p) {
  n <- ncmt + as.integer(depot)
  M <- mpfrArray(0, PREC, c(n, n))
  ic <- if (depot) 2L else 1L
  ke <- mp(p[["ke"]])
  if (depot) {
    ka <- mp(p[["ka"]])
    M[1, 1] <- -ka
    M[ic, 1] <- ka
  }
  out <- -ke
  if (ncmt >= 2L) {
    k12 <- mp(p[["k12"]]); k21 <- mp(p[["k21"]])
    out <- out - k12
    M[ic + 1L, ic] <- k12
    M[ic, ic + 1L] <- k21
    M[ic + 1L, ic + 1L] <- -k21
  }
  if (ncmt == 3L) {
    k13 <- mp(p[["k13"]]); k31 <- mp(p[["k31"]])
    out <- out - k13
    M[ic + 2L, ic] <- k13
    M[ic, ic + 2L] <- k31
    M[ic + 2L, ic + 2L] <- -k31
  }
  M[ic, ic] <- out
  M
}

# amount in `out` at lag u after a unit bolus into `cin`
bolus_mp <- function(M, u, cin, out) {
  if (u <= 0) return(if (cin == out) mp(1) else mp(0))
  E <- expm_mp(M * mp(u))
  E[out, cin]
}

# amount in `out` at lag u after a unit-AMOUNT infusion of length dur
inf_mp <- function(M, u, dur, cin, out) {
  n <- nrow(M)
  A <- mpfrArray(0, PREC, c(n + 1L, n + 1L))
  A[1:n, 1:n] <- M
  A[cin, n + 1L] <- mp(1)
  w <- min(u, dur)
  if (w <= 0) return(mp(0))
  E <- expm_mp(A * mp(w))
  v <- E[1:n, n + 1L] / mp(dur)
  if (u > dur) {
    F <- expm_mp(M * mp(u - dur))
    s <- mp(0)
    for (j in 1:n) s <- s + F[out, j] * v[j]
    return(s)
  }
  v[out]
}

resp_mp <- function(M, u, dur, cin, out) {
  if (dur > 0) inf_mp(M, u, dur, cin, out) else
    bolus_mp(M, u, cin, out)
}

# the reference for one schedule: a list of doses, each with amount,
# lag, compartment and duration, plus an optional steady-state block
ref_mp <- function(ncmt, depot, p, terms, out) {
  M <- mk_mat(ncmt, depot, p)
  n <- nrow(M)
  s <- mp(0)
  for (z in terms) {
    if (is.null(z$ii)) {
      s <- s + mp(z$amt) * resp_mp(M, z$u, z$dur, z$cin, out)
    } else {
      # one exponential for the first lag and one for the interval,
      # then repeated matrix-vector products: 400 cycles cost two
      # matrix exponentials rather than 400 of them
      E <- expm_mp(M * mp(z$ii))
      v <- expm_mp(M * mp(z$u))[, z$cin]
      acc <- mp(0)
      for (j in 0:(z$n - 1L)) {
        acc <- acc + v[out]
        w <- rep(mp(0), n)
        for (a in 1:n) {
          t <- mp(0)
          for (b in 1:n) t <- t + E[a, b] * v[b]
          w[a] <- t
        }
        v <- w
      }
      s <- s + mp(z$amt) * acc
    }
  }
  s
}

relerr <- function(a, b) as.numeric(abs(mp(a) - b) / abs(b))

# ---------------------------------------------------------------------
# 1. The absorption singularity: ka approaching ke, one dose.

cat("== ka -> ke, single oral dose, 1 compartment ==\n")
flush(stdout()); cat(sprintf("%12s %6s %16s %16s %10s\n", "ka - ke", "u", "frm_lincmt",
            "240-bit", "rel"))
for (d in c(0, 10^-(16:0))) {
  p <- list(ke = 0.2, ka = 0.2 + d)
  for (u in c(1, 24)) {
    a <- frm_lincmt(parms = c(p, list(V = 1)), times = u, ncmt = 1,
                    depot = TRUE, init = list(depot = 100),
                    output = "central")
    b <- ref_mp(1L, TRUE, p, list(list(amt = 100, u = u, dur = 0,
                                       cin = 1L)), 2L)
    flush(stdout()); cat(sprintf("%12.1e %6g %16.10e %16.10e %10.2e\n", d, u, a,
                as.numeric(b), relerr(a, b)))
  }
}

cat("\n== ka -> ke, steady state (exact limit), 1 compartment ==\n")
for (d in c(0, 10^-(12:0))) {
  p <- list(ke = 0.2, ka = 0.2 + d)
  u0 <- 3
  a <- frm_lincmt(parms = c(p, list(V = 1)), times = u0, ncmt = 1,
                  depot = TRUE, output = "central",
                  events = data.frame(time = 0, state = "depot",
                                      value = 100, ii = 8, ss = TRUE))
  # the ss trough at t = 0 plus the row's own dose at t = 0
  b <- ref_mp(1L, TRUE, p,
              list(list(amt = 100, u = u0, dur = 0, cin = 1L),
                   list(amt = 100, u = u0 + 8, dur = 0, cin = 1L,
                        ii = 8, n = 400L)), 2L)
  flush(stdout()); cat(sprintf("%12.1e %16.10e %16.10e %10.2e\n", d, a, as.numeric(b),
              relerr(a, b)))
}

cat("\n== two-compartment double root: k12 -> 0 with k21 == ke ==\n")
for (d in c(0, 10^-(16:0))) {
  p <- list(ke = 0.2, k12 = d, k21 = 0.2, ka = 1.1)
  u <- 6
  a <- frm_lincmt(parms = c(p, list(V = 1)), times = u, ncmt = 2,
                  depot = TRUE, init = list(depot = 100),
                  output = "central")
  b <- ref_mp(2L, TRUE, p, list(list(amt = 100, u = u, dur = 0,
                                     cin = 1L)), 2L)
  flush(stdout()); cat(sprintf("%12.1e %16.10e %16.10e %10.2e\n", d, a, as.numeric(b),
              relerr(a, b)))
}

cat("\n== three-compartment double root: k31 -> k21 ==\n")
for (d in c(0, 10^-(16:0))) {
  p <- list(ke = 0.2, k12 = 0.4, k21 = 0.1, k13 = 0.4, k31 = 0.1 + d,
            ka = 1.1)
  u <- 6
  a <- frm_lincmt(parms = c(p, list(V = 1)), times = u, ncmt = 3,
                  depot = TRUE, init = list(depot = 100),
                  output = "central")
  b <- ref_mp(3L, TRUE, p, list(list(amt = 100, u = u, dur = 0,
                                     cin = 1L)), 2L)
  flush(stdout()); cat(sprintf("%12.1e %16.10e %16.10e %10.2e\n", d, a, as.numeric(b),
              relerr(a, b)))
}

cat("\n== three-compartment triple root: every rate equal ==\n")
for (d in c(0, 10^-(12:1))) {
  q <- 0.3
  p <- list(ke = q, k12 = q, k21 = q + d, k13 = q, k31 = q - d,
            ka = q)
  for (u in c(2, 12)) {
    a <- frm_lincmt(parms = c(p, list(V = 1)), times = u, ncmt = 3,
                    depot = TRUE, init = list(depot = 100),
                    output = "central")
    b <- ref_mp(3L, TRUE, p, list(list(amt = 100, u = u, dur = 0,
                                       cin = 1L)), 2L)
    flush(stdout()); cat(sprintf("%12.1e %6g %16.10e %16.10e %10.2e\n", d, u, a,
                as.numeric(b), relerr(a, b)))
  }
}

cat("\n== infusion into the central compartment, 2 cmt ==\n")
for (dur in c(0.001, 0.1, 1, 4)) {
  for (u in c(dur / 2, dur, dur + 1e-9, 8)) {
    p <- list(ke = 0.2, k12 = 0.4, k21 = 0.1)
    a <- frm_lincmt(parms = c(p, list(V = 1)), times = u, ncmt = 2,
                    depot = FALSE, output = "central",
                    events = data.frame(time = 0, state = "central",
                                        value = 100, duration = dur))
    b <- ref_mp(2L, FALSE, p, list(list(amt = 100, u = u, dur = dur,
                                        cin = 1L)), 1L)
    flush(stdout()); cat(sprintf("dur %8.3f u %10.6f %16.10e %16.10e %10.2e\n", dur, u,
                a, as.numeric(b), relerr(a, b)))
  }
}

# ---------------------------------------------------------------------
# 2. Random sweep, worst case over the whole parameter box.

cat("\n== random sweep, worst relative error ==\n")
set.seed(11)
for (nc in 1:3) {
  worst <- 0; wat <- NULL
  for (i in 1:60) {
    p <- as.list(exp(runif(6, log(1e-3), log(5))))
    names(p) <- c("ke", "k12", "k21", "k13", "k31", "ka")
    u <- exp(runif(1, log(0.05), log(72)))
    a <- frm_lincmt(parms = c(p[c("ke", "k12", "k21", "k13",
                                  "k31", "ka")][
      c(TRUE, nc >= 2, nc >= 2, nc == 3, nc == 3, TRUE)],
      list(V = 1)),
      times = u, ncmt = nc, depot = TRUE,
      init = list(depot = 100), output = "central")
    b <- ref_mp(nc, TRUE, p, list(list(amt = 100, u = u, dur = 0,
                                       cin = 1L)), 2L)
    r <- relerr(a, b)
    if (is.finite(r) && r > worst) { worst <- r; wat <- c(unlist(p), u) }
  }
  cat(nc, "cmt + depot: worst rel", format(worst), "at",
      paste(format(wat, digits = 4), collapse = " "), "\n")
}
