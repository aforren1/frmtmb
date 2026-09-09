# Round 2, item 1. The lane's Boundary section argues the NaN gradient
# at a three-compartment double root is CORRECT, "because at a double
# root the eigenvalues have a square-root branch point, so the
# derivative with respect to any parameter that SPLITS the root is
# genuinely unbounded".
#
# That is true of the EIGENVALUES. It is not true of what
# frm_lincmt() RETURNS. The central compartment's impulse response is a
# symmetric function of the three roots, so it is a function of the
# characteristic polynomial's COEFFICIENTS, which are polynomials in
# the rate constants. The residue formula's poles at a double root are
# removable. An analytic function of an analytic function is analytic,
# so the derivative of the trajectory with respect to k31 exists and is
# finite AT the tangency, and the branch point cancels.
#
# This script tests that claim rather than the prose.
#
# Script path: dev/rev-lincmt-f6.R. Every point is constructed; the
# one seeded part names its seed.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")
suppressPackageStartupMessages(library(Rmpfr))

tt <- c(0.5, 1, 2, 4, 8, 12, 24, 48)
ev <- data.frame(time = 0, state = "depot", value = 100, ii = 8,
                 addl = 3L)
KE <- 0.2; K12 <- 0.4; K21 <- 0.1
b <- KE + K12 + K21
qlo <- (b - sqrt(b * b - 4 * KE * K21)) / 2

mk <- function(k13, k31) function(th) sum(frm_lincmt(
  parms = list(ke = exp(th[1]), k12 = K12, k21 = K21, k13 = k13,
               k31 = exp(th[2]), ka = exp(th[3]), V = 10),
  times = tt, ncmt = 3, depot = TRUE, events = ev))
fd1 <- function(f, x, j, h) {
  xp <- x; xp[[j]] <- xp[[j]] + h
  xm <- x; xm[[j]] <- xm[[j]] - h
  (f(xp) - f(xm)) / (2 * h)
}

cat("\n=== A. which component is NaN, and is the FD stable there? ===\n")
cat("k13 = 1e-300, k31 exactly on the slow root. The tape's gradient\n")
cat("against central differences at four step sizes.\n\n")
k13 <- 1e-300
th <- c(log(KE), log(qlo), log(1.1))
f <- mk(k13, qlo)
g <- as.numeric(MakeTape(f, th)$jacobian(th))
cat("  tape gradient   ",
    paste(format(g, digits = 8), collapse = "  "), "\n")
for (h in c(1e-4, 1e-5, 1e-6, 1e-7)) {
  gf <- vapply(1:3, function(j) fd1(f, th, j, h), 0)
  cat(sprintf("  central h=%.0e ", h),
      paste(format(gf, digits = 8), collapse = "  "), "\n")
}

cat("\n=== B. is the trajectory smooth THROUGH the tangency? ===\n")
cat("The tape's d/d(log k31) either side, and the central difference\n")
cat("AT it. A jump would mean the derivative really is unbounded.\n\n")
cat(sprintf("%14s %14s %14s\n", "k31/slow - 1", "eigen split",
            "d/d(log k31)"))
for (dk in c(-1e-4, -1e-6, -1e-8, 0, 1e-8, 1e-6, 1e-4)) {
  k31 <- qlo * (1 + dk)
  d <- frmtmb.ode:::lincmt_disp(3L, KE, K12, K21, k13, k31)
  lam <- sort(vapply(d[["lam"]], as.numeric, 0))
  ff <- mk(k13, k31)
  x <- c(log(KE), log(k31), log(1.1))
  gt <- as.numeric(MakeTape(ff, x)$jacobian(x))[[2L]]
  cat(sprintf("%14.0e %14.3e %14.8g\n", dk, min(diff(lam)), gt))
}
cat("\n  central difference of the same quantity AT the tangency:\n")
for (h in c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8))
  cat(sprintf("    h = %.0e   %.10g\n", h, fd1(f, th, 2L, h)))

cat("\n=== C. the true derivative, from a 300-bit reference ===\n")
cat("The reference never forms an eigenvalue, so it cannot inherit\n")
cat("the branch point. Its own central difference in log k31, taken\n")
cat("at 300 bits with a step of 1e-30, is the true derivative.\n\n")
PB <- 300L
mp <- function(x) mpfr(x, PB)
expm_ss <- function(A) {
  n <- nrow(A)
  nrm <- max(as.numeric(apply(abs(A), 1, sum)))
  s <- max(0L, as.integer(ceiling(log2(max(nrm, 1e-300)))) + 8L)
  B <- A / mp(2)^s
  I <- mpfrArray(0, PB, c(n, n))
  for (i in seq_len(n)) I[i, i] <- mp(1)
  T <- I; P <- I
  for (k in 1:60) { P <- P %*% B / mp(k); T <- T + P }
  for (k in seq_len(s)) T <- T %*% T
  T
}
# same schedule: four doses 8 apart into the depot, amount 100, V = 10
resp_mp <- function(ke, k31, ka, k13) {
  M <- mpfrArray(0, PB, c(4L, 4L))
  M[1, 1] <- -mp(ka); M[2, 1] <- mp(ka)
  M[2, 2] <- -(mp(ke) + mp(K12) + mp(k13))
  M[3, 2] <- mp(K12); M[2, 3] <- mp(K21); M[3, 3] <- -mp(K21)
  M[4, 2] <- mp(k13); M[2, 4] <- mp(k31); M[4, 4] <- -mp(k31)
  s <- mp(0)
  for (u in tt) {
    for (j in 0:3) {
      lag <- u - 8 * j
      if (lag > 0) s <- s + mp(100) * expm_ss(M * mp(lag))[2, 1]
    }
  }
  s / mp(10)
}
h <- mp(1e-30)
lk <- log(mp(qlo))
hi <- resp_mp(KE, exp(lk + h), 1.1, 1e-300)
lo <- resp_mp(KE, exp(lk - h), 1.1, 1e-300)
cat("  300-bit central difference d/d(log k31):",
    format((hi - lo) / (2 * h), digits = 20), "\n")
cat("  tape at the tangency:", format(g[[2L]], digits = 10), "\n")

cat("\n=== D. how far does the NaN region reach in k13? ===\n")
cat("k31 pinned exactly on the slow root; k13 raised until the tape\n")
cat("resolves the split.\n\n")
cat(sprintf("%10s %14s %10s\n", "k13", "eigen split", "finite"))
for (k in c(0, 10^-c(300, 40, 30, 22, 20, 18, 17, 16, 14, 12, 8))) {
  d <- frmtmb.ode:::lincmt_disp(3L, KE, K12, K21, k, qlo)
  lam <- sort(vapply(d[["lam"]], as.numeric, 0))
  ff <- mk(k, qlo)
  x <- c(log(KE), log(qlo), log(1.1))
  gg <- as.numeric(MakeTape(ff, x)$jacobian(x))
  cat(sprintf("%10.0e %14.3e %10s\n", k, min(diff(lam)),
              all(is.finite(gg))))
}
