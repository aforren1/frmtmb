# Degeneracy sweep: how the closed forms behave as rate constants
# converge. Reference is Matrix::expm on the same rate matrix, which is
# a scaling-and-squaring Pade and does not use the eigenvalues.
source("lincmt-proto.R")

K2 <- function(k10, k12, k21) matrix(c(-(k10 + k12), k21, k12, -k21), 2, 2)
K3 <- function(k10, k12, k21, k13, k31)
  matrix(c(-(k10 + k12 + k13), k12, k13, k21, -k21, 0, k31, 0, -k31), 3, 3)

h2 <- function(k10, k12, k21, u) {
  d <- lin_disp(2L, k10, k12, k21)
  sum(mapply(function(cc, ll) cc * exp(-ll * u), d$coef, d$lam))
}
h3 <- function(k10, k12, k21, k13, k31, u) {
  d <- lin_disp(3L, k10, k12, k21, k13, k31)
  sum(mapply(function(cc, ll) cc * exp(-ll * u), d$coef, d$lam))
}

cat("\n== 2 cmt: k12 -> 0 with k21 == k10 (the double root) ==\n")
cat(sprintf("%10s %14s %14s %10s %12s\n", "k12", "alpha-beta",
            "closed", "expm", "rel"))
for (e in c(0, 10^-(16:1))) {
  k10 <- 0.15; k21 <- 0.15; k12 <- e
  d <- lin_disp(2L, k10, k12, k21)
  gap <- d$lam[[1]] - d$lam[[2]]
  cf <- h2(k10, k12, k21, 3)
  rf <- as.matrix(Matrix::expm(K2(k10, k12, k21) * 3))[1, 1]
  cat(sprintf("%10.1e %14.6e %14.10f %10.6f %12.2e\n", e, gap, cf, rf,
              abs(cf - rf) / abs(rf)))
}

cat("\n== 2 cmt: k21 -> k10 with k12 = 1e-6 ==\n")
for (e in c(0, 10^-(12:1))) {
  k10 <- 0.15; k21 <- 0.15 + e; k12 <- 1e-6
  cf <- h2(k10, k12, k21, 3)
  rf <- as.matrix(Matrix::expm(K2(k10, k12, k21) * 3))[1, 1]
  cat(sprintf("%10.1e %14.10f %14.10f %12.2e\n", e, cf, rf,
              abs(cf - rf) / abs(rf)))
}

cat("\n== 3 cmt: k31 -> k21 (peripheral rates converge) ==\n")
cat(sprintf("%10s %14s %14s %14s %12s\n", "eps", "min gap", "closed",
            "expm", "rel"))
for (e in c(0, 10^-(14:1))) {
  k10 <- 0.15; k12 <- 0.4; k21 <- 0.2; k13 <- 0.4; k31 <- 0.2 + e
  d <- lin_disp(3L, k10, k12, k21, k13, k31)
  ll <- sort(unlist(d$lam))
  gap <- min(diff(ll))
  cf <- h3(k10, k12, k21, k13, k31, 3)
  rf <- as.matrix(Matrix::expm(K3(k10, k12, k21, k13, k31) * 3))[1, 1]
  cat(sprintf("%10.1e %14.4e %14.10f %14.10f %12.2e\n", e, gap, cf, rf,
              abs(cf - rf) / abs(rf)))
}

cat("\n== 3 cmt: k13 -> 0 (third compartment vanishes) ==\n")
for (e in c(0, 10^-(14:1))) {
  k10 <- 0.15; k12 <- 0.4; k21 <- 0.2; k13 <- e; k31 <- 0.05
  d <- lin_disp(3L, k10, k12, k21, k13, k31)
  ll <- sort(unlist(d$lam))
  cf <- h3(k10, k12, k21, k13, k31, 3)
  rf <- as.matrix(Matrix::expm(K3(k10, k12, k21, k13, k31) * 3))[1, 1]
  cat(sprintf("%10.1e %14.4e %14.10f %14.10f %12.2e\n", e, min(diff(ll)),
              cf, rf, abs(cf - rf) / abs(rf)))
}

cat("\n== 3 cmt: random sweep, worst relative error ==\n")
set.seed(4)
worst <- 0; where <- NULL
for (i in 1:20000) {
  p <- exp(runif(5, log(1e-4), log(10)))
  u <- exp(runif(1, log(0.05), log(200)))
  cf <- h3(p[1], p[2], p[3], p[4], p[5], u)
  rf <- as.matrix(Matrix::expm(K3(p[1], p[2], p[3], p[4], p[5]) * u))[1, 1]
  if (rf > 1e-12) {
    r <- abs(cf - rf) / abs(rf)
    if (is.finite(r) && r > worst) { worst <- r; where <- c(p, u) }
  }
}
cat("worst rel:", format(worst), "at", format(where), "\n")

cat("\n== 2 cmt: random sweep ==\n")
set.seed(5)
worst <- 0; where <- NULL
for (i in 1:20000) {
  p <- exp(runif(3, log(1e-4), log(10)))
  u <- exp(runif(1, log(0.05), log(200)))
  cf <- h2(p[1], p[2], p[3], u)
  rf <- as.matrix(Matrix::expm(K2(p[1], p[2], p[3]) * u))[1, 1]
  if (rf > 1e-12) {
    r <- abs(cf - rf) / abs(rf)
    if (is.finite(r) && r > worst) { worst <- r; where <- c(p, u) }
  }
}
cat("worst rel:", format(worst), "at", format(where), "\n")
