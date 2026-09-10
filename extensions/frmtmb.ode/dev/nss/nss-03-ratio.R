# M1/M3. Does the contraction ratio estimate exp(-lambda_z * ii), and
# what does a geometric-tail extrapolation of the run-in buy?
#
# The cycle map of a LINEAR compartment system is affine,
# y -> M (y + a) with M = expm(A * ii), so every quantity here has a
# closed form and the arithmetic is exact to machine precision. That
# separates two error sources the ODE path mixes: truncation of the
# run-in, and the integrator's own tolerance.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-03-ratio.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages(library(Matrix))

# amounts, two compartments plus a depot, matching frm_lincmt()'s
# micro-constant parameterization
Amat <- function(ke, k12, k21, ka) {
  matrix(c(-ka, 0, 0,
           ka, -(ke + k12), k21,
           0, k12, -k21), 3L, 3L, byrow = TRUE)
}
lamz <- function(ke, k12, k21) {
  b <- ke + k12 + k21
  (b - sqrt(b * b - 4 * ke * k21)) / 2
}
lam1 <- function(ke, k12, k21) {
  b <- ke + k12 + k21
  (b + sqrt(b * b - 4 * ke * k21)) / 2
}

# the run-in exactly as ode_run_in() performs it: start at zero, add
# the dose to the depot, advance one interval, repeat
cycles <- function(ke, k12, k21, ka, ii, amt, n) {
  M <- as.matrix(expm(Amat(ke, k12, k21, ka) * ii))
  a <- c(amt, 0, 0)
  y <- c(0, 0, 0)
  ys <- matrix(0, n + 1L, 3L)
  ys[1L, ] <- y
  for (k in seq_len(n)) {
    y <- as.numeric(M %*% (y + a))
    ys[k + 1L, ] <- y
  }
  # the exact fixed point of y -> M (y + a)
  ystar <- as.numeric(solve(diag(3) - M, M %*% a))
  list(ys = ys, ystar = ystar, M = M)
}

# two estimators of the contraction ratio from three iterates
r_norm <- function(d1, d2) sqrt(sum(d2^2)) / sqrt(sum(d1^2))
r_dot  <- function(d1, d2) sum(d1 * d2) / sum(d1 * d1)

grid <- list(
  c(ke = 0.2,  k12 = 0.4,  k21 = 0.1,   ka = 1.1, ii = 8),
  c(ke = 0.2,  k12 = 0.4,  k21 = 0.1,   ka = 1.1, ii = 24),
  c(ke = 0.5,  k12 = 1.0,  k21 = 0.05,  ka = 1.1, ii = 12),
  c(ke = 0.15, k12 = 0.3,  k21 = 0.02,  ka = 1.0, ii = 24),
  c(ke = 0.15, k12 = 0.3,  k21 = 0.02,  ka = 1.0, ii = 12),
  c(ke = 0.1,  k12 = 0.2,  k21 = 0.008, ka = 1.0, ii = 24),
  c(ke = 0.08, k12 = 0.15, k21 = 0.003, ka = 1.0, ii = 24),
  c(ke = 0.05, k12 = 0.1,  k21 = 0.001, ka = 1.0, ii = 168))

cat("=== A. the ratio estimators against exp(-lambda_z * ii) ===\n")
cat("r_true is the dominant eigenvalue of M = expm(A ii), which is\n",
    "exp(-lambda_z ii) unless ka is slower. r_norm and r_dot are read\n",
    "off the last three iterates of an n-cycle run-in.\n\n", sep = "")
cat(sprintf("%9s %5s %9s %9s %9s %9s %9s %9s\n",
            "t_half_z", "ii", "r_true", "n=3 norm", "n=3 dot",
            "n=5 norm", "n=5 dot", "n=20 dot"))
for (g in grid) {
  z <- cycles(g[["ke"]], g[["k12"]], g[["k21"]], g[["ka"]], g[["ii"]],
              100, 22L)
  rt <- max(Mod(eigen(z$M, only.values = TRUE)$values))
  d <- diff(z$ys)                      # d[k, ] is y_k - y_{k-1}
  get <- function(n, f) f(d[n - 1L, ], d[n, ])
  cat(sprintf("%9.1f %5g %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n",
              log(2) / lamz(g[["ke"]], g[["k12"]], g[["k21"]]),
              g[["ii"]], rt,
              get(3L, r_norm), get(3L, r_dot),
              get(5L, r_norm), get(5L, r_dot), get(20L, r_dot)))
}

cat("\n=== B. n_ss the 'auto' rule would pick, at tol 1e-9 ===\n")
cat(sprintf("%9s %5s %9s %8s %8s %8s %8s\n", "t_half_z", "ii",
            "r_true", "exact", "n=3 dot", "n=5 dot", "n=8 dot"))
pick <- function(r) if (r <= 0 || r >= 1) NA_integer_ else
  as.integer(ceiling(log(1e-9) / log(r)))
for (g in grid) {
  z <- cycles(g[["ke"]], g[["k12"]], g[["k21"]], g[["ka"]], g[["ii"]],
              100, 10L)
  rt <- max(Mod(eigen(z$M, only.values = TRUE)$values))
  d <- diff(z$ys)
  cat(sprintf("%9.1f %5g %9.6f %8d %8d %8d %8d\n",
              log(2) / lamz(g[["ke"]], g[["k12"]], g[["k21"]]),
              g[["ii"]], rt, pick(rt),
              pick(r_dot(d[2L, ], d[3L, ])),
              pick(r_dot(d[4L, ], d[5L, ])),
              pick(r_dot(d[7L, ], d[8L, ]))))
}

cat("\n=== C. what a geometric-tail extrapolation buys at n_ss = 20 ===\n")
cat("shortfall is max|y_n - y*| / max|y*| over the three states.\n",
    "'extrap' adds d_n * r/(1 - r) with r from the last two",
    " differences.\n\n", sep = "")
cat(sprintf("%9s %5s %11s %11s %11s %11s\n", "t_half_z", "ii",
            "plain 20", "extrap 20", "plain 5", "extrap 5"))
short <- function(y, ys) max(abs(y - ys)) / max(abs(ys))
ext <- function(ys, n) {
  d <- diff(ys)
  r <- r_dot(d[n - 1L, ], d[n, ])
  ys[n + 1L, ] + d[n, ] * r / (1 - r)
}
for (g in grid) {
  z <- cycles(g[["ke"]], g[["k12"]], g[["k21"]], g[["ka"]], g[["ii"]],
              100, 20L)
  cat(sprintf("%9.1f %5g %11.3e %11.3e %11.3e %11.3e\n",
              log(2) / lamz(g[["ke"]], g[["k12"]], g[["k21"]]),
              g[["ii"]],
              short(z$ys[21L, ], z$ystar),
              short(ext(z$ys, 20L), z$ystar),
              short(z$ys[6L, ], z$ystar),
              short(ext(z$ys, 5L), z$ystar)))
}

cat("\n=== D. the same, through deSolve at the package's tolerances ===\n")
cat("atol = rtol = 1e-8, which is frm_ode()'s default, so this is the\n",
    "noise floor a real run-in carries.\n\n", sep = "")
dyn <- function(t, y, p) {
  list(c(-p[["ka"]] * y[1L],
         p[["ka"]] * y[1L] - (p[["ke"]] + p[["k12"]]) * y[2L] +
           p[["k21"]] * y[3L],
         p[["k12"]] * y[2L] - p[["k21"]] * y[3L]))
}
cycles_num <- function(g, amt, n, atol = 1e-8, rtol = 1e-8) {
  p <- as.list(g)
  y <- c(0, 0, 0)
  ys <- matrix(0, n + 1L, 3L)
  for (k in seq_len(n)) {
    y[1L] <- y[1L] + amt
    s <- deSolve::ode(y, c(0, g[["ii"]]), dyn, p, method = "lsoda",
                      atol = atol, rtol = rtol)
    y <- as.numeric(s[2L, 2:4])
    ys[k + 1L, ] <- y
  }
  ys
}
cat(sprintf("%9s %5s %11s %11s %11s %11s\n", "t_half_z", "ii",
            "plain 20", "extrap 20", "plain 5", "extrap 5"))
for (g in grid) {
  z <- cycles(g[["ke"]], g[["k12"]], g[["k21"]], g[["ka"]], g[["ii"]],
              100, 20L)
  yn <- cycles_num(g, 100, 20L)
  cat(sprintf("%9.1f %5g %11.3e %11.3e %11.3e %11.3e\n",
              log(2) / lamz(g[["ke"]], g[["k12"]], g[["k21"]]),
              g[["ii"]],
              short(yn[21L, ], z$ystar),
              short(ext(yn, 20L), z$ystar),
              short(yn[6L, ], z$ystar),
              short(ext(yn, 5L), z$ystar)))
}
