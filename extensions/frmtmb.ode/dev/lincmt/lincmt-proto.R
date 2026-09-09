# lincmt lane: prototype of the analytic linear-compartment primitives.
# Plain doubles first; the AD check comes after.

TINY <- 1e-300

# (1 - exp(-x)) / x for x >= 0, with phi(0) = 1.
lin_phi <- function(x) {
  z <- x + TINY
  -expm1(-z) / z
}

# (exp(-p) - exp(-q)) / (q - p), symmetric, overflow-free for p, q >= 0.
lin_diff <- function(p, q) {
  d <- abs(p - q)
  mn <- (p + q - d) / 2
  exp(-mn) * lin_phi(d)
}

# E2(a, b; u) = (exp(-a u) - exp(-b u)) / (b - a) = conv of the two.
lin_e2 <- function(a, b, u) u * lin_diff(a * u, b * u)

# T(a, b; u1, ii) = sum_{j >= 0} E2(a, b; u1 + j ii)
lin_e2_ss <- function(a, b, u1, ii) {
  t1 <- u1 * lin_diff(a * u1, b * u1)
  P <- a * u1 + b * ii
  Q <- b * u1 + a * ii
  t2 <- (ii - u1) * lin_diff(P, Q)
  (t1 + t2) / ((-expm1(-a * ii)) * (-expm1(-b * ii)))
}

# --- disposition: eigenvalues and central-response coefficients -------

lin_disp <- function(ncmt, k10, k12 = 0, k21 = 0, k13 = 0, k31 = 0) {
  if (ncmt == 1L) return(list(lam = list(k10), coef = list(1)))
  if (ncmt == 2L) {
    g <- (k10 + k12 - k21) / 2
    D <- sqrt(g * g + k12 * k21 + TINY)
    m <- (k10 + k12 + k21) / 2
    lam1 <- m + D
    lam2 <- (k10 * k21) / lam1          # product of roots: stable
    c1 <- (D + g) / (2 * D)
    return(list(lam = list(lam1, lam2), coef = list(c1, 1 - c1)))
  }
  a2 <- k10 + k12 + k13 + k21 + k31
  a1 <- k10 * k21 + k10 * k31 + k21 * k31 + k12 * k31 + k13 * k21
  a0 <- k10 * k21 * k31
  p <- a1 - a2 * a2 / 3
  q <- -2 * a2^3 / 27 + a1 * a2 / 3 - a0
  r <- sqrt(pmax(-p / 3, 0) + TINY)
  arg <- (3 * q) / (2 * p) * sqrt(pmax(-3 / p, 0))
  arg <- arg / sqrt(1 + pmax(arg * arg - 1, 0))
  th <- acos(arg) / 3
  s <- a2 / 3
  l1 <- 2 * r * cos(th) + s
  l2 <- 2 * r * cos(th - 2 * pi / 3) + s
  l3 <- 2 * r * cos(th - 4 * pi / 3) + s
  den <- function(x, y, z) (x - y) * (x - z)
  c1 <- (l1 - k21) * (l1 - k31) / den(l1, l2, l3)
  c2 <- (l2 - k21) * (l2 - k31) / den(l2, l1, l3)
  c3 <- (l3 - k21) * (l3 - k31) / den(l3, l1, l2)
  list(lam = list(l1, l2, l3), coef = list(c1, c2, c3))
}

# --- checks ----------------------------------------------------------

cat("phi(0) =", lin_phi(0), " phi(1e-20) =", lin_phi(1e-20),
    " phi(1) =", lin_phi(1), " ref", (1 - exp(-1)), "\n")
cat("diff(2,2) =", lin_diff(2, 2), " exp(-2) =", exp(-2), "\n")
cat("diff(1000,1) =", lin_diff(1000, 1), " ref",
    (exp(-1000) - exp(-1)) / (1 - 1000), "\n")

# E2 against the naive form away from the singularity
set.seed(1)
a <- 1.3; b <- 0.21; u <- 3.7
cat("E2 stable", lin_e2(a, b, u), " naive",
    (exp(-a * u) - exp(-b * u)) / (b - a), "\n")

# E2 at coincidence
cat("E2(a,a,u) =", lin_e2(1.3, 1.3, 3.7), " ref", 3.7 * exp(-1.3 * 3.7),
    "\n")

# steady state sum by brute force
brute_ss <- function(a, b, u1, ii, n = 200000) {
  j <- 0:n
  sum((exp(-a * (u1 + j * ii)) - exp(-b * (u1 + j * ii))) / (b - a))
}
cat("ss E2: closed", lin_e2_ss(1.0, 0.15, 0.5, 12),
    " brute", brute_ss(1.0, 0.15, 0.5, 12), "\n")
cat("ss E2 u1>ii: closed", lin_e2_ss(1.0, 0.15, 25, 12),
    " brute", brute_ss(1.0, 0.15, 25, 12), "\n")
cat("ss E2 a==b: closed", lin_e2_ss(0.15, 0.15, 0.5, 12),
    " brute", brute_ss(0.15, 0.15 + 1e-9, 0.5, 12), "\n")

# disposition: 2 cmt against eigen()
K2 <- function(k10, k12, k21) {
  matrix(c(-(k10 + k12), k21, k12, -k21), 2, 2)
}
chk2 <- function(k10, k12, k21, u) {
  d <- lin_disp(2L, k10, k12, k21)
  cf <- sum(mapply(function(cc, ll) cc * exp(-ll * u), d$coef, d$lam))
  ref <- (as.matrix(Matrix::expm(K2(k10, k12, k21) * u)))[1, 1]
  c(closed = cf, ref = ref, rel = abs(cf - ref) / abs(ref))
}
if (requireNamespace("Matrix", quietly = TRUE)) {
  print(chk2(0.15, 0.4, 0.2, 3))
  print(chk2(0.15, 1e-12, 0.15, 3))   # degenerate: alpha == beta
  print(chk2(0.15, 1e-6, 0.150001, 7))
}

K3 <- function(k10, k12, k21, k13, k31) {
  matrix(c(-(k10 + k12 + k13), k12, k13,
           k21, -k21, 0,
           k31, 0, -k31), 3, 3)
}
chk3 <- function(k10, k12, k21, k13, k31, u) {
  d <- lin_disp(3L, k10, k12, k21, k13, k31)
  cf <- sum(mapply(function(cc, ll) cc * exp(-ll * u), d$coef, d$lam))
  ref <- (as.matrix(Matrix::expm(K3(k10, k12, k21, k13, k31) * u)))[1, 1]
  c(closed = cf, ref = ref, rel = abs(cf - ref) / abs(ref))
}
if (requireNamespace("Matrix", quietly = TRUE)) {
  print(chk3(0.15, 0.4, 0.2, 0.05, 0.01, 3))
  print(chk3(0.15, 0.4, 0.2, 0.05, 0.01, 40))
  cat("coefs 3cmt:", unlist(lin_disp(3L, .15, .4, .2, .05, .01)$coef), "\n")
  cat("lams  3cmt:", unlist(lin_disp(3L, .15, .4, .2, .05, .01)$lam), "\n")
}
