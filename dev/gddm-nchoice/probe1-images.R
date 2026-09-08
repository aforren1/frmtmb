# Probe 1: is the three-alternative image series right?
#
# Run from the worktree root:
#   Rscript dev/gddm-nchoice/probe1-images.R
#
# Checks, four of the five against something that does not go through the
# series at all:
#
#  A. the image sum vanishes on all three walls, and its integral over the
#     triangle is a survival function falling from one;
#  B. the three defective densities integrate to one, for any drift;
#  C. zero drift from the centre gives exactly 1/3 per choice;
#  D. zero drift gives mean first-passage time (2/3r) prod_k (r - <n_k,z0>),
#     the torsion function of the equilateral triangle. That is exact:
#     l_k = r - <n_k,z> are affine with unit gradients meeting at 120
#     degrees and summing to 3r, so Laplacian(l1 l2 l3) = -3r exactly, and
#     u = 2 l1 l2 l3 / 3r solves (1/2)Lap u = -1 with u = 0 on the walls;
#  E. against Brownian-bridge simulation, choice proportions and response
#     time quantiles, with Monte Carlo error.

source("dev/gddm-nchoice/common.R")
source("dev/gddm-nchoice/images.R")

tl <- gd_tiles(400L)
gg <- tl$geom
cat("tiles:", length(tl$sgn), " first radii:",
    paste(sprintf("%.2f", sqrt(rowSums(tl$w[1:10, ]^2))), collapse = " "), "\n")

trim <- function(tl, m) list(A = tl$A[seq_len(m)],
                             w = tl$w[seq_len(m), , drop = FALSE],
                             sgn = tl$sgn[seq_len(m)], geom = tl$geom)

# quadrature on a fixed fine grid: the density is smooth and vanishes at
# both ends, so the trapezoid rule is enough and, unlike integrate(), it
# cannot silently mis-sample a spike.
tgrid <- c(seq(1e-4, 4, by = 5e-4), seq(4.01, 40, by = 0.01))
trap <- function(v, t) sum(diff(t) * (v[-1L] + v[-length(v)]) / 2)

summ <- function(r, a, z0, tl) {
  f <- lapply(1:3, function(k) gd_img_density(tgrid, a, z0, r, tl, k))
  P <- vapply(f, trap, numeric(1), t = tgrid)
  M <- vapply(f, function(v) trap(v * tgrid, tgrid), numeric(1))
  list(P = P, mass = sum(P), mean = sum(M) / sum(P), f = f)
}

tor <- function(r, z0)
  2 / (3 * r) * prod(vapply(1:3, function(k) r - sum(gg$n[, k] * z0), numeric(1)))

r1 <- sqrt(1.5) * 0.7071   # c = 0.7071, the scale of a two-choice bs = 2

# --- A: the kernel itself -------------------------------------------------
qq <- function(z, t, z0, r, tl) {
  s <- 0
  for (i in seq_along(tl$sgn)) {
    y <- as.numeric(tl$A[[i]] %*% z0) + r * tl$w[i, ]
    s <- s + tl$sgn[i] * exp(-sum((z - y)^2) / (2 * t)) / (2 * pi * t)
  }
  s
}
z0t <- c(0.1, -0.05)
wv <- max(abs(vapply(1:3, function(k)
  max(abs(vapply(c(-1, 0, 1), function(u)
    qq(r1 * gg$n[, k] + u * gg$m[, k], 0.3, z0t, r1, tl), numeric(1)))),
  numeric(1))))
cat(sprintf("A. worst |q| on the three walls at t = 0.3: %.2e\n", wv))

h <- 0.005
gx <- seq(-2 * r1, 2 * r1, by = h)
G <- as.matrix(expand.grid(x = gx, y = gx))
ins <- rep(TRUE, nrow(G))
for (k in 1:3) ins <- ins & (G[, 1] * gg$n[1, k] + G[, 2] * gg$n[2, k] < r1)
Gi <- G[ins, , drop = FALSE]
surv <- function(t, z0, r, tl) {
  s <- 0
  for (i in seq_along(tl$sgn)) {
    y <- as.numeric(tl$A[[i]] %*% z0) + r * tl$w[i, ]
    s <- s + tl$sgn[i] *
      sum(exp(-((Gi[, 1] - y[1])^2 + (Gi[, 2] - y[2])^2) / (2 * t))) / (2 * pi * t)
  }
  s * h * h
}
cat(sprintf("%8s %12s %14s %14s %10s\n",
            "t", "S(t) grid", "-dS/dt", "sum_k f_k(t)", "rel"))
for (t in c(0.1, 0.3, 0.6, 1.0, 1.5)) {
  dd <- 0.01
  ds <- -(surv(t + dd, z0t, r1, tl) - surv(t - dd, z0t, r1, tl)) / (2 * dd)
  fs <- sum(vapply(1:3, function(k) gd_img_density(t, c(0, 0), z0t, r1, tl, k),
                   numeric(1)))
  cat(sprintf("%8.2f %12.6f %14.6f %14.6f %10.2e\n",
              t, surv(t, z0t, r1, tl), ds, fs, abs(ds - fs) / fs))
}

# --- how many images, for a horizon --------------------------------------
cat("\n-- images needed, drift (1.2, -0.4), r =", sprintf("%.4f", r1), "--\n")
ref <- vapply(1:3, function(k)
  gd_img_density(seq(0.05, 3, by = 0.05), c(1.2, -0.4), c(0, 0), r1, tl, k),
  numeric(60))
cat(sprintf("%8s %14s\n", "images", "max rel err"))
for (m in c(4L, 10L, 22L, 46L, 76L, 118L, 166L, 250L)) {
  v <- vapply(1:3, function(k)
    gd_img_density(seq(0.05, 3, by = 0.05), c(1.2, -0.4), c(0, 0), r1,
                   trim(tl, m), k), numeric(60))
  cat(sprintf("%8d %14.3e\n", m, max(abs(v - ref) / pmax(abs(ref), 1e-14))))
}

# --- B/C/D ----------------------------------------------------------------
cat("\n-- B/C/D: mass, symmetry, torsion-function mean --\n")
cat(sprintf("%-20s %12s %10s %10s %10s %11s %11s\n",
            "case", "mass", "P1", "P2", "P3", "mean", "exact mean"))
cases <- list(
  list(lab = "zero drift, centre", a = c(0, 0),      z0 = c(0, 0)),
  list(lab = "zero drift, offset", a = c(0, 0),      z0 = c(0.25, -0.15)),
  list(lab = "drift (1.2,-0.4)",   a = c(1.2, -0.4), z0 = c(0, 0)),
  list(lab = "drift (2.5, 1.0)",   a = c(2.5, 1.0),  z0 = c(0.2, 0.1)),
  list(lab = "drift (0,0), r=1.6", a = c(0, 0),      z0 = c(0, 0), r = 1.6))
for (cs in cases) {
  rr <- if (is.null(cs$r)) r1 else cs$r
  s <- summ(rr, cs$a, cs$z0, tl)
  ex <- if (all(cs$a == 0)) sprintf("%.7f", tor(rr, cs$z0)) else "-"
  cat(sprintf("%-20s %12.9f %10.6f %10.6f %10.6f %11.7f %11s\n",
              cs$lab, s$mass, s$P[1], s$P[2], s$P[3], s$mean, ex))
}

# --- E: simulation --------------------------------------------------------
cat("\n-- E: Brownian-bridge simulation --\n")
mu <- c(0.9, 0.1, -1.0); cc <- 0.7071
a <- as.numeric(t(gd_basis(3L)) %*% mu)
r <- sqrt(1.5) * cc
cat(sprintf("mu = (%s), c = %.4f -> a = (%.4f, %.4f), r = %.4f\n",
            paste(mu, collapse = ", "), cc, a[1], a[2], r))
for (nsim in c(100000L)) for (sdt in c(1e-3, 2.5e-4)) {
  t0 <- proc.time()[["elapsed"]]
  sim <- gd_sim_nchoice(nsim, mu, cc, dt = sdt, t_max = 8, seed = 7)
  el <- proc.time()[["elapsed"]] - t0
  s <- summ(r, a, c(0, 0), tl)
  n <- sum(!is.na(sim$choice))
  cat(sprintf("\nn = %d, step %.1e, %.1f s, %d unresolved\n",
              nsim, sdt, el, sum(is.na(sim$rt))))
  cat(sprintf("%8s %11s %11s %10s %7s\n", "choice", "series P", "sim P", "se", "z"))
  for (k in 1:3) {
    ph <- mean(sim$choice == k, na.rm = TRUE)
    se <- sqrt(ph * (1 - ph) / n)
    cat(sprintf("%8d %11.6f %11.6f %10.6f %7.2f\n", k, s$P[k], ph, se,
                (ph - s$P[k]) / se))
  }
  pq <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  cat(sprintf("%8s %5s %10s %10s %9s %7s\n",
              "choice", "q", "series", "sim", "se", "z"))
  for (k in 1:3) {
    fk <- s$f[[k]]
    cdf <- c(0, cumsum(diff(tgrid) * (fk[-1L] + fk[-length(fk)]) / 2)) / s$P[k]
    qs <- stats::approx(cdf, tgrid, xout = pq)$y
    rts <- sim$rt[!is.na(sim$choice) & sim$choice == k]
    qh <- stats::quantile(rts, pq, names = FALSE)
    dens <- stats::approx(tgrid, fk / s$P[k], xout = qs)$y
    se <- sqrt(pq * (1 - pq) / length(rts)) / dens
    for (i in seq_along(pq))
      cat(sprintf("%8d %5.2f %10.5f %10.5f %9.5f %7.2f\n",
                  k, pq[i], qs[i], qh[i], se[i], (qh[i] - qs[i]) / se[i]))
  }
}
