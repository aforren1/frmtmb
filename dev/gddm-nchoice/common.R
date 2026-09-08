# Shared setup for the n-choice gddm feasibility probes.
#
# Prototype code. Nothing here is loaded by any package; the probes reach
# into frmtmb.eam with ::: where they want the shipped 1-D solver as a
# baseline.

suppressPackageStartupMessages(library(RTMB))

gd_time <- function(expr, reps = 1L) {
  e <- substitute(expr); p <- parent.frame()
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(reps)) v <- eval(e, p)
  list(sec = (proc.time()[["elapsed"]] - t0) / reps, value = v)
}

# --- the n-alternative geometry -------------------------------------------
#
# n accumulators x_1..x_n, dx_i = mu_i dt + dW_i, unit independent noise.
# Only relative evidence decides, so project onto the sum-zero plane. With U
# an orthonormal basis of that plane (n by n-1), z = U'x has Cov(dz) = U'U =
# I: the reduced process is ISOTROPIC Brownian motion in n-1 dimensions with
# drift U'mu. This is Roxin's reduction, and the isotropy is what makes the
# reflection-group machinery in the image probe apply at all.
#
# Stopping rule: alternative k wins when its evidence beats the mean of all n
# by c, i.e. <v_k, z> >= c with v_k = U'(e_k - 1/n). The n vectors v_k have
# equal length sqrt((n-1)/n) and equal pairwise angles, so the continuation
# region is a REGULAR (n-1)-simplex centred on the origin.
#
# At n = 2 the rule is x_1 - x_2 = +/- 2c, so z = (x_1-x_2)/sqrt(2) is a
# unit-diffusion Wiener process with drift (mu_1-mu_2)/sqrt(2) absorbed at
# +/- c*sqrt(2). In wiener() parameterization: bs = 2*sqrt(2)*c,
# v = (mu_1-mu_2)/sqrt(2), bias = 0.5. That is the closed-form check.

# Helmert basis: deterministic, so every probe sees the same coordinates.
gd_basis <- function(n) {
  U <- matrix(0, n, n - 1L)
  for (k in seq_len(n - 1L)) {
    U[seq_len(k), k] <- 1
    U[k + 1L, k] <- -k
    U[, k] <- U[, k] / sqrt(k * (k + 1))
  }
  U
}

gd_normals <- function(n) {
  U <- gd_basis(n)
  list(U = U, V = t(U) %*% (diag(n) - matrix(1 / n, n, n)))
}

# wiener() parameters equivalent to the n = 2 case of the above.
gd_wiener_equiv <- function(mu, cc) {
  list(v = (mu[1] - mu[2]) / sqrt(2), bs = 2 * sqrt(2) * cc, bias = 0.5)
}

# --- direct simulation ----------------------------------------------------
#
# Euler-Maruyama on the reduced process with a Brownian-bridge absorption
# correction, which shares no code with any of the densities.
#
# WHY THE BRIDGE. Plain Euler only notices a crossing when a monitoring point
# lands outside, so it misses excursions between points and reports first
# passages LATE, with a bias of order sqrt(dt) in the boundary distance. At a
# step fine enough for that bias to sit below the Monte Carlo error, 200k
# trials cost tens of minutes. Conditional on the two endpoints the normal
# component of the increment is a Brownian bridge, and the probability that
# it touched a flat wall is exp(-2 d0 d1 / dt) in closed form, so absorbing
# with that probability removes the leading bias and lets the step be 25
# times coarser. The three walls are treated as independent within a step,
# which is exact except within a corner.
gd_sim_nchoice <- function(ntrial, mu, cc, dt = 5e-4, t_max = 8,
                           z0 = NULL, seed = NULL, bridge = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(mu)
  nv <- gd_normals(n)
  a <- as.numeric(t(nv$U) %*% mu)
  d <- n - 1L
  vn <- sqrt(sum(nv$V[, 1L]^2))          # |v_k|, common to all k
  z <- matrix(if (is.null(z0)) 0 else z0, nrow = ntrial, ncol = d, byrow = TRUE)
  alive <- rep(TRUE, ntrial)
  tt <- rep(NA_real_, ntrial)
  ch <- rep(NA_integer_, ntrial)
  sdt <- sqrt(dt)
  d0 <- (cc - z %*% nv$V) / vn            # normal distances to each wall
  for (k in seq_len(as.integer(round(t_max / dt)))) {
    m <- sum(alive)
    if (!m) break
    ia <- which(alive)
    zz <- z[ia, , drop = FALSE] +
      rep(a * dt, each = m) + matrix(stats::rnorm(m * d, 0, sdt), m, d)
    d1 <- (cc - zz %*% nv$V) / vn
    da <- d0[ia, , drop = FALSE]
    if (bridge) {
      pc <- exp(-2 * pmax(da, 0) * pmax(d1, 0) / dt)
      pc[d1 <= 0] <- 1
      fire <- matrix(stats::runif(m * n), m, n) < pc
    } else {
      fire <- d1 <= 0
    }
    z[ia, ] <- zz
    d0[ia, ] <- d1
    hit <- rowSums(fire) > 0
    if (any(hit)) {
      idx <- ia[hit]
      # the step is the resolution, so date the crossing at its midpoint
      ch[idx] <- max.col(fire[hit, , drop = FALSE] * (1 - d1[hit, , drop = FALSE]),
                         ties.method = "first")
      tt[idx] <- (k - 0.5) * dt
      alive[idx] <- FALSE
    }
  }
  data.frame(choice = ch, rt = tt)
}
