# Shared GDDM machinery for the eLife 56938 feasibility probes.
# Sourced by probe-04 (tape cost and first fit) and probe-05 (solver accuracy and a
# recovery study against the solver's own density). See probe-04's header for why each
# piece is written the way it is.

CMAX <- 0.512
COH <- c(0, 0.032, 0.064, 0.128, 0.256, 0.512)
T_DUR <- 2
NDT_MAX <- 0.4
LAPSE <- 0.05

# ---- the solver, written so RTMB can tape it ---------------------------------------
# Returns defective densities at both walls on the fixed time grid, for one condition.
gddm_solve <- function(mu, leak, B0, tau, nt, ny, dt) {
  "[<-" <- ADoverload("[<-")
  h <- 2 / (ny + 1)
  y <- seq(-1 + h, 1 - h, length.out = ny)

  rho <- numeric(ny) * B0 * 0            # seeded from a parameter so the tape sees AD
  rho[(ny + 1L) %/% 2L] <- 1 / h
  p_up <- numeric(nt + 1L) * B0 * 0
  p_lo <- numeric(nt + 1L) * B0 * 0

  Bt <- function(t) B0 * exp(-t / tau)

  for (k in seq_len(nt)) {
    B0t <- Bt((k - 1) * dt); B1t <- Bt(k * dt)
    a0 <- mu / B0t - leak * y + y / tau;  D0 <- 0.5 / (B0t * B0t)
    a1 <- mu / B1t - leak * y + y / tau;  D1 <- 0.5 / (B1t * B1t)

    # explicit half-step, M(t0) %*% rho, tridiagonal
    lo0 <- a0[seq_len(ny - 1L)] / (2 * h) + D0 / h^2
    up0 <- -a0[2:ny] / (2 * h) + D0 / h^2
    mv <- (-2 * D0 / h^2) * rho
    mv[seq_len(ny - 1L)] <- mv[seq_len(ny - 1L)] + up0 * rho[2:ny]
    mv[2:ny] <- mv[2:ny] + lo0 * rho[seq_len(ny - 1L)]
    rhs <- rho + (dt / 2) * mv

    # implicit half-step: Thomas sweep on (I - dt/2 M(t1))
    lo <- -(dt / 2) * (a1[seq_len(ny - 1L)] / (2 * h) + D1 / h^2)
    up <- -(dt / 2) * (-a1[2:ny] / (2 * h) + D1 / h^2)
    di <- (1 + dt * D1 / h^2) * (numeric(ny) + 1)
    cp <- numeric(ny) * B0 * 0; dp <- numeric(ny) * B0 * 0
    cp[1L] <- up[1L] / di[1L]; dp[1L] <- rhs[1L] / di[1L]
    for (i in 2:ny) {
      m <- di[i] - lo[i - 1L] * cp[i - 1L]
      cp[i] <- if (i < ny) up[i] / m else 0
      dp[i] <- (rhs[i] - lo[i - 1L] * dp[i - 1L]) / m
    }
    new <- numeric(ny) * B0 * 0
    new[ny] <- dp[ny]
    for (i in (ny - 1L):1L) new[i] <- dp[i] - cp[i] * new[i + 1L]

    # second-order one-sided flux against the zero-density wall
    p_up[k + 1L] <- 0.5 * (D0 * (4 * rho[ny] - rho[ny - 1L]) / (2 * h) +
                           D1 * (4 * new[ny] - new[ny - 1L]) / (2 * h))
    p_lo[k + 1L] <- 0.5 * (D0 * (4 * rho[1L] - rho[2L]) / (2 * h) +
                           D1 * (4 * new[1L] - new[2L]) / (2 * h))
    rho <- new
  }
  list(up = p_up, lo = p_lo)
}

# Interpolation kernel. The obvious choice, a Catmull-Rom cubic, is defined piecewise on
# |u| and so branches on t_nd; RTMB rejects it outright ("Comparison is generally unsafe
# for AD types"). The cubic B-spline has an equivalent truncated-power form that needs no
# branch at all, because max(x, 0) = (x + |x|)/2 and abs() is overloaded:
#
#   B3(u) = [ r(u+2)^3 - 4 r(u+1)^3 + 6 r(u)^3 - 4 r(u-1)^3 + r(u-2)^3 ] / 6
#
# It is C2 in t_nd and a partition of unity, so the shift conserves mass exactly. It
# smooths rather than interpolates, but the smoothing error is O(dt^2 p'') and dt^2 is
# 1e-4 here against a density that varies on a 0.1 s scale.
relu <- function(x) 0.5 * (x + abs(x))
b3_kernel <- function(u)
  (relu(u + 2)^3 - 4 * relu(u + 1)^3 + 6 * relu(u)^3 -
     4 * relu(u - 1)^3 + relu(u - 2)^3) / 6

# Shift a density on the time grid by t_nd using a fixed-length parameter-weighted
# convolution, so no integer index ever depends on a parameter.
shift_density <- function(p, t_nd, dt, wmax) {
  "[<-" <- ADoverload("[<-")
  n <- length(p)
  s <- t_nd / dt
  out <- numeric(n) * t_nd * 0
  for (j in 0:wmax) {
    w <- b3_kernel(s - j)
    out[(j + 1L):n] <- out[(j + 1L):n] + w * p[1:(n - j)]
  }
  out
}

# ---- simulate GDDM data at Roitman dimensions --------------------------------------

sim_gddm <- function(n_per, p, dt = 2e-4, seed = 7) {
  set.seed(seed)
  out <- list()
  for (cc in COH) {
    drift <- p$mu0 * (cc / CMAX)^p$alpha
    n <- n_per
    x <- numeric(n); alive <- rep(TRUE, n)
    rt <- rep(NA_real_, n); resp <- rep(NA_integer_, n)
    for (k in seq_len(as.integer(T_DUR / dt))) {
      t <- k * dt
      idx <- which(alive); if (!length(idx)) break
      x[idx] <- x[idx] + (drift - p$leak * x[idx]) * dt + sqrt(dt) * rnorm(length(idx))
      B <- p$B0 * exp(-t / p$tau)
      hu <- idx[x[idx] >= B]; hl <- idx[x[idx] <= -B]
      if (length(hu)) { rt[hu] <- t; resp[hu] <- 1L; alive[hu] <- FALSE }
      if (length(hl)) { rt[hl] <- t; resp[hl] <- 0L; alive[hl] <- FALSE }
    }
    rt <- rt + p$t_nd
    # uniform contaminant, exactly the paper's 5 percent split evenly across responses
    lap <- runif(n) < LAPSE
    rt[lap] <- runif(sum(lap), 0, T_DUR)
    resp[lap] <- rbinom(sum(lap), 1, 0.5)
    out[[length(out) + 1L]] <- data.frame(coh = cc, rt = rt, resp = resp)
  }
  d <- do.call(rbind, out)
  d[!is.na(d$rt) & d$rt < T_DUR, ]
}

# ---- build the aggregated likelihood ------------------------------------------------
# The density depends on the trial only through its condition and its RT bin, so trials
# collapse to counted cells. This is the same economy PyDDM gets from solving once per
# condition, and it keeps the tape small.
make_obj <- function(dat, nt, ny) {
  dt <- T_DUR / nt
  dat$bin <- pmin(pmax(as.integer(round(dat$rt / dt)), 1L), nt)
  cells <- aggregate(list(n = rep(1, nrow(dat))),
                     by = list(coh = dat$coh, bin = dat$bin, resp = dat$resp), FUN = sum)
  wmax <- as.integer(ceiling(NDT_MAX / dt)) + 2L
  ci <- match(cells$coh, COH)

  f <- function(par) {
    getAll(par)
    mu0 <- exp(lmu0); alpha <- exp(lalpha); B0 <- exp(lB0)
    tau <- exp(ltau); t_nd <- NDT_MAX * plogis(qt_nd); leak <- leak_raw
    nll <- 0
    for (j in seq_along(COH)) {
      # d/d(alpha) of C^alpha is C^alpha log(C), which is NaN at C = 0, and the Roitman
      # design has a 0 percent coherence condition. PyDDM never trips over this because
      # differential evolution takes no derivatives. COH[j] is data, so branching on it
      # here happens at tape-build time and is legal.
      drift <- if (COH[j] == 0) 0 * mu0 else mu0 * (COH[j] / CMAX)^alpha
      s <- gddm_solve(drift, leak, B0, tau, nt, ny, dt)
      pu <- shift_density(s$up, t_nd, dt, wmax)
      pl <- shift_density(s$lo, t_nd, dt, wmax)
      # paper Eq. 13 overlay: 95 percent process, 5 percent uniform split evenly
      pu <- (1 - LAPSE) * pu + LAPSE * 0.5 / T_DUR
      pl <- (1 - LAPSE) * pl + LAPSE * 0.5 / T_DUR
      k <- which(ci == j)
      if (length(k)) {
        ku <- k[cells$resp[k] == 1L]; kl <- k[cells$resp[k] == 0L]
        if (length(ku)) nll <- nll - sum(cells$n[ku] * log(pu[cells$bin[ku] + 1L]))
        if (length(kl)) nll <- nll - sum(cells$n[kl] * log(pl[cells$bin[kl] + 1L]))
      }
    }
    nll
  }
  list(f = f, ncell = nrow(cells), dt = dt)
}
