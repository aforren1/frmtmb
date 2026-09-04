# Probe 1: is the GDDM first-passage density reachable from RTMB at all?
#
# PyDDM solves the Fokker-Planck equation on a grid whose spatial extent shrinks with the
# collapsing bound, and interpolates between two integer grids when the bound falls
# between nodes. That makes its objective continuous but kinked in B0 and tau, which is
# why the paper reaches for differential evolution. A gradient-based host like frmtmb
# needs a formulation that is smooth in every parameter.
#
# The fix is a change of variable. With y = x / B(t) the absorbing boundaries sit at
# +/-1 for all t, so the grid is fixed and no index depends on a parameter:
#
#   dx = (mu0 f(C) - l x) dt + dW,   B(t) = B0 exp(-r t)
#   y  = x / B(t)
#   dy = [mu0 f(C) / B(t) - l y + r y] dt + (1 / B(t)) dW
#
# Every coefficient is then a smooth function of the parameters, and the only remaining
# non-smoothness is the non-decision-time lookup, handled by interpolation.
#
# This script contains no frmtmb code. It only establishes that the density is
# computable and correct, so that the feasibility verdict rests on measurement.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Crank-Nicolson solve of the Fokker-Planck equation in the rescaled coordinate.
# Returns defective densities at the upper and lower boundary on the time grid.
gddm_fpt <- function(mu, leak, B0, rate, T_dur = 2, dt = 0.005, ny = 201,
                     ops = list()) {
  # ops lets the caller swap in AD-aware primitives later; base R by default.
  exp_ <- ops$exp %||% exp

  h <- 2 / (ny + 1)                 # interior nodes only; boundaries carry zero density
  y <- seq(-1 + h, 1 - h, length.out = ny)
  nt <- as.integer(round(T_dur / dt))

  rho <- numeric(ny)
  i0 <- (ny + 1L) %/% 2L             # ny odd puts a node exactly at y = 0
  rho[i0] <- 1 / h                   # delta at x = 0, normalized on the grid

  p_up <- numeric(nt + 1L)
  p_lo <- numeric(nt + 1L)

  coefs <- function(t) {
    Bt <- B0 * exp_(-rate * t)
    list(a = mu / Bt - leak * y + rate * y, D = 0.5 / (Bt * Bt))
  }

  tri_mult <- function(lo, di, up, v) {
    out <- di * v
    out[-ny] <- out[-ny] + up[-ny] * v[-1L]
    out[-1L] <- out[-1L] + lo[-1L] * v[-ny]
    out
  }

  # Thomas algorithm. Written as an explicit scalar sweep because that is what an AD
  # tape can record; solve() on a dense matrix would work too but costs O(ny^3).
  tri_solve <- function(lo, di, up, rhs) {
    cp <- numeric(ny); dp <- numeric(ny)
    cp[1L] <- up[1L] / di[1L]
    dp[1L] <- rhs[1L] / di[1L]
    for (i in 2:ny) {
      m <- di[i] - lo[i] * cp[i - 1L]
      cp[i] <- if (i < ny) up[i] / m else 0
      dp[i] <- (rhs[i] - lo[i] * dp[i - 1L]) / m
    }
    x <- numeric(ny)
    x[ny] <- dp[ny]
    for (i in (ny - 1L):1L) x[i] <- dp[i] - cp[i] * x[i + 1L]
    x
  }

  build <- function(cf) {
    a <- cf$a; D <- cf$D
    lo <- c(0, (a[-ny] / (2 * h) + D / h^2)[-0])
    lo <- c(0, a[seq_len(ny - 1L)] / (2 * h) + D / h^2)
    up <- c(-a[2:ny] / (2 * h) + D / h^2, 0)
    di <- rep(-2 * D / h^2, ny)
    list(lo = lo, di = di, up = up)
  }

  for (k in seq_len(nt)) {
    t0 <- (k - 1) * dt
    t1 <- k * dt
    c0 <- coefs(t0); c1 <- coefs(t1)
    M0 <- build(c0);  M1 <- build(c1)

    rhs <- rho + (dt / 2) * tri_mult(M0$lo, M0$di, M0$up, rho)
    new <- tri_solve(-(dt / 2) * M1$lo, 1 - (dt / 2) * M1$di, -(dt / 2) * M1$up, rhs)

    # Absorbed flux. With zero density at the wall the advective term drops out and the
    # flux is purely diffusive, D d(rho)/ds with s measured inward from the wall. A
    # one-sided first difference loses a percent or two of the absorbed mass at high
    # drift, so use the second-order stencil against the zero-density wall.
    fu <- function(r) (4 * r[ny] - r[ny - 1L]) / (2 * h)
    fl <- function(r) (4 * r[1L] - r[2L]) / (2 * h)
    p_up[k + 1L] <- 0.5 * (c0$D * fu(rho) + c1$D * fu(new))
    p_lo[k + 1L] <- 0.5 * (c0$D * fl(rho) + c1$D * fl(new))
    rho <- new
  }

  list(t = seq(0, T_dur, by = dt), up = p_up, lo = p_lo, dt = dt)
}

# ---- validation against closed-form DDM identities -------------------------------
# For constant drift v, symmetric bounds +/-B, unit noise, unbiased start:
#   P(upper) = 1 / (1 + exp(-2 v B))
#   E[decision time] = (B / v) tanh(v B)

cat("== constant-bound, no-leak checks (analytic targets) ==\n")
grid <- expand.grid(v = c(0.5, 1, 2), B = c(0.8, 1.5))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  v <- grid$v[i]; B <- grid$B[i]
  s <- gddm_fpt(mu = v, leak = 0, B0 = B, rate = 0, T_dur = 12, dt = 0.002, ny = 301)
  p_up_num <- sum(s$up) * s$dt
  p_lo_num <- sum(s$lo) * s$dt
  mrt_num <- sum(s$t * (s$up + s$lo)) * s$dt / (p_up_num + p_lo_num)
  data.frame(v = v, B = B,
             p_up = p_up_num, p_up_exact = 1 / (1 + exp(-2 * v * B)),
             mrt = mrt_num, mrt_exact = (B / v) * tanh(v * B),
             mass = p_up_num + p_lo_num)
}))
res[] <- lapply(res, function(x) if (is.numeric(x)) round(x, 5) else x)
print(res)
cat("max |p_up error| =", max(abs(res$p_up - res$p_up_exact)), "\n")
cat("max |mean RT error| =", max(abs(res$mrt - res$mrt_exact)), "\n\n")

# ---- validation of the leaky + collapsing case against Euler-Maruyama ------------
cat("== leaky + exponentially collapsing bound vs Monte Carlo ==\n")
set.seed(1)
mc_gddm <- function(mu, leak, B0, rate, T_dur = 2, dt = 2e-4, n = 2e5) {
  x <- numeric(n); alive <- rep(TRUE, n)
  rt <- rep(NA_real_, n); resp <- rep(NA_integer_, n)
  nt <- as.integer(T_dur / dt)
  sq <- sqrt(dt)
  for (k in seq_len(nt)) {
    t <- k * dt
    idx <- which(alive)
    if (!length(idx)) break
    x[idx] <- x[idx] + (mu - leak * x[idx]) * dt + sq * rnorm(length(idx))
    Bt <- B0 * exp(-rate * t)
    hit_u <- idx[x[idx] >= Bt]; hit_l <- idx[x[idx] <= -Bt]
    if (length(hit_u)) { rt[hit_u] <- t; resp[hit_u] <- 1L; alive[hit_u] <- FALSE }
    if (length(hit_l)) { rt[hit_l] <- t; resp[hit_l] <- 0L; alive[hit_l] <- FALSE }
  }
  data.frame(rt = rt, resp = resp)
}

pars <- list(mu = 10.49091 * 0.128, leak = -0.482, B0 = 1.811, rate = 1.992)
s <- do.call(gddm_fpt, c(pars, list(T_dur = 2, dt = 0.002, ny = 301)))
sim <- do.call(mc_gddm, c(pars, list(T_dur = 2)))

cat(sprintf("solver: P(upper)=%.4f P(lower)=%.4f mean RT|decided=%.4f\n",
            sum(s$up) * s$dt, sum(s$lo) * s$dt,
            sum(s$t * (s$up + s$lo)) * s$dt / (sum(s$up + s$lo) * s$dt)))
cat(sprintf("  MC  : P(upper)=%.4f P(lower)=%.4f mean RT|decided=%.4f  (n=%d)\n",
            mean(sim$resp == 1, na.rm = TRUE) * mean(!is.na(sim$resp)),
            mean(sim$resp == 0, na.rm = TRUE) * mean(!is.na(sim$resp)),
            mean(sim$rt, na.rm = TRUE), nrow(sim)))

qs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
cdf_u <- cumsum(s$up) * s$dt / (sum(s$up) * s$dt)
solver_q <- approx(cdf_u, s$t, xout = qs)$y
mc_q <- quantile(sim$rt[which(sim$resp == 1)], qs, na.rm = TRUE)
cat("upper-boundary RT quantiles\n")
print(round(rbind(solver = solver_q, montecarlo = as.numeric(mc_q)), 4))
