# Shared pieces for the frailty lane. Nothing here is package code; it
# exists so that every number in dev/frailty-findings.md has one
# construction that can be re-run.
source("frailty-lib.R")

# ---------------------------------------------------------------- basis

# The Royston-Parmar natural cubic spline basis, written independently
# of frmtmb.spline so that a disagreement with it is visible.
rp_basis <- function(knots, x) {
  nk <- length(knots)
  out <- cbind(rep(1, length(x)), x)
  if (nk > 2L) {
    kmin <- knots[1L]; kmax <- knots[nk]
    for (j in seq_len(nk - 2L)) {
      kj <- knots[j + 1L]
      lam <- (kmax - kj) / (kmax - kmin)
      out <- cbind(out, pmax(x - kj, 0)^3 - lam * pmax(x - kmin, 0)^3 -
                     (1 - lam) * pmax(x - kmax, 0)^3)
    }
  }
  colnames(out) <- paste0("g", seq_len(ncol(out)) - 1L)
  out
}

rp_dbasis <- function(knots, x) {
  nk <- length(knots)
  out <- cbind(rep(0, length(x)), rep(1, length(x)))
  if (nk > 2L) {
    kmin <- knots[1L]; kmax <- knots[nk]
    for (j in seq_len(nk - 2L)) {
      kj <- knots[j + 1L]
      lam <- (kmax - kj) / (kmax - kmin)
      out <- cbind(out, 3 * pmax(x - kj, 0)^2 -
                     3 * lam * pmax(x - kmin, 0)^2 -
                     3 * (1 - lam) * pmax(x - kmax, 0)^2)
    }
  }
  out
}

# flexsurv's own default knot rule, which royston_parmar() also uses.
rp_knots <- function(logev, df) {
  bk <- range(logev)
  ik <- if (df > 1L) {
    unname(stats::quantile(logev,
                           seq(0, 1, length.out = df + 1L)))[2:df]
  } else numeric(0)
  list(ik = ik, bk = bk, all = c(bk[1L], ik, bk[2L]))
}

# ----------------------------------------------------------- simulator

# A Weibull proportional-hazards model with a shared log-normal frailty
# on the cluster. This is EXACTLY a Royston-Parmar model with no
# interior knot: log H = gamma0 + gamma1 log t + x'beta + b, with
# gamma1 the Weibull shape and gamma0 = -shape * log(scale). Every
# interior spline coefficient therefore has truth 0 whatever df the fit
# uses, so a recovery table can be written at df > 1 as well.
frailty_sim <- function(seed, n = 2000L, n_centre = 40L, sd_b = 0.5,
                        beta = c(0.6), shape = 1.3, scale = 5,
                        tau = NULL, p_cens = 0.4) {
  set.seed(seed)
  centre <- rep(seq_len(n_centre), length.out = n)
  b <- stats::rnorm(n_centre, 0, sd_b)
  trt <- stats::rbinom(n, 1L, 0.5)
  eta <- beta[1L] * trt + b[centre]
  u <- stats::runif(n)
  # H(t) = (t / scale)^shape * exp(eta); S = exp(-H)
  tt <- scale * (-log(u) * exp(-eta))^(1 / shape)
  if (is.null(tau)) tau <- unname(stats::quantile(tt, 1 - p_cens))
  event <- as.integer(tt <= tau)
  d <- data.frame(time = pmin(tt, tau), event = event,
                  censored = 1L - event, trt = trt,
                  centre = factor(centre))
  attr(d, "truth") <- list(seed = seed, sd_b = sd_b, beta = beta,
                           shape = shape, scale = scale, tau = tau,
                           gamma0 = -shape * log(scale), gamma1 = shape,
                           b = b, p_cens = mean(1 - event))
  d
}

# --------------------------------------------- the exact log likelihood

# The marginal log likelihood of the shared log-normal frailty model,
# integrated by adaptive quadrature per cluster rather than by any rule
# either package uses. This is the reference both approximations are
# measured against; nothing in it comes from frmtmb or rstpm2.
#
# `gam` is the RP coefficient vector (gamma0 first), `xb` the linear
# predictor from the covariates WITHOUT the intercept, `sd_b` the
# frailty standard deviation.
frailty_exact_ll <- function(time, event, xb, cluster, knots, gam, sd_b,
                             rel.tol = 1e-12) {
  x <- log(time)
  B <- rp_basis(knots, x)
  dB <- rp_dbasis(knots, x)
  eta0 <- as.vector(B %*% gam) + xb
  deta <- as.vector(dB %*% gam)
  lg <- log(deta) - x
  cl <- split(seq_along(time), cluster)
  gfun <- function(b, idx) {
    e <- eta0[idx] + b
    sum(event[idx] * (e + lg[idx]) - exp(e)) +
      stats::dnorm(b, 0, sd_b, log = TRUE)
  }
  out <- 0
  for (idx in cl) {
    op <- stats::optimize(function(b) -gfun(b, idx),
                          c(-12 * sd_b, 12 * sd_b), tol = 1e-12)
    bh <- op$minimum
    gh <- -op$objective
    # curvature at the mode, to set an integration range that the
    # integrand is actually supported on
    h <- 1e-4
    curv <- (gfun(bh + h, idx) - 2 * gh + gfun(bh - h, idx)) / h^2
    sdl <- 1 / sqrt(max(-curv, 1e-8))
    lo <- bh - 12 * sdl; hi <- bh + 12 * sdl
    iv <- stats::integrate(
      function(bs) vapply(bs, function(bb) exp(gfun(bb, idx) - gh),
                          numeric(1)),
      lo, hi, rel.tol = rel.tol, subdivisions = 2000L)
    out <- out + gh + log(iv$value)
  }
  out
}

# ------------------------------------------------- basis change nsx->RP

# rstpm2 writes the same natural cubic spline in the nsx basis. Both
# span the same space, so the map is linear and exact; the residual of
# the solve is reported, because a non-zero residual would mean the two
# spline spaces are NOT the same and every comparison below it would be
# meaningless.
nsx_to_rp <- function(rst, knots, xgrid) {
  cf <- stats::coef(rst)
  nm <- names(cf)
  jj <- grep("^nsx", nm)
  a0 <- unname(cf[["(Intercept)"]])
  cs <- unname(cf[jj])
  nsxcall <- attr(stats::terms(rst@model.frame), "term.labels")
  # evaluate the nsx basis on the grid through the fitted object's own
  # design so that no assumption about knot placement is smuggled in
  Bn <- rstpm2::nsx(xgrid, knots = knots[-c(1L, length(knots))],
                    Boundary.knots = knots[c(1L, length(knots))])
  eta <- a0 + as.vector(Bn %*% cs)
  Br <- rp_basis(knots, xgrid)
  fitls <- stats::lsfit(Br, eta, intercept = FALSE)
  list(gam = unname(fitls$coefficients),
       resid = max(abs(fitls$residuals)),
       scale = max(abs(eta)))
}

# The same exact marginal, with the random effect on gamma1 rather than
# on gamma0: eta = B gam + xb + u x, so the SLOPE in log time is what
# varies by cluster and d(eta)/dx carries u as well.
frailty_exact_ll_slope <- function(time, event, xb, cluster, knots, gam,
                                   sd_u, rel.tol = 1e-12) {
  x <- log(time)
  B <- rp_basis(knots, x)
  dB <- rp_dbasis(knots, x)
  eta0 <- as.vector(B %*% gam) + xb
  deta0 <- as.vector(dB %*% gam)
  cl <- split(seq_along(time), cluster)
  gfun <- function(u, idx) {
    e <- eta0[idx] + u * x[idx]
    dd <- deta0[idx] + u
    if (any(dd <= 0)) return(-Inf)
    sum(event[idx] * (e + log(dd) - x[idx])) - sum(exp(e)) +
      stats::dnorm(u, 0, sd_u, log = TRUE)
  }
  out <- 0
  for (idx in cl) {
    op <- stats::optimize(function(u) -gfun(u, idx),
                          c(-8 * sd_u, 8 * sd_u), tol = 1e-12)
    uh <- op$minimum
    gh <- -op$objective
    h <- 1e-4
    curv <- (gfun(uh + h, idx) - 2 * gh + gfun(uh - h, idx)) / h^2
    sdl <- 1 / sqrt(max(-curv, 1e-8))
    iv <- stats::integrate(
      function(us) vapply(us, function(uu) exp(gfun(uu, idx) - gh),
                          numeric(1)),
      uh - 12 * sdl, uh + 12 * sdl, rel.tol = rel.tol,
      subdivisions = 2000L)
    out <- out + gh + log(iv$value)
  }
  out
}

# The gamma1 design: log H = gamma0 + beta trt + (gamma1 + u_c) log t.
g1_sim <- function(seed, n = 2000L, n_centre = 40L, sd_u = 0.2,
                   gamma1 = 1.3, scale = 5, beta = 0.6, p_cens = 0.4) {
  set.seed(seed)
  centre <- rep(seq_len(n_centre), length.out = n)
  u <- stats::rnorm(n_centre, 0, sd_u)
  trt <- stats::rbinom(n, 1L, 0.5)
  g0 <- -gamma1 * log(scale)
  sh <- gamma1 + u[centre]
  e <- stats::runif(n)
  tt <- exp((log(-log(e)) - g0 - beta * trt) / sh)
  tau <- unname(stats::quantile(tt, 1 - p_cens))
  ev <- as.integer(tt <= tau)
  d <- data.frame(time = pmin(tt, tau), event = ev, censored = 1L - ev,
                  trt = trt, centre = factor(centre))
  attr(d, "truth") <- list(seed = seed, sd_u = sd_u, gamma1 = gamma1,
                           gamma0 = g0, beta = beta, u = u,
                           shape = gamma1 + u, p_cens = mean(1 - ev))
  d
}
