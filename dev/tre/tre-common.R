# Shared helpers for the t-distributed random-effects feasibility probes.
#
# The model everywhere below is the simplest one that isolates the open
# question: a gaussian LMM with a SCALAR random intercept whose latent
# density is a scaled Student-t. A scalar latent is what makes an exact
# reference possible: the marginal likelihood of each group is a
# one-dimensional integral, so adaptive quadrature is the truth and the
# Laplace approximation can be measured against it rather than against
# another approximation.
#
# Parameterization matches brms's gr(dist = "student") (see
# dev/tre-feasibility.md section 2): b_j ~ t_nu(0, s), where `s` is the
# t's SCALE, not its standard deviation. sd(b) = s * sqrt(nu/(nu-2)).

suppressMessages({
  library(RTMB)
})

# ---------------------------------------------------------------- data

#' Simulate a gaussian LMM with t-distributed random intercepts.
#'
#' @param G number of groups, `n` observations per group.
#' @param nu latent degrees of freedom; `Inf` draws gaussian latents.
#' @param s latent SCALE (not SD).
sim_tre <- function(G = 40, n = 8, beta = c(1, 0.5), sigma = 1, s = 1,
                    nu = 3, seed = 1, b_override = NULL) {
  set.seed(seed)
  x <- rnorm(G * n)
  g <- rep(seq_len(G), each = n)
  b <- if (!is.null(b_override)) {
    b_override
  } else if (is.finite(nu)) {
    s * rt(G, df = nu)
  } else {
    s * rnorm(G)
  }
  X <- cbind(1, x)
  y <- as.vector(X %*% beta) + b[g] + rnorm(G * n, 0, sigma)
  list(y = y, X = X, x = x, g = g, G = G, n = n, b = b,
       truth = c(beta = beta, sigma = sigma, s = s, nu = nu))
}

#' Per-group sufficient statistics of the residuals for a given beta.
#'
#' The conditional log-likelihood of group j is quadratic in the latent,
#' so `n_j`, `sum(r)` and `sum(r^2)` are all the exact reference needs.
suff_stats <- function(dat, beta) {
  r <- dat$y - as.vector(dat$X %*% beta)
  list(n = as.vector(tapply(r, dat$g, length)),
       S = as.vector(tapply(r, dat$g, sum)),
       SS = as.vector(tapply(r, dat$g, function(z) sum(z^2))))
}

# ------------------------------------------------- the exact reference

#' log density of a scaled t, vectorized, no normalizing shortcuts.
ldt <- function(b, s, nu) {
  if (!is.finite(nu)) return(dnorm(b, 0, s, log = TRUE))
  stats::dt(b / s, df = nu, log = TRUE) - log(s)
}

#' d/db of `ldt`.
ldt_grad <- function(b, s, nu) {
  if (!is.finite(nu)) return(-b / s^2)
  -(nu + 1) * b / (s^2 * nu + b^2)
}

#' d2/db2 of `ldt`.
ldt_hess <- function(b, s, nu) {
  if (!is.finite(nu)) return(rep(-1 / s^2, length(b)))
  d <- s^2 * nu + b^2
  -(nu + 1) * (d - 2 * b^2) / d^2
}

#' The per-group log integrand h_j(b) = log p(y_j | b) + log p(b).
#'
#' Vectorized over a G x K matrix of latent values (one row per group).
h_mat <- function(B, st, sigma, s, nu) {
  n <- st$n
  -(st$SS - 2 * B * st$S + n * B^2) / (2 * sigma^2) -
    n / 2 * log(2 * pi * sigma^2) + ldt(B, s, nu)
}

#' Global mode and curvature of h_j, per group.
#'
#' Two Newton starts (the latent-free residual mean, and zero) cover the
#' bimodality a conflicting outlier group produces: with a t latent the
#' prior mode at 0 and the data mode near `S/n` can both be local
#' maxima. `all_modes` reports whether the two starts disagreed.
h_mode <- function(st, sigma, s, nu) {
  n <- st$n
  grad <- function(b) -(n * b - st$S) / sigma^2 + ldt_grad(b, s, nu)
  hess <- function(b) -n / sigma^2 + ldt_hess(b, s, nu)
  newton <- function(b0) {
    b <- b0
    for (it in 1:200) {
      gh <- hess(b)
      # damp away from a non-concave region rather than stepping uphill
      gh <- pmin(gh, -1e-8)
      step <- -grad(b) / gh
      step <- pmax(pmin(step, 5 * (abs(b) + 1)), -5 * (abs(b) + 1))
      bn <- b + step
      if (max(abs(bn - b)) < 1e-12) { b <- bn; break }
      b <- bn
    }
    b
  }
  m1 <- newton(st$S / n)
  m2 <- newton(rep(0, length(n)))
  h1 <- h_mat(cbind(m1), st, sigma, s, nu)[, 1]
  h2 <- h_mat(cbind(m2), st, sigma, s, nu)[, 1]
  pick <- ifelse(h1 >= h2, 1, 2)
  m <- ifelse(pick == 1, m1, m2)
  list(mode = m, other = ifelse(pick == 1, m2, m1),
       curv = -hess(m),
       bimodal = abs(m1 - m2) > 1e-6,
       h_gap = abs(h1 - h2))
}

#' Adaptive Gauss-Hermite marginal log-likelihood, exact reference.
#'
#' The latent density is NOT gaussian, so the rule is the general
#' adaptive one: center and scale at the conditional mode and curvature
#' and undo the gaussian weight with exp(z^2). Convergence in `K` is
#' the thing to watch with fat tails, which is why every probe reports
#' more than one `K`.
aghq_loglik <- function(dat, beta, sigma, s, nu, K = 101, st = NULL,
                        modes = NULL) {
  if (is.null(st)) st <- suff_stats(dat, beta)
  if (is.null(modes)) modes <- h_mode(st, sigma, s, nu)
  gh <- statmod::gauss.quad(K, kind = "hermite")
  sdj <- 1 / sqrt(pmax(modes$curv, 1e-10))
  B <- modes$mode + outer(sdj, sqrt(2) * gh$nodes)
  H <- h_mat(B, st, sigma, s, nu)
  lw <- rep(log(gh$weights) + gh$nodes^2, each = length(sdj))
  dim(lw) <- dim(H)
  L <- H + lw + log(sqrt(2) * sdj)
  mx <- apply(L, 1, max)
  sum(mx + log(rowSums(exp(L - mx))))
}

#' Gold-standard per-group marginal log-likelihood by adaptive
#' quadrature over the whole line (R's `integrate`, which does the
#' tan-substitution). Slow; used to certify `aghq_loglik`.
integrate_loglik <- function(dat, beta, sigma, s, nu, st = NULL) {
  if (is.null(st)) st <- suff_stats(dat, beta)
  modes <- h_mode(st, sigma, s, nu)
  tot <- 0
  for (j in seq_along(st$n)) {
    stj <- list(n = st$n[j], S = st$S[j], SS = st$SS[j])
    hmax <- h_mat(cbind(modes$mode[j]), stj, sigma, s, nu)[1, 1]
    f <- function(b) exp(h_mat(cbind(b), stj, sigma, s, nu)[, 1] - hmax)
    v <- stats::integrate(f, -Inf, Inf, rel.tol = 1e-12,
                          subdivisions = 2000L)
    tot <- tot + hmax + log(v$value)
  }
  tot
}

#' The Laplace approximation computed by hand, so its error can be
#' attributed to the approximation and not to TMB.
laplace_loglik <- function(dat, beta, sigma, s, nu, st = NULL) {
  if (is.null(st)) st <- suff_stats(dat, beta)
  modes <- h_mode(st, sigma, s, nu)
  hmax <- h_mat(cbind(modes$mode), st, sigma, s, nu)[, 1]
  sum(hmax + 0.5 * log(2 * pi) - 0.5 * log(modes$curv))
}

# ------------------------------------------------------- the two fits

#' Exact ML fit: maximize the AGHQ marginal likelihood directly.
fit_exact <- function(dat, nu, K = 101, start = NULL, fit_nu = FALSE) {
  # unname: nlminb propagates the start's names onto op$par, and the
  # named c() below would then build "beta0.beta0"
  p0 <- unname(start %||% c(0, 0, 0, 0))
  nll <- function(p) {
    beta <- p[1:2]
    sigma <- exp(p[3])
    s <- exp(p[4])
    nuu <- if (fit_nu) 2 + exp(p[5]) else nu
    -aghq_loglik(dat, beta, sigma, s, nuu, K = K)
  }
  if (fit_nu) p0 <- c(p0, log(nu - 2))
  op <- stats::nlminb(p0, nll, control = list(eval.max = 2000,
                                              iter.max = 1000,
                                              rel.tol = 1e-12))
  list(par = c(beta0 = op$par[1], beta1 = op$par[2],
               sigma = exp(op$par[3]), s = exp(op$par[4]),
               nu = if (fit_nu) 2 + exp(op$par[5]) else nu),
       loglik = -op$objective, conv = op$convergence, opt = op)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Laplace fit through RTMB, `random = "b"`.
#'
#' @param nu fixed degrees of freedom, or NA to estimate `nu` too
#'   (parameterized as `nu = lower + exp(log_nu_free)`).
fit_laplace <- function(dat, nu = 3, start = NULL, estimate_nu = FALSE,
                        nu_lower = 2, silent = TRUE, inner.control = list()) {
  ddat <- list(y = dat$y, X = dat$X, g = dat$g, nu_fix = nu,
               nu_lower = nu_lower, est = as.integer(estimate_nu))
  par <- list(beta = start$beta %||% c(0, 0),
              log_sigma = start$log_sigma %||% 0,
              log_s = start$log_s %||% 0,
              log_nu = start$log_nu %||% log(max(nu - nu_lower, 1)),
              b = rep(0, dat$G))
  fn <- function(p) {
    getAll(ddat, p)
    nuu <- if (est == 1L) nu_lower + exp(log_nu) else nu_fix
    sd_e <- exp(log_sigma)
    s <- exp(log_s)
    mu <- X %*% beta + b[g]
    nll <- -sum(dnorm(y, mu, sd_e, log = TRUE))
    nll <- nll - sum(dt(b / s, nuu, log = TRUE) - log_s)
    nll
  }
  map <- if (estimate_nu) list() else list(log_nu = factor(NA))
  obj <- MakeADFun(fn, par, random = "b", map = map, silent = silent,
                   inner.control = inner.control)
  op <- stats::nlminb(obj$par, obj$fn, obj$gr,
                      control = list(eval.max = 2000, iter.max = 1000,
                                     rel.tol = 1e-12))
  pl <- obj$env$parList(op$par)
  list(par = c(beta0 = pl$beta[1], beta1 = pl$beta[2],
               sigma = exp(pl$log_sigma), s = exp(pl$log_s),
               nu = if (estimate_nu) nu_lower + exp(pl$log_nu) else nu),
       b = obj$env$parList(op$par)$b,
       loglik = -op$objective, conv = op$convergence, obj = obj, opt = op)
}

#' Gaussian-latent Laplace fit, for the robustness comparison.
fit_gaussian <- function(dat, start = NULL) {
  ddat <- list(y = dat$y, X = dat$X, g = dat$g)
  par <- list(beta = c(0, 0), log_sigma = 0, log_s = 0,
              b = rep(0, dat$G))
  fn <- function(p) {
    getAll(ddat, p)
    mu <- X %*% beta + b[g]
    -sum(dnorm(y, mu, exp(log_sigma), log = TRUE)) -
      sum(dnorm(b, 0, exp(log_s), log = TRUE))
  }
  obj <- MakeADFun(fn, par, random = "b", silent = TRUE)
  op <- stats::nlminb(obj$par, obj$fn, obj$gr,
                      control = list(eval.max = 2000, iter.max = 1000))
  pl <- obj$env$parList(op$par)
  list(par = c(beta0 = pl$beta[1], beta1 = pl$beta[2],
               sigma = exp(pl$log_sigma), s = exp(pl$log_s), nu = Inf),
       b = pl$b, loglik = -op$objective, conv = op$convergence, obj = obj)
}

fmt <- function(x, d = 6) formatC(x, format = "f", digits = d)
