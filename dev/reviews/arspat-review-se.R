# Review of lane arspat: se() beside the ar() family. Draws against the
# summed covariance (gaussian, student, heteroscedastic), the pointwise
# log-density against a dense leave-one-out conditional, refusals on the
# latent-residual families, the nu -> large gap. Output to stdout.
source("dev/reviews/arspat-review-common.R")

mk <- function(seed = 7, G = 40, K = 5, tdf = NULL) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  n <- nrow(d)
  d$x <- rnorm(n)
  d$x2 <- rnorm(n)
  d$s <- runif(n, 0.2, 0.9)
  e <- as.vector(vapply(seq_len(G), function(i)
    as.vector(arima.sim(list(ar = 0.6), K, sd = 0.6)), numeric(K)))
  if (!is.null(tdf)) e <- e * rep(sqrt(tdf / rchisq(G, tdf)), each = K)
  d$y <- 1 + 0.6 * d$x + e * exp(0.3 * d$x2) + rnorm(n, 0, d$s)
  d$yp <- exp(d$y / 3)
  d
}

group_S <- function(fit, d, i) {
  dp <- ns$eval_dpars(fit)[[1L]]
  sig <- rep_len(dp$sigma, nrow(d))
  nat <- ac_nat(fit)
  ar <- unname(nat[grep("^ar", names(nat))])
  ma <- unname(nat[grep("^ma", names(nat))])
  gm <- arma_gamma(ar, ma, 1, length(i))
  tt <- d$t[i]
  matrix(gm[abs(outer(tt, tt, "-")) + 1L], length(i)) *
    outer(sig[i], sig[i]) + diag(d$s[i]^2, length(i))
}

# empirical covariance of simulate() draws of the first groups against S
sim_check <- function(fit, d, nsim = 20000, scale = 1) {
  sims <- as.matrix(simulate(fit, nsim = nsim, seed = 3))
  mu <- rep_len(ns$eval_dpars(fit)[[1L]]$mu, nrow(d))
  worst <- 0
  for (g in 1:4) {
    i <- which(as.integer(d$g) == g)
    S <- group_S(fit, d, i) * scale
    C <- stats::cov(t(sims[i, ]))
    # the entrywise error in units of S's diagonal scale
    worst <- max(worst, max(abs(C - S) / sqrt(outer(diag(S), diag(S)))))
    if (g == 1) {
      cat("group 1 S diag:", format(diag(S), digits = 4), "\n")
      cat("group 1 C diag:", format(diag(C), digits = 4), "\n")
      cat("group 1 S[1,2] C[1,2]:", format(c(S[1, 2], C[1, 2]), digits = 4),
          "\n")
      cat("mean offset max:", format(max(abs(rowMeans(sims[i, ]) - mu[i])),
                                     digits = 3), "\n")
    }
  }
  cat(sprintf("max |C - S| / sqrt(S_ii S_jj) over 4 groups: %.4f (nsim %d)\n",
              worst, nsim))
}

d <- mk()
chk("gaussian het se: simulate covariance", {
  fit <- frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE),
                sigma ~ x2), data = d)
  sim_check(fit, d)
  # predict(): Est.Error against sqrt(diag(S))
  p <- predict(fit)
  sdr <- unlist(lapply(1:4, function(g) {
    i <- which(as.integer(d$g) == g)
    sqrt(diag(group_S(fit, d, i)))
  }))
  cat("predict Est.Error / sqrt(S_ii), first 20 rows: range",
      format(range(p[1:20, "Est.Error"] / sdr), digits = 4), "\n")
  # pointwise log-density against brms's leave-one-out conditional
  ll <- frmtmb::autocor_row_loglik(fit, ns$eval_dpars(fit), names(fit$spec$responses)[1])
  mu <- rep_len(ns$eval_dpars(fit)[[1L]]$mu, nrow(d))
  ref <- unlist(lapply(split(seq_len(nrow(d)), d$g), function(i) {
    S <- group_S(fit, d, i)
    Q <- solve(S)
    e <- d$y[i] - mu[i]
    g <- as.vector(Q %*% e)
    dnorm(d$y[i], d$y[i] - g / diag(Q), sqrt(1 / diag(Q)), log = TRUE)
  }))
  cat(sprintf("log_lik pointwise max abs diff to dense LOO: %.3e\n",
              max(abs(ll - ref))))
  NULL
})

dt <- mk(seed = 8, tdf = 6)
chk("student se: simulate covariance is nu/(nu-2) S", {
  fit <- frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)),
             family = student(), data = dt)
  nu <- ns$eval_dpars(fit)[[1L]]$nu[1]
  cat("nu", nu, "\n")
  if (nu > 2.5) sim_check(fit, dt, scale = nu / (nu - 2))
  NULL
})

chk("se + ar on skew_normal", frm(bf(y | se(s, sigma = TRUE) ~ x +
                                     ar(t, g, cov = TRUE)),
                                  family = skew_normal(), data = d))
chk("se + ar on lognormal", frm(bf(yp | se(s, sigma = TRUE) ~ x +
                                     ar(t, g, cov = TRUE)),
                                  family = lognormal(), data = d))
chk("se + ar cov = FALSE", frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g)),
                               data = d))
chk("se + cosy sigma = FALSE", frm(bf(y | se(s) ~ x + cosy(t, g)), data = d))
chk("se + ar + mixture refused?", frm(bf(y | se(s, sigma = TRUE) ~ x +
                                         ar(t, g, cov = TRUE)),
                                      family = mixture(gaussian, gaussian),
                                      data = d))

# simulate() with newdata under se(): reads newdata's s
chk("simulate newdata se", {
  fit <- frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)), data = d)
  nd <- d[d$g %in% c("1", "2"), ]
  nd$s <- 5
  s1 <- as.matrix(simulate(fit, newdata = nd, nsim = 4000, seed = 1))
  cat("sd of draws with s = 5 (expect >= 5):",
      format(range(apply(s1, 1, sd)), digits = 4), "\n")
  nd$s <- NULL
  simulate(fit, newdata = nd, nsim = 2, seed = 1)
})

# the nu -> large gap: lgamma difference rounding at nu = 1e9
chk("lgamma rounding at large nu", {
  for (nu in c(1e5, 1e7, 1e9)) {
    for (k in c(3, 6)) {
      naive <- lgamma((nu + k) / 2) - lgamma(nu / 2)
      shift <- ns$lgamma_shift_diff(nu / 2, k / 2)
      # reference: sum_{j} log(nu/2 + j) for integer k/2, else Rmpfr-free
      # asymptotic: k/2 log(nu/2) + (k/2)(k/2 - 1)/(nu) ... via digamma
      a <- nu / 2; b <- k / 2
      ref <- b * log(a) + b * (b - 1) / (2 * a) -
        b * (b - 1) * (2 * b - 1) / (12 * a^2)
      cat(sprintf("nu %.0e k %d naive-ref %.3e shift-ref %.3e\n", nu, k,
                  naive - ref, shift - ref))
    }
  }
  NULL
})
