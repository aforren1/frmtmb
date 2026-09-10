# Lane `latent`, item 2.3: the simulator for the plan's hmm design.
#
# 50 sequences x 500 steps, K = 3 gaussian, one random intercept on the
# 1 -> 2 transition. This is the Phase 0 scale row's construction
# (extensions/frmtmb.latent/tests/testthat/test-scale.R,
# `latent_hmm_data()`), lifted here so a replicate can be drawn at an
# arbitrary seed and so the truth travels with the data.
#
# The transition matrix is a row-wise multinomial logit with STATE 1 as
# the reference cell in every row, which is frmtmb's convention and
# hmmTMB's when it is given `ref = rep(1, K)`.

hmm_truth <- list(
  mu = c(-2, 0, 3),
  sigma = 0.7,
  eta = rbind(c(-1.5, -2.5), c(2.0, -1.0), c(-1.0, 2.0)),
  sd_tr12 = 0.6)

hmm_tpm <- function(eta) {
  t(apply(cbind(0, eta), 1L, function(z) exp(z) / sum(exp(z))))
}

hmm_sim <- function(seed, ns = 50L, tl = 500L, tr = hmm_truth) {
  set.seed(seed)
  u <- stats::rnorm(ns, 0, tr$sd_tr12)
  out <- vector("list", ns)
  for (i in seq_len(ns)) {
    e <- tr$eta
    e[1L, 1L] <- e[1L, 1L] + u[i]
    G <- hmm_tpm(e)
    s <- integer(tl)
    s[1L] <- sample.int(3L, 1L)
    for (t in seq_len(tl)[-1L]) {
      s[t] <- sample.int(3L, 1L, prob = G[s[t - 1L], ])
    }
    out[[i]] <- data.frame(id = i, t = seq_len(tl), state = s,
                           y = stats::rnorm(tl, tr$mu[s], tr$sigma))
  }
  d <- do.call(rbind, out)
  d$id <- factor(d$id)
  list(d = d, u = u, truth = tr, seed = seed)
}

# The fixed-transition version of the same design: no random effect, so
# depmixS4 fits exactly this model.
hmm_sim_fixed <- function(seed, ns = 50L, tl = 500L, tr = hmm_truth) {
  set.seed(seed)
  G <- hmm_tpm(tr$eta)
  out <- vector("list", ns)
  for (i in seq_len(ns)) {
    s <- integer(tl)
    s[1L] <- sample.int(3L, 1L)
    for (t in seq_len(tl)[-1L]) {
      s[t] <- sample.int(3L, 1L, prob = G[s[t - 1L], ])
    }
    out[[i]] <- data.frame(id = i, t = seq_len(tl), state = s,
                           y = stats::rnorm(tl, tr$mu[s], tr$sigma))
  }
  d <- do.call(rbind, out)
  d$id <- factor(d$id)
  list(d = d, truth = tr, G = G, seed = seed)
}

# The states come back in whatever order the optimizer found them, so
# every comparison scores the best matching and says which it used.
hmm_perms <- function(K) {
  g <- as.matrix(expand.grid(rep(list(seq_len(K)), K)))
  g <- unname(g[apply(g, 1L, function(p) length(unique(p)) == K), ,
                drop = FALSE])
  g
}

# Match fitted state means to the true ones. `perm[i]` is the fitted
# state that plays true state `i`.
hmm_align <- function(mu_hat, mu_true) {
  perms <- hmm_perms(length(mu_true))
  d <- apply(perms, 1L, function(p) max(abs(mu_hat[p] - mu_true)))
  list(perm = perms[which.min(d), ], dist = min(d))
}
