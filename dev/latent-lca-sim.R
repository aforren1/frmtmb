# Lane `latent`, item 2.4: the simulator for the plan's lca design.
#
# n = 2000, ten binary items, K = 4, two covariates on membership. This
# is the same construction as the Phase 0 scale row
# (extensions/frmtmb.latent/tests/testthat/test-scale.R,
# `latent_lca_data()`), lifted here so that a replicate can be drawn at
# an arbitrary seed. The truth it returns is what item 2.4's recovery
# table is scored against.

# Class 1 is the reference row of the gating matrix, which is poLCA's
# convention; frmtmb's theta dpars reference the LAST class, so the
# comparison re-references before it compares.
lca_truth_gam <- rbind(c(0, 0, 0), c(-0.4, 0.8, -0.5),
                       c(0.2, -0.6, 0.9), c(-0.1, 0.3, 0.4))

lca_truth_base <- function(K = 4L, J = 10L) {
  base <- matrix(0.2, K, J)
  for (k in seq_len(K)) {
    base[k, ((k - 1L) * 2L + 1L):((k - 1L) * 2L + 3L)] <- 0.85
  }
  base
}

lca_sim <- function(seed, n = 2000L, J = 10L, K = 4L) {
  gam <- lca_truth_gam
  base <- lca_truth_base(K, J)
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rbinom(n, 1L, 0.5)
  eta <- cbind(1, x1, x2) %*% t(gam)
  pr <- exp(eta) / rowSums(exp(eta))
  cl <- apply(pr, 1L, function(p) sample.int(K, 1L, prob = p))
  Y <- matrix(0L, n, J)
  for (j in seq_len(J)) Y[, j] <- 1L + stats::rbinom(n, 1L, base[cl, j])
  dd <- data.frame(x1 = x1, x2 = factor(x2))
  dd$Y <- Y
  list(dd = dd, Y = Y, cl = cl, base = base, gam = gam, seed = seed)
}

# Every permutation of seq_len(K), as rows of a matrix. A finite mixture
# is identified only up to a relabeling, so every comparison in this lane
# scores the BEST matching and reports which permutation it used.
lca_perms <- function(K) {
  g <- as.matrix(expand.grid(rep(list(seq_len(K)), K)))
  g <- g[apply(g, 1L, function(p) length(unique(p)) == K), , drop = FALSE]
  unname(g)
}

# Match the rows of `a` to the rows of `b` by minimizing the largest
# absolute cell difference. Returns the permutation `p` with
# `a[p, ] ~ b`.
lca_align <- function(a, b, perms = lca_perms(nrow(a))) {
  d <- apply(perms, 1L, function(p) max(abs(a[p, , drop = FALSE] - b)))
  list(perm = perms[which.min(d), ], dist = min(d))
}
