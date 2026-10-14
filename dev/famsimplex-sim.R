# Simulators for the simplex-response families, shared by the brms fit
# script, the recovery study and the agreement check. Not a package file.

# softmax over the columns of an n x K matrix of predictors
fs_softmax <- function(E) {
  E <- E - apply(E, 1L, max)
  P <- exp(E)
  P / rowSums(P)
}

# Dirichlet draws by gamma normalization: alpha is n x K
fs_rdirichlet <- function(alpha) {
  G <- matrix(stats::rgamma(length(alpha), shape = alpha), nrow(alpha))
  G / rowSums(G)
}

# Logistic-normal draws: an MVN draw on the additive log-ratio scale
# against the reference column `ref`, mapped back by the softmax.
# mu is n x (K - 1), Sigma is (K - 1) x (K - 1).
fs_rlogistic_normal <- function(mu, Sigma, ref = 1L) {
  n <- nrow(mu)
  L <- t(chol(Sigma))
  Z <- matrix(stats::rnorm(n * ncol(mu)), n) %*% t(L)
  X <- mu + Z
  K <- ncol(mu) + 1L
  E <- matrix(0, n, K)
  E[, -ref] <- X
  fs_softmax(E)
}

# The Dirichlet data set: K = 3 categories a, b, c, reference a.
# mub = 0.4 + 0.6 x (+ u_b[g]), muc = -0.3 - 0.5 x (+ u_c[g]),
# log phi = 2 + 0.5 z.
fs_data_dirichlet <- function(seed, n = 300, ng = 0, sd_re = c(0.5, 0.3),
                              cor_re = 0.4, phi_formula = TRUE) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  E <- cbind(0, 0.4 + 0.6 * dd$x, -0.3 - 0.5 * dd$x)
  if (ng > 0) {
    dd$g <- factor(rep(seq_len(ng), length.out = n))
    S <- diag(sd_re) %*% matrix(c(1, cor_re, cor_re, 1), 2) %*% diag(sd_re)
    U <- matrix(stats::rnorm(ng * 2), ng) %*% chol(S)
    E[, 2:3] <- E[, 2:3] + U[as.integer(dd$g), ]
  }
  phi <- exp(2 + if (phi_formula) 0.5 * dd$z else 0)
  Y <- fs_rdirichlet(fs_softmax(E) * phi)
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd
}

# The logistic-normal data set: K = 3 categories a, b, c, reference a,
# mub = 0.4 + 0.6 x, muc = -0.3 - 0.5 x (+ group effects),
# log sigmab = -0.5 + 0.4 z, sigmac = exp(-0.2), lncor = 0.3.
fs_data_logistic_normal <- function(seed, n = 300, ng = 0,
                                    sd_re = c(0.5, 0.3), cor_re = 0.4,
                                    sigma_formula = TRUE, lncor = 0.3) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  M <- cbind(0.4 + 0.6 * dd$x, -0.3 - 0.5 * dd$x)
  if (ng > 0) {
    dd$g <- factor(rep(seq_len(ng), length.out = n))
    S <- diag(sd_re) %*% matrix(c(1, cor_re, cor_re, 1), 2) %*% diag(sd_re)
    U <- matrix(stats::rnorm(ng * 2), ng) %*% chol(S)
    M <- M + U[as.integer(dd$g), ]
  }
  sb <- exp(-0.5 + if (sigma_formula) 0.4 * dd$z else numeric(n))
  sc <- rep(exp(-0.2), n)
  Y <- matrix(0, n, 3)
  R <- matrix(c(1, lncor, lncor, 1), 2)
  for (i in seq_len(n)) {
    Si <- diag(c(sb[i], sc[i])) %*% R %*% diag(c(sb[i], sc[i]))
    Y[i, ] <- fs_rlogistic_normal(M[i, , drop = FALSE], Si)
  }
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd
}

# The dirichlet2 data set: alpha_k = exp(eta_k), all three predicted.
fs_data_dirichlet2 <- function(seed, n = 300) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  A <- exp(cbind(1 + 0.3 * dd$x, 1.5 - 0.4 * dd$x, 0.8 + 0.2 * dd$x))
  Y <- fs_rdirichlet(A)
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd
}
