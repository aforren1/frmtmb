# Registry of random-effect covariance structures. Each entry supplies:
#   npar(dim):        number of parameters in `theta` for a block of size dim
#   nll(b, theta, blk): AD log-density of the block's random effects,
#                       vectorized over levels (never loop over levels)
#   vcov(theta, blk): numeric covariance matrix for VarCorr()
#   start(dim):       starting values for theta
#
# Parameterizations follow glmmTMB's covstruct vignette so that fitted
# covariances are directly comparable. `b` for a block is level-major: the
# `dim` coefficients of each level are contiguous (mkReTrms ordering), so
# matrix(b, nrow = dim) has one level per column.

us_chol_cor <- function(theta_cor, d) {
  # Lower-triangular L with unit diagonal, rows normalized: C = Lr Lr'.
  "[<-" <- RTMB::ADoverload("[<-")
  L <- diag(d)
  L[lower.tri(L)] <- theta_cor
  rs <- sqrt((L * L) %*% rep(1, d))
  Lr <- L / as.vector(rs)
  Lr %*% t(Lr)
}

# Unstructured d x d covariance from its theta segment (AD-safe).
us_sigma <- function(theta, d) {
  if (d == 1L) return(RTMB::matrix(exp(2 * theta[1]), 1, 1))
  sdv <- exp(theta[seq_len(d)])
  C <- us_chol_cor(theta[-seq_len(d)], d)
  C * (RTMB::matrix(sdv, ncol = 1) %*% RTMB::matrix(sdv, nrow = 1))
}

covstruct_registry <- list(
  us = list(
    npar = function(dim) dim + dim * (dim - 1L) / 2L,
    sd_idx = function(dim) seq_len(dim),
    nll = function(b, theta, blk) {
      d <- blk$dim
      if (d == 1L) {
        return(sum(RTMB::dnorm(b, 0, exp(theta), log = TRUE)))
      }
      sdv <- exp(theta[seq_len(d)])
      C <- us_chol_cor(theta[-seq_len(d)], d)
      # RTMB::matrix, not base::matrix: base strips the advector class.
      Sigma <- C * (RTMB::matrix(sdv, ncol = 1) %*% RTMB::matrix(sdv, nrow = 1))
      # dim<- (not matrix()) reshapes b: it preserves simref objects, so
      # obj$simulate() can draw the block
      dim(b) <- c(d, length(b) %/% d)
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      d <- blk$dim
      if (d == 1L) {
        V <- matrix(exp(theta)^2, 1, 1)
      } else {
        sdv <- exp(theta[seq_len(d)])
        C <- us_chol_cor(theta[-seq_len(d)], d)
        V <- C * (sdv %o% sdv)
      }
      dimnames(V) <- list(blk$cnms, blk$cnms)
      V
    },
    start = function(dim) numeric(dim + dim * (dim - 1L) / 2L)
  ),
  diag = list(
    npar = function(dim) dim,
    sd_idx = function(dim) seq_len(dim),
    nll = function(b, theta, blk) {
      sdv <- rep(exp(theta), times = blk$n_levels)
      sum(RTMB::dnorm(b, 0, sdv, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- diag(exp(theta)^2, nrow = blk$dim)
      dimnames(V) <- list(blk$cnms, blk$cnms)
      V
    },
    start = function(dim) numeric(dim)
  ),
  homdiag = list(
    npar = function(dim) 1L,
    sd_idx = function(dim) 1L,
    nll = function(b, theta, blk) {
      sum(RTMB::dnorm(b, 0, exp(theta), log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- diag(exp(theta)^2, nrow = blk$dim)
      dimnames(V) <- list(blk$cnms, blk$cnms)
      V
    },
    start = function(dim) 0
  ),
  # Homogeneous AR(1) over the term's factor levels, glmmTMB
  # parameterization: theta = (log sd, phi) with rho = phi/sqrt(1+phi^2).
  ar1 = list(
    npar = function(dim) 2L,
    sd_idx = function(dim) 1L,
    nll = function(b, theta, blk) {
      "[<-" <- RTMB::ADoverload("[<-")
      d <- blk$dim
      sd1 <- exp(theta[1])
      rho <- theta[2] / sqrt(1 + theta[2]^2)
      # rho^|i-j| via sequential products: safe for negative rho on the
      # tape (pow would go through exp/log), and d is small
      pows <- rep(rho, d)
      pows[1] <- 1
      for (k in seq_len(d - 1L) + 1L) pows[k] <- pows[k - 1L] * rho
      M <- abs(outer(seq_len(d), seq_len(d), "-")) + 1L
      C <- RTMB::matrix(pows[as.vector(M)], d, d)
      Sigma <- sd1^2 * C
      dim(b) <- c(d, length(b) %/% d)
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      d <- blk$dim
      rho <- theta[2] / sqrt(1 + theta[2]^2)
      V <- exp(theta[1])^2 * rho^abs(outer(seq_len(d), seq_len(d), "-"))
      dimnames(V) <- list(blk$cnms, blk$cnms)
      V
    },
    start = function(dim) c(0, 0)
  ),
  # Heterogeneous compound symmetry, glmmTMB-style parameter count (d+1):
  # d log-sds plus one correlation mapped onto (-1/(d-1), 1).
  cs = list(
    npar = function(dim) dim + 1L,
    sd_idx = function(dim) seq_len(dim),
    nll = function(b, theta, blk) {
      d <- blk$dim
      sdv <- exp(theta[seq_len(d)])
      a <- 1 / (d - 1)
      rho <- -a + (1 + a) / (1 + exp(-theta[d + 1L]))
      C <- diag(d) * (1 - rho) + rho
      Sigma <- C * (RTMB::matrix(sdv, ncol = 1) %*% RTMB::matrix(sdv, nrow = 1))
      dim(b) <- c(d, length(b) %/% d)
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      d <- blk$dim
      sdv <- exp(theta[seq_len(d)])
      a <- 1 / (d - 1)
      rho <- -a + (1 + a) / (1 + exp(-theta[d + 1L]))
      C <- diag(d) * (1 - rho) + rho
      V <- C * (sdv %o% sdv)
      dimnames(V) <- list(blk$cnms, blk$cnms)
      V
    },
    start = function(dim) numeric(dim + 1L)
  )
)

# Continuous-position AR / exponential covariance over num_factor()
# levels: Sigma_ij = sd^2 * exp(-rate * |t_i - t_j|). theta = (log sd,
# log rate); blk$aux_D holds the distance matrix.
covstruct_registry$ou <- list(
  npar = function(dim) 2L,
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    Sigma <- exp(2 * theta[1]) * exp(-exp(theta[2]) * blk$aux_D)
    dim(b) <- c(blk$dim, length(b) %/% blk$dim)
    sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- exp(theta[1])^2 * exp(-exp(theta[2]) * blk$aux_D)
    dimnames(V) <- list(blk$cnms, blk$cnms)
    V
  },
  start = function(dim) c(0, 0)
)

# Heterogeneous Toeplitz (glmmTMB parameter count 2 * dim - 1): d
# log-sds plus banded correlations rho_k = phi_k / sqrt(1 + phi_k^2).
# PD is not guaranteed for every parameter value (same as glmmTMB); the
# optimizer stays in the feasible region.
covstruct_registry$toep <- list(
  npar = function(dim) 2L * dim - 1L,
  sd_idx = function(dim) seq_len(dim),
  nll = function(b, theta, blk) {
    "[<-" <- RTMB::ADoverload("[<-")
    d <- blk$dim
    sdv <- exp(theta[seq_len(d)])
    phi <- theta[d + seq_len(d - 1L)]
    cvec <- rep(phi[1], d)
    cvec[1] <- 1
    for (k in seq_len(d - 1L)) {
      cvec[k + 1L] <- phi[k] / sqrt(1 + phi[k]^2)
    }
    M <- abs(outer(seq_len(d), seq_len(d), "-")) + 1L
    C <- RTMB::matrix(cvec[as.vector(M)], d, d)
    Sigma <- C * (RTMB::matrix(sdv, ncol = 1) %*% RTMB::matrix(sdv, nrow = 1))
    dim(b) <- c(d, length(b) %/% d)
    sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
  },
  vcov = function(theta, blk) {
    d <- blk$dim
    sdv <- exp(theta[seq_len(d)])
    phi <- theta[d + seq_len(d - 1L)]
    cvec <- c(1, phi / sqrt(1 + phi^2))
    C <- matrix(cvec[abs(outer(seq_len(d), seq_len(d), "-")) + 1L], d, d)
    V <- C * (sdv %o% sdv)
    dimnames(V) <- list(blk$cnms, blk$cnms)
    V
  },
  start = function(dim) numeric(2L * dim - 1L)
)

# Known-covariance random effects: level-major b ~ N(0, A (x) Sigma)
# with A a fixed matrix over the grouping levels (phylogenetic,
# pedigree, neighbor structures) and Sigma an unstructured d x d
# within-level covariance - brms (x | gr(g, cov = A)). Correlation runs
# ACROSS levels, so the block is one multivariate observation of
# dimension d * n_levels; the Kronecker product is assembled from
# precomputed index maps (blk$aux_kron). Dense solve: fine into the
# hundreds of levels, revisit (sparse precision) for thousands.
covstruct_registry$gr_cov <- list(
  npar = function(dim) dim + dim * (dim - 1L) / 2L,
  sd_idx = function(dim) seq_len(dim),
  nll = function(b, theta, blk) {
    if (blk$dim == 1L) {
      Sigma <- exp(2 * theta[1]) * blk$aux_A
      return(sum(RTMB::dmvnorm(b, 0, Sigma, log = TRUE)))
    }
    S <- us_sigma(theta, blk$dim)
    D <- blk$dim * blk$n_levels
    K <- RTMB::matrix(as.vector(blk$aux_A)[blk$aux_kron$ia] *
                        as.vector(S)[blk$aux_kron$is], D, D)
    sum(RTMB::dmvnorm(b, 0, K, log = TRUE))
  },
  vcov = function(theta, blk) {
    d <- blk$dim
    V <- if (d == 1L) {
      matrix(exp(theta[1])^2, 1, 1)
    } else {
      sdv <- exp(theta[seq_len(d)])
      us_chol_cor(theta[-seq_len(d)], d) * (sdv %o% sdv)
    }
    dimnames(V) <- list(blk$cnms, blk$cnms)
    V
  },
  start = function(dim) numeric(dim + dim * (dim - 1L) / 2L)
)

# Smooth wiggly blocks are iid-Gaussian with one variance (the inverse
# smoothing parameter); reuse the homdiag machinery under its own name so
# blocks stay self-describing.
covstruct_registry$smooth <- covstruct_registry$homdiag
