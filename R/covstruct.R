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

# Known-covariance intercepts: b ~ N(0, sd^2 * A) with A a fixed matrix
# over the grouping levels (phylogenetic, pedigree, spatial neighbor
# structures) - brms (1 | gr(g, cov = A)). Correlation is ACROSS levels,
# so the block is one multivariate observation. Dense solve: fine for
# hundreds of levels, revisit (sparse precision) for thousands.
covstruct_registry$gr_cov <- list(
  npar = function(dim) 1L,
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    Sigma <- exp(2 * theta[1]) * blk$aux_A
    sum(RTMB::dmvnorm(b, 0, Sigma, log = TRUE))
  },
  vcov = function(theta, blk) {
    matrix(exp(theta[1])^2, 1, 1,
           dimnames = list(blk$cnms, blk$cnms))
  },
  start = function(dim) 0
)

# Smooth wiggly blocks are iid-Gaussian with one variance (the inverse
# smoothing parameter); reuse the homdiag machinery under its own name so
# blocks stay self-describing.
covstruct_registry$smooth <- covstruct_registry$homdiag
