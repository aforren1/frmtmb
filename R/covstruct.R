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

# Heterogeneous AR(1): d log-sds plus one phi (rho = phi/sqrt(1+phi^2)).
covstruct_registry$hetar1 <- list(
  npar = function(dim) dim + 1L,
  sd_idx = function(dim) seq_len(dim),
  nll = function(b, theta, blk) {
    "[<-" <- RTMB::ADoverload("[<-")
    d <- blk$dim
    sdv <- exp(theta[seq_len(d)])
    rho <- theta[d + 1L] / sqrt(1 + theta[d + 1L]^2)
    pows <- rep(rho, d)
    pows[1] <- 1
    for (k in seq_len(d - 1L) + 1L) pows[k] <- pows[k - 1L] * rho
    M <- abs(outer(seq_len(d), seq_len(d), "-")) + 1L
    C <- RTMB::matrix(pows[as.vector(M)], d, d)
    Sigma <- C * (RTMB::matrix(sdv, ncol = 1) %*% RTMB::matrix(sdv, nrow = 1))
    dim(b) <- c(d, length(b) %/% d)
    sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
  },
  vcov = function(theta, blk) {
    d <- blk$dim
    sdv <- exp(theta[seq_len(d)])
    rho <- theta[d + 1L] / sqrt(1 + theta[d + 1L]^2)
    V <- rho^abs(outer(seq_len(d), seq_len(d), "-")) * (sdv %o% sdv)
    dimnames(V) <- list(blk$cnms, blk$cnms)
    V
  },
  start = function(dim) numeric(dim + 1L)
)

# Homogeneous compound symmetry: one log-sd plus one correlation on
# (-1/(d-1), 1).
covstruct_registry$homcs <- list(
  npar = function(dim) 2L,
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    d <- blk$dim
    a <- 1 / (d - 1)
    rho <- -a + (1 + a) / (1 + exp(-theta[2]))
    C <- diag(d) * (1 - rho) + rho
    Sigma <- exp(2 * theta[1]) * C
    dim(b) <- c(d, length(b) %/% d)
    sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
  },
  vcov = function(theta, blk) {
    d <- blk$dim
    a <- 1 / (d - 1)
    rho <- -a + (1 + a) / (1 + exp(-theta[2]))
    V <- exp(theta[1])^2 * (diag(d) * (1 - rho) + rho)
    dimnames(V) <- list(blk$cnms, blk$cnms)
    V
  },
  start = function(dim) c(0, 0)
)

# Homogeneous Toeplitz: one log-sd plus d-1 banded correlations.
covstruct_registry$homtoep <- list(
  npar = function(dim) dim,
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    "[<-" <- RTMB::ADoverload("[<-")
    d <- blk$dim
    phi <- theta[1L + seq_len(d - 1L)]
    cvec <- rep(phi[1], d)
    cvec[1] <- 1
    for (k in seq_len(d - 1L)) {
      cvec[k + 1L] <- phi[k] / sqrt(1 + phi[k]^2)
    }
    M <- abs(outer(seq_len(d), seq_len(d), "-")) + 1L
    C <- RTMB::matrix(cvec[as.vector(M)], d, d)
    Sigma <- exp(2 * theta[1]) * C
    dim(b) <- c(d, length(b) %/% d)
    sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
  },
  vcov = function(theta, blk) {
    d <- blk$dim
    phi <- theta[1L + seq_len(d - 1L)]
    cvec <- c(1, phi / sqrt(1 + phi^2))
    C <- matrix(cvec[abs(outer(seq_len(d), seq_len(d), "-")) + 1L], d, d)
    V <- exp(theta[1])^2 * C
    dimnames(V) <- list(blk$cnms, blk$cnms)
    V
  },
  start = function(dim) numeric(dim)
)

# Spatial structures over num_factor() coordinates (blk$aux_D holds the
# distance matrix), glmmTMB parameterizations: theta = (log sd,
# log range[, log shape]).
spatial_entry <- function(corr_fn, npar_k) {
  list(
    npar = function(dim) npar_k,
    sd_idx = function(dim) 1L,
    nll = function(b, theta, blk) {
      Sigma <- exp(2 * theta[1]) * corr_fn(blk$aux_D, theta)
      dim(b) <- c(blk$dim, length(b) %/% blk$dim)
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- exp(theta[1])^2 * corr_fn(blk$aux_D, theta)
      V <- as.matrix(V)
      dimnames(V) <- list(blk$cnms, blk$cnms)
      V
    },
    start = function(dim) numeric(npar_k)
  )
}

covstruct_registry$exp <- spatial_entry(
  function(D, theta) exp(-D / exp(theta[2])), 2L
)

covstruct_registry$gau <- spatial_entry(
  function(D, theta) exp(-(D / exp(theta[2]))^2), 2L
)

# Matern correlation 2^(1-nu)/Gamma(nu) (d/range)^nu K_nu(d/range); the
# zero-distance diagonal is masked in (data-only mask, branch-free).
# Bounded internal transforms - smoothness in (0.1, 5), range floored at
# 5% of the median distance - keep every intermediate representable, so
# the optimizer cannot step into 0 * Inf territory (which is where both
# we and glmmTMB otherwise die). Compare fits at the covariance level,
# not the theta level.
covstruct_registry$mat <- spatial_entry(
  function(D, theta) {
    d <- nrow(D)
    mask <- as.numeric(D > 0)
    Dp <- as.vector(D) + (1 - mask)   # zeros -> 1, masked out below
    rng <- 0.05 * stats::median(D[D > 0]) + exp(theta[2])
    u <- Dp / rng
    nu <- 0.1 + 4.9 / (1 + exp(-theta[3]))
    cr <- 2^(1 - nu) / exp(lgamma(nu)) * u^nu * RTMB::besselK(u, nu)
    RTMB::matrix(mask * cr + (1 - mask), d, d)
  }, 3L
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

# Known-precision intercepts: b ~ N(0, sd^2 * Q^-1) with Q a (sparse)
# precision matrix over the levels - gr(g, prec = Q). The sparse GMRF
# density keeps large structures (big phylogenies, spatial graphs)
# tractable where the dense gr(cov=) path is not.
covstruct_registry$gr_prec <- list(
  npar = function(dim) 1L,
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    Qs <- exp(-2 * theta[1]) * blk$aux_Q
    sum(RTMB::dgmrf(b, 0, Qs, log = TRUE))
  },
  vcov = function(theta, blk) {
    matrix(exp(theta[1])^2, 1, 1,
           dimnames = list(blk$cnms, blk$cnms))
  },
  start = function(dim) 0
)

# Exact Gaussian process over observed positions (gp(x) without k=):
# squared-exponential kernel sd^2 exp(-d^2 / (2 rho^2)), matching the
# Hilbert-space approximation's kernel so gp(x) and gp(x, k=) estimate
# the same quantity. theta = (log sd, log rho); dense over the unique
# positions (blk$aux_D).
# the 1e-6 nugget keeps the notoriously ill-conditioned SE kernel
# Cholesky-factorizable as the range grows (standard GP practice)
covstruct_registry$gp <- spatial_entry(
  function(D, theta) {
    exp(-(D / exp(theta[2]))^2 / 2) + diag(1e-6, nrow(D))
  }, 2L
)

# Hilbert-space GP approximation (Riutort-Mayol et al.): sine basis in
# Z, independent coefficients whose prior SDs follow the SE-kernel
# spectral density at the basis frequencies (blk$aux_omega).
hsgp_sds <- function(theta, omega) {
  sd_ <- exp(theta[1])
  rho <- exp(theta[2])
  sqrt(sd_^2 * rho * sqrt(2 * pi) * exp(-0.5 * (rho * omega)^2))
}

covstruct_registry$hsgp <- list(
  npar = function(dim) 2L,
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    sum(RTMB::dnorm(b, 0, hsgp_sds(theta, blk$aux_omega), log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- diag(as.numeric(hsgp_sds(theta, blk$aux_omega))^2,
              nrow = blk$dim)
    dimnames(V) <- list(blk$cnms, blk$cnms)
    V
  },
  start = function(dim) c(0, 0)
)

# Known, fully fixed within-level covariance (glmmTMB equalto): b ~
# N(0, V) with V supplied by the user - zero parameters (meta-analysis
# sampling covariances and similar).
covstruct_registry$equalto <- list(
  npar = function(dim) 0L,
  sd_idx = function(dim) integer(0),
  nll = function(b, theta, blk) {
    dim(b) <- c(blk$dim, length(b) %/% blk$dim)
    sum(RTMB::dmvnorm(t(b), 0, blk$aux_A, log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- blk$aux_A
    dimnames(V) <- list(blk$cnms, blk$cnms)
    V
  },
  start = function(dim) numeric(0)
)

# Reduced-rank / factor-analytic (glmmTMB rr): per-level coefficients
# c = Lambda f with f iid standard normal factors (the block's b segment)
# and Lambda the D x rank loadings matrix, columnwise lower-triangular
# (theta). Covariance = Lambda Lambda', rank-deficient by design. The
# expansion from factors to coefficients happens in expand_b(); npar and
# start are handled at the frame call sites (they need the rank).
rr_npar <- function(dim, rank) {
  as.integer(dim * rank - rank * (rank - 1L) / 2L)
}

# start at Lambda = [I_rank; 0]: identified, unit factor scale
rr_start <- function(dim, rank) {
  th <- numeric(rr_npar(dim, rank))
  pos <- 0L
  for (j in seq_len(rank)) {
    th[pos + 1L] <- 1
    pos <- pos + (dim - j + 1L)
  }
  th
}

rr_loadings <- function(theta, dim, rank) {
  "[<-" <- RTMB::ADoverload("[<-")
  L <- RTMB::matrix(0, dim, rank)
  pos <- 0L
  for (j in seq_len(rank)) {
    len <- dim - j + 1L
    L[seq.int(j, dim), j] <- theta[pos + seq_len(len)]
    pos <- pos + len
  }
  L
}

covstruct_registry$rr <- list(
  npar = function(dim) {
    stop("rr npar needs the rank; handled at the frame call site",
         call. = FALSE)
  },
  sd_idx = function(dim) integer(0),
  nll = function(b, theta, blk) {
    sum(RTMB::dnorm(b, 0, 1, log = TRUE))
  },
  vcov = function(theta, blk) {
    L <- rr_loadings(theta, blk$dim, blk$rank)
    V <- as.matrix(L %*% t(L))
    dimnames(V) <- list(blk$cnms, blk$cnms)
    V
  },
  start = function(dim) {
    stop("rr start needs the rank; handled at the frame call site",
         call. = FALSE)
  }
)

# Coefficient-space vector the Z matrices multiply: identical to b
# except for rr blocks, whose factors expand through the loadings.
# AD-safe (RTMB::matrix + [<- overload) and numeric-safe.
expand_b <- function(frame, b, theta) {
  if (!isTRUE(frame$has_rr)) return(b)
  "[<-" <- RTMB::ADoverload("[<-")
  cvec <- rep(b[1] * 0, frame$n_c)   # keeps the advector class if taped
  for (bk in frame$re_blocks) {
    if (bk$covstruct == "rr") {
      L <- rr_loadings(theta[bk$theta_idx], bk$dim, bk$rank)
      Fm <- RTMB::matrix(b[bk$b_idx], bk$rank, bk$n_levels)
      cvec[bk$c_idx] <- as.vector(L %*% Fm)
    } else {
      cvec[bk$c_idx] <- b[bk$b_idx]
    }
  }
  cvec
}

# Smooth wiggly blocks are iid-Gaussian with one variance (the inverse
# smoothing parameter); reuse the homdiag machinery under its own name so
# blocks stay self-describing.
covstruct_registry$smooth <- covstruct_registry$homdiag
