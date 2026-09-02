# The rung-1 recipe: a K = 2 gaussian HMM expressed as a frmtmb
# custom_family(). Sourced by probeC / probeD / probeF.
#
# Contract used:
#   response      y
#   vint(g, t)    sequence id and within-sequence order
#   dpars         mu, mu2      state means      (identity)
#                 sigma1, sigma2 state SDs      (log)
#                 tr12, tr22   transition logits (identity), one per
#                              row of the transition matrix, state 1
#                              the reference cell
#   extra_pars    hmm_ldel     initial-distribution logit
#
# The lpdf runs the forward recursion once per sequence and puts each
# sequence's log likelihood on its FIRST row, zero elsewhere, because
# the objective only ever forms sum(weights * lpdf).

hmm2_lpdf <- function(y, dpars, aterms, extra = list(hmm_ldel = 0)) {
  "c" <- RTMB::ADoverload("c")
  K <- 2L
  n <- length(y)
  rows_by_g <- hmm_seq_index(aterms$vint1, aterms$vint2)

  lp <- list(RTMB::dnorm(y, dpars$mu, dpars$sigma1, log = TRUE),
             RTMB::dnorm(y, dpars$mu2, dpars$sigma2, log = TRUE))
  eta <- list(list(0, dpars$tr12), list(0, dpars$tr22))
  lg <- tpm_logs_ad(eta, K)
  d <- softmax0_ad(extra$hmm_ldel, K)
  ld <- log(d)

  llv <- NULL
  for (gi in seq_along(rows_by_g)) {
    v <- fwd_ad_log_tv(lp, lg, rows_by_g[[gi]], ld, K)
    llv <- if (is.null(llv)) v else c(llv, v)
  }
  # scatter to the first row of each sequence with a constant sparse
  # matrix: no elementwise assignment on the tape at all
  first <- vapply(rows_by_g, function(r) r[1], integer(1))
  S <- Matrix::sparseMatrix(i = first, j = seq_along(first), x = 1,
                            dims = c(n, length(first)))
  as.vector(S %*% llv)
}

# Stationary-initial-distribution variant. The transition matrix must be
# constant (row 1's values are read off), and delta solves
# delta (I - Gamma + 1 1') = 1' on the tape. RTMB::solve, not base's:
# the S4 advector method is not imported by a probe script.
hmm2_lpdf_stat <- function(y, dpars, aterms, extra = list()) {
  "c" <- RTMB::ADoverload("c")
  K <- 2L
  n <- length(y)
  rows_by_g <- hmm_seq_index(aterms$vint1, aterms$vint2)

  lp <- list(RTMB::dnorm(y, dpars$mu, dpars$sigma1, log = TRUE),
             RTMB::dnorm(y, dpars$mu2, dpars$sigma2, log = TRUE))
  eta <- list(list(0, dpars$tr12), list(0, dpars$tr22))
  lg <- tpm_logs_ad(eta, K)

  # constant transitions: read row 1 and build the K x K matrix
  g11 <- exp(lg[[1]][[1]][1]); g12 <- exp(lg[[1]][[2]][1])
  g21 <- exp(lg[[2]][[1]][1]); g22 <- exp(lg[[2]][[2]][1])
  G <- RTMB::matrix(c(g11, g21, g12, g22), K, K)   # column-major
  A <- RTMB::matrix(c(1, 1, 1, 1), K, K) + diag(K) - G
  d <- as.vector(RTMB::solve(t(A), c(1, 1)))
  ld <- log(d)

  llv <- NULL
  for (gi in seq_along(rows_by_g)) {
    v <- fwd_ad_log_tv(lp, lg, rows_by_g[[gi]], ld, K)
    llv <- if (is.null(llv)) v else c(llv, v)
  }
  first <- vapply(rows_by_g, function(r) r[1], integer(1))
  S <- Matrix::sparseMatrix(i = first, j = seq_along(first), x = 1,
                            dims = c(n, length(first)))
  as.vector(S %*% llv)
}

hmm2_family_stat <- function() {
  frmtmb::custom_family(
    "hmm2_gaussian_stat",
    dpars = c("mu", "mu2", "sigma1", "sigma2", "tr12", "tr22"),
    links = list(mu = "identity", mu2 = "identity",
                 sigma1 = "log", sigma2 = "log",
                 tr12 = "identity", tr22 = "identity"),
    lpdf = hmm2_lpdf_stat,
    type = "continuous",
    init_dpars = list(
      mu = function(y, aterms) unname(stats::quantile(y, 0.25)),
      mu2 = function(y, aterms) unname(stats::quantile(y, 0.75)),
      sigma1 = function(y, aterms) stats::sd(y),
      sigma2 = function(y, aterms) stats::sd(y),
      tr12 = function(y, aterms) -1.5,
      tr22 = function(y, aterms) 1.5
    )
  )
}

hmm2_family <- function(post = list(), sim = NULL) {
  frmtmb::custom_family(
    "hmm2_gaussian",
    dpars = c("mu", "mu2", "sigma1", "sigma2", "tr12", "tr22"),
    links = list(mu = "identity", mu2 = "identity",
                 sigma1 = "log", sigma2 = "log",
                 tr12 = "identity", tr22 = "identity"),
    lpdf = hmm2_lpdf,
    type = "continuous",
    init_dpars = list(
      # mixture()'s quantile convention: spread the state means over the
      # response, start both SDs at the pooled SD, transitions sticky
      mu = function(y, aterms) unname(stats::quantile(y, 0.25)),
      mu2 = function(y, aterms) unname(stats::quantile(y, 0.75)),
      sigma1 = function(y, aterms) stats::sd(y),
      sigma2 = function(y, aterms) stats::sd(y),
      tr12 = function(y, aterms) -1.5,
      tr22 = function(y, aterms) 1.5
    ),
    extra_pars = function(y, aterms) list(hmm_ldel = 0),
    post = post, sim = sim
  )
}
