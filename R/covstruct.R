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

#' Lower-triangular Cholesky factor of the correlation matrix behind the
#' correlation segment of a `us` theta: `C = Lr Lr'`.
#'
#' The unit-diagonal lower-triangular `L` gets row-normalized, so every
#' value of `theta_cor` maps onto a valid correlation matrix. AD-safe.
#'
#' The factor, and not just the product, is what the Student-t blocks
#' need: their density is a quadratic form and a log-determinant, both of
#' which come off a triangular factor without ever forming an inverse.
#'
#' @noRd
us_chol_L <- function(theta_cor, d) {
  "[<-" <- RTMB::ADoverload("[<-")
  L <- diag(d)
  L[lower.tri(L)] <- theta_cor
  rs <- sqrt((L * L) %*% rep(1, d))
  L / as.vector(rs)
}

#' Correlation matrix behind the correlation segment of a `us` theta.
#'
#' @noRd
us_chol_cor <- function(theta_cor, d) {
  Lr <- us_chol_L(theta_cor, d)
  Lr %*% t(Lr)
}

#' Inverse of us_chol_cor(): the theta segment whose row-normalized unit
#' lower-triangular L reproduces the correlation matrix C. C = Lc Lc'
#' with Lc lower-triangular has unit-norm rows (diag(C) = 1), so Lc is
#' already the normalized Lr; dividing row i by `Lc[i, i]` undoes the
#' normalization and recovers L. Numeric only - this runs off the tape.
#'
#' @noRd
us_theta_cor <- function(C) {
  d <- nrow(C)
  Lc <- tryCatch(t(chol(C)), error = function(e) {
    stop("The requested correlation matrix is not positive definite",
         call. = FALSE)
  })
  # column-major recycling divides row i by diag(Lc)[i]
  L <- Lc / diag(Lc)
  L[lower.tri(L)]
}

#' Unstructured d x d covariance from its theta segment (AD-safe).
#'
#' @noRd
us_sigma <- function(theta, d) {
  if (d == 1L) return(RTMB::matrix(exp(2 * theta[1]), 1, 1))
  sdv <- exp(theta[seq_len(d)])
  C <- us_chol_cor(theta[-seq_len(d)], d)
  C * (RTMB::matrix(sdv, ncol = 1) %*% RTMB::matrix(sdv, nrow = 1))
}

# A structure that carries `from_natural(sds, C, blk)` can be set from
# natural-scale standard deviations (length blk$dim) and a correlation
# matrix (dim x dim, or NULL when the structure has no correlation
# parameters); it returns the WHOLE theta segment, so no entry of the
# block is left at a stale value. frm_simulate()'s natural-scale
# newparams path refuses structures without it.
#' The one log standard deviation a homogeneous block keeps in `theta`.
#' Natural-scale standard deviations that differ have no representation in
#' such a block, so the call stops instead.
#'
#' @noRd
homogeneous_sd <- function(sds, what) {
  if (length(unique(signif(sds, 12))) > 1L) {
    stop("A '", what, "' block has one shared standard deviation; got ",
         paste(signif(sds, 4), collapse = ", "), call. = FALSE)
  }
  log(sds[1L])
}

# ------------------------------------------------- Student-t latents
#
# brms `gr(g, dist = "student")` builds a level's effects as
#
#   b_j = sqrt(nu * u_j) * W z_j,   u_j ~ inv-chi2(nu),  z_j ~ N(0, I),
#
# with `W` the Cholesky factor of the scale matrix `Sigma = W W'` and ONE
# mixing variable `u_j` per LEVEL, shared across that level's
# coefficients whether or not the block is correlated. That is exactly
# the standard multivariate t with `nu` degrees of freedom and scale
# matrix `Sigma` (verified against the mixture integral in
# dev/tre/probeF1-mvt-and-limits.R), so the marginal density is closed
# form and TMB can integrate it the same way it integrates a gaussian
# one.
#
# `theta` therefore holds the SCALE, not the standard deviation. The two
# differ by sqrt(nu/(nu-2)); see student_var_factor() and ?gr.

#' Student-t distributed random effects
#'
#' `(x | gr(g, dist = "student"))` gives a grouping term a Student-t
#' latent instead of a gaussian one, which is brms's spelling
#' (`brms::gr()`, argument `dist`). A group far from the others then
#' costs the variance component much less than it does under a gaussian
#' latent, because the t's tail can hold it.
#'
#' @section Parameterization:
#'
#' A level's coefficients are drawn as
#' `b_j = sqrt(nu * u_j) W z_j` with `u_j ~ inv-chi2(nu)` and
#' `z_j ~ N(0, I)`, which is the multivariate t with `nu` degrees of
#' freedom and scale matrix `Sigma = W W'`. The mixing variable `u_j` is
#' per LEVEL and shared across that level's coefficients, so a
#' correlated or a `diag()` block is one multivariate t and not several
#' independent univariate ones. brms builds it the same way.
#'
#' **`Sigma` is the SCALE, not the covariance.** The variance is
#' `Sigma * nu / (nu - 2)`, so a standard deviation is
#' `scale * sqrt(nu / (nu - 2))`. `VarCorr()` stores the scale matrix,
#' tags it with `nu`, and prints both columns; `confint()`, `variables()`
#' and `frm_simulate(newparams = )` speak of it as `sd_<group>__<term>`,
#' which is the name brms gives the same quantity. The correlations are
#' the same either way.
#'
#' @section Why `nu` is fixed:
#'
#' brms estimates `nu` under a `gamma(2, 0.1)` prior truncated at 1, and
#' that prior is carrying the parameter. Maximum likelihood has no such
#' help: profiling a simulated fit with 20 groups leaves the whole grid
#' from 2.1 to 500 inside the 95% profile interval, and joint ML sends
#' `nu` to a boundary in 24% to 41% of replicates at 20 groups and still
#' 5% to 14% at 100 (`dev/tre-feasibility.md`). So `nu` is a constant
#' here, set by `dist_nu` and defaulting to 5. That is the frequentist
#' analogue of brms's own `prior(constant(3), class = "df")`. To ask
#' what the data say about it, fit two or three values and compare
#' `logLik()`.
#'
#' `dist_nu` must exceed 2, so that the latent has a variance for
#' `VarCorr()` to convert and for new-level prediction to use.
#'
#' @section Accuracy:
#'
#' The t density is not log-concave, so the Laplace approximation is not
#' exact over it as it is over a gaussian latent. Measured against
#' adaptive quadrature (`dev/tre-feasibility.md`), the approximation
#' pushes the estimated latent SCALE UPWARD, and how much depends on how
#' much the data say about each level:
#'
#' \itemize{
#'   \item With the latent scale near the residual SD, the bias is
#'     under 2% of one standard error at 8 observations per group and
#'     0.2% at 25, at every `nu` from 2.5 up. It reaches 12% of a
#'     standard error in the worst case tested, `nu = 2.5` with 3
#'     observations per group.
#'   \item It becomes material only where the variance component is
#'     small AND the groups are tiny: at 2 observations per group and a
#'     true scale a quarter of the residual SD, the Laplace estimate of
#'     the scale came back three times the exact one.
#' }
#'
#' Two checks, both already in the package. `quadrature = TRUE`
#' marginalizes a scalar random intercept by Gauss-Kronrod quadrature
#' instead, which over a t latent is EXACT, not merely better; it is the
#' one to run when a t block's variance component matters and the groups
#' are small. `frmtmb.sample::check_laplace()` measures the same thing without
#' refitting, by NUTS on the objective: its `z_shift` for
#' `theta` reproduced the displacement above to within a percentage
#' point in the probe.
#'
#' A last consequence of the constant: `logLik()`, `AIC()` and `BIC()`
#' carry roughly `G * c(nu)` where `G` is the number of levels and
#' `c(nu) = lgamma((nu+1)/2) - lgamma(nu/2) + log(2/(nu+1))/2` (-0.226 at
#' `nu = 3`, -0.141 at `nu = 5`). Comparing two t fits with the same
#' `nu` and the same grouping is fine, because the offset cancels;
#' comparing a t fit against a gaussian one by AIC is not.
#'
#' @section What is refused:
#'
#' \describe{
#'   \item{`gr(cov = )` / `gr(prec = )`}{A relationship matrix
#'     correlates the LEVELS, and the t's mixing variable is per level,
#'     so the joint density over the field is not a multivariate t and
#'     has no closed form to hand the Laplace machinery. brms writes
#'     this combination, but as a hierarchical construction Stan samples
#'     rather than a density.}
#'   \item{Other covariance structures}{`us` (the default) and `diag`
#'     only. `ar1()`, `cs()`, `toep()` and the rest describe a
#'     covariance over the block's LEVELS, which the per-level mixing
#'     variable does not compose with.}
#'   \item{`mm()`}{A multi-membership row loads several levels at once,
#'     so the per-level mixing variable has no single value on it.}
#'   \item{`|ID|` keys}{Merged blocks are assembled as gaussian ones.
#'     Write the merged coefficients as one term instead -
#'     `(x1 + x2 | gr(g, dist = "student"))` is the same multivariate-t
#'     block.}
#' }
#'
#' @section Downstream:
#'
#' `ranef()` and its conditional variances, `sdreport()` standard
#' errors, `REML = TRUE` and `frm_bootstrap()` all work: the Laplace
#' machinery does not care which density the latent has. `simulate()`
#' and `frm_simulate()` draw a multivariate t with one chi-square per
#' level. `predict(allow_new_levels = TRUE)` inflates the unseen level's
#' variance by `nu / (nu - 2)`, so the interval has the right variance
#' around a heavier-tailed truth. It is still built as a gaussian
#' interval, so it is not the right quantile far into the tail. Use
#' `simulate()` for that.
#'
#' @return `gr(g, dist = "student")` is a formula term, not a
#'   free-standing function: `bf()` reads it at parse time, and the
#'   value it contributes is the heavy-tailed random-effect block of
#'   the model, reachable through [ranef()] and [VarCorr()]. This page
#'   itself documents the term grammar and returns nothing.
#' @name frmtmb-student-re
#' @seealso [VarCorr()] for the scale matrix, [frm_compat()] for what a
#'   `dist = "student"` block may be combined with, and
#'   `vignette("brms-migration")`.
#' @examples
#' set.seed(1)
#' n <- 12
#' d <- data.frame(x = rnorm(20 * n), g = factor(rep(1:20, each = n)))
#' b <- rnorm(20)
#' b[20] <- b[20] + 6          # one outlying group
#' d$y <- 1 + 0.5 * d$x + b[d$g] + rnorm(20 * n)
#'
#' fit_t <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))),
#'              family = gaussian(), data = d)
#' fit_n <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
#'
#' # the gaussian latent has to widen to cover the outlying group
#' VarCorr(fit_t)
#' VarCorr(fit_n)
#'
#' # heavier tails, at the cost of a fixed nu
#' frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 3))),
#'     family = gaussian(), data = d)
NULL

#' The degrees-of-freedom default and floor for `gr(dist = "student")`.
#'
#' The floor is 2 and not brms's 1 because everything downstream that
#' reports a variance needs one to exist: at `nu <= 2` the latent has no
#' finite variance, so `VarCorr()`'s documented conversion, the
#' new-level prediction variance and the simulation of a marginal
#' response all lose their meaning. The default is settled in
#' dev/tre-feasibility.md section 5 on the probe D2 sweep.
#'
#' @noRd
student_nu_default <- 5
student_nu_floor <- 2

#' Variance inflation from the t's scale to its variance.
#'
#' @noRd
student_var_factor <- function(nu) nu / (nu - 2)

#' Whether a random-effect block carries a Student-t latent.
#'
#' @noRd
is_student_block <- function(bk) !is.null(bk[["dist_nu"]])

#' Log-density of `n` scalar Student-t latents with log scale
#' `log_scale`, vectorized. The `d = 1` case of student_lpdf_core(),
#' written out because a scalar block needs no factor and no solve.
#'
#' @noRd
student_lpdf_scalar <- function(b, log_scale, nu) {
  z <- b * exp(-log_scale)
  student_lpdf_core(z * z, log_scale, nu, 1L, length(b))
}

#' Multivariate-t log-density from the per-level Mahalanobis forms `q`
#' and the log-determinant of the scale matrix's Cholesky factor.
#'
#' `log1p`, not `log(1 + .)`: the tail term is multiplied by
#' `(nu + d)/2`, so at a large `nu` the rounding of `1 + q/nu` is
#' amplified by the same factor and the density comes back wrong in the
#' first decimal. A large `nu` is how the gaussian limit is reached, and
#' how a user checks that a t block reduces to one.
#'
#' @noRd
student_lpdf_core <- function(q, ldet_W, nu, d, n_levels) {
  n_levels * (lgamma((nu + d) / 2) - lgamma(nu / 2) -
                d / 2 * log(nu * pi) - ldet_W) -
    (nu + d) / 2 * sum(log1p(q / nu))
}

#' Log-density of a level-major `b` under a multivariate t with `nu`
#' degrees of freedom and scale matrix `W W'`.
#'
#' `W` is lower triangular, so the quadratic form comes from one
#' triangular solve and the log-determinant from its diagonal: no
#' inverse and no second factorization. `dim<-` (not `matrix()`)
#' reshapes `b`, which keeps simref objects intact.
#'
#' @noRd
student_lpdf <- function(b, W, nu, d, n_levels) {
  dim(b) <- c(d, n_levels)
  z <- RTMB::solve(W, b)
  student_lpdf_core(RTMB::colSums(z * z), sum(log(RTMB::diag(W))),
                    nu, d, n_levels)
}

# ------------------------------------------------- AR(1), at O(d) ---
#
# THE SCALE QUESTION, ANSWERED ON PAPER, because getting it wrong is a
# silently wrong likelihood.
#
# `RTMB::dautoreg(x, phi = rho)` is the density of a STATIONARY AR(1)
# with UNIT MARGINAL variance: it draws `x[1] ~ N(0, 1)` from the
# stationary law and then `x[i] ~ N(rho x[i-1], sqrt(1 - rho^2))`. The
# innovation standard deviation is NOT 1; it is the value that keeps the
# marginal variance at 1. frmtmb's `ar1`/`hetar1` theta carries MARGINAL
# standard deviations, which is the same convention, so the two differ
# only by that scale and by nothing about the recursion.
#
# The relation. Write `C` for the AR(1) correlation matrix
# `C_ij = rho^|i-j|` and `D = diag(sd_1, ..., sd_d)` for the marginal
# standard deviations (all equal for `ar1`, free for `hetar1`). The
# block's covariance is `Sigma = D C D`, so with `z = D^-1 x`,
#
#   log N(x; 0, D C D) = log N(z; 0, C) - sum_i log(sd_i)
#
# and `log N(z; 0, C)` is what `dautoreg` computes. The `- sum log sd`
# is the Jacobian of the standardization, one term per TIME point (not
# per level), which is why `hetar1`'s per-time sds enter here and
# nowhere else. Equivalently: `dautoreg(x, phi = rho, scale = sd)`,
# whose `scale` argument does exactly this elementwise and accepts the
# length-d vector `hetar1` needs.
#
# The correlation part in closed form. `C^-1` is tridiagonal, so
#
#   z' C^-1 z = z_1^2 + sum_{i=2}^d (z_i - rho z_{i-1})^2 / (1 - rho^2)
#   log|C|    = (d - 1) log(1 - rho^2)
#
# which is the innovation form `dautoreg` recurses over, written out.
# Cost is O(d) per level, and every term is a whole-vector operation, so
# the tape sees no loop at all.
#
# One more simplification the parameterization hands us: theta's
# correlation is `rho = t / sqrt(1 + t^2)`, so
#
#   1 - rho^2 = 1 / (1 + t^2)   exactly.
#
# The reciprocal `1 + t^2` is therefore formed by addition rather than
# by the cancelling subtraction `1 - rho^2`, which is where a rho near
# +/-1 would otherwise lose its significant digits.
#
# WHY NOT `dautoreg` ITSELF, and why not `dgmrf`. All three routes agree
# to 7e-15, and the gates in test-sparsear1.R hold the other two against
# this one. `dautoreg` loops over time points in R and takes ONE vector,
# so it has to be called once per level: measured at 0.28 s against
# 0.01 s for d = 4 over 500 levels, and 0.14 s against under 0.005 s at
# d = 2000. A `dgmrf` route with the tridiagonal precision matches the
# speed but has to assemble a parameter-dependent sparse matrix and pay
# a sparse Cholesky for a log-determinant that is already available in
# closed form, and `hetar1`'s per-time scaling would have to be folded
# into that precision as well. The closed form is the cleaner of the two
# and the faster of the three, so it is what runs.
#
# WHY ONLY AR(1). `ou`, `exp`, `gau` and `mat` are genuinely dense
# kernels over arbitrary positions, with no banded inverse to exploit,
# and `toep`/`homtoep` have no positive-definite parameterization to
# factor in the first place; all of them stay on the dense `dmvnorm`
# path on purpose.

#' Log-density of a level-major `b` under AR(1) with correlation
#' `rho = tphi / sqrt(1 + tphi^2)` and marginal standard deviations
#' `sdv` (length `d`, recycled over levels, or length 1 for the
#' homogeneous block). O(d) per level, vectorized over levels.
#'
#' @noRd
ar1_lpdf <- function(b, sdv, tphi, d) {
  n_levels <- length(b) %/% d
  # the homogeneous block divides by a scalar, so it never builds the
  # recycled vector the heterogeneous one needs
  hom <- length(sdv) == 1L
  z <- if (hom) b / sdv else b / rep(sdv, length.out = d * n_levels)
  ldet_D <- if (hom) d * log(sdv) else sum(log(sdv))
  if (d == 1L) {
    return(n_levels * (-0.5 * log(2 * pi) - ldet_D) - 0.5 * sum(z * z))
  }
  rho <- tphi / sqrt(1 + tphi^2)
  # 1 / (1 - rho^2) under this map, without the cancelling subtraction
  inv_omr2 <- 1 + tphi^2
  # index arithmetic on the flat vector, so no advector matrix is
  # formed: level j of the level-major layout starts at (j - 1) * d
  offs <- (seq_len(n_levels) - 1L) * d
  cur <- as.vector(outer(seq.int(2L, d), offs, "+"))
  e <- z[cur] - rho * z[cur - 1L]
  z1 <- z[offs + 1L]
  n_levels * (-0.5 * d * log(2 * pi) + 0.5 * (d - 1) * log(inv_omr2) -
                ldet_D) -
    0.5 * (sum(z1 * z1) + inv_omr2 * sum(e * e))
}

covstruct_registry <- list(
  us = list(
    npar = function(dim) dim + dim * (dim - 1L) / 2L,
    sd_idx = function(dim) seq_len(dim),
    from_natural = function(sds, C, blk) {
      if (blk[["dim"]] == 1L) return(log(sds[1L]))
      c(log(sds), us_theta_cor(C))
    },
    nll = function(b, theta, blk) {
      d <- blk[["dim"]]
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
      d <- blk[["dim"]]
      if (d == 1L) {
        V <- matrix(exp(theta)^2, 1, 1)
      } else {
        sdv <- exp(theta[seq_len(d)])
        C <- us_chol_cor(theta[-seq_len(d)], d)
        V <- C * (sdv %o% sdv)
      }
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    start = function(dim) numeric(dim + dim * (dim - 1L) / 2L)
  ),
  diag = list(
    npar = function(dim) dim,
    sd_idx = function(dim) seq_len(dim),
    from_natural = function(sds, C, blk) log(sds),
    nll = function(b, theta, blk) {
      sdv <- rep(exp(theta), times = blk[["n_levels"]])
      sum(RTMB::dnorm(b, 0, sdv, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- diag(exp(theta)^2, nrow = blk[["dim"]])
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    start = function(dim) numeric(dim)
  ),
  # Student-t latents, brms `gr(g, dist = "student")`. Same theta as the
  # gaussian `us` block, plus a FIXED `blk$dist_nu`. `theta` holds the
  # SCALE, not the standard deviation. See student_lpdf().
  us_t = list(
    npar = function(dim) dim + dim * (dim - 1L) / 2L,
    sd_idx = function(dim) seq_len(dim),
    from_natural = function(sds, C, blk) {
      if (blk[["dim"]] == 1L) return(log(sds[1L]))
      c(log(sds), us_theta_cor(C))
    },
    nll = function(b, theta, blk) {
      d <- blk[["dim"]]
      if (d == 1L) return(student_lpdf_scalar(b, theta[1], blk[["dist_nu"]]))
      W <- us_chol_L(theta[-seq_len(d)], d) * exp(theta[seq_len(d)])
      student_lpdf(b, W, blk[["dist_nu"]], d, blk[["n_levels"]])
    },
    vcov = function(theta, blk) {
      d <- blk[["dim"]]
      if (d == 1L) {
        V <- matrix(exp(theta)^2, 1, 1)
      } else {
        sdv <- exp(theta[seq_len(d)])
        C <- us_chol_cor(theta[-seq_len(d)], d)
        V <- C * (sdv %o% sdv)
      }
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    start = function(dim) numeric(dim + dim * (dim - 1L) / 2L)
  ),
  diag_t = list(
    npar = function(dim) dim,
    sd_idx = function(dim) seq_len(dim),
    from_natural = function(sds, C, blk) log(sds),
    nll = function(b, theta, blk) {
      d <- blk[["dim"]]
      if (d == 1L) return(student_lpdf_scalar(b, theta[1], blk[["dist_nu"]]))
      # a diagonal SCALE matrix, but still ONE mixing variable per level
      # (probe E: brms shares `dfm` across the coefficients of a level
      # even under cor = FALSE), so this is a multivariate t and not a
      # product of d univariate ones. The quadratic form couples the
      # level's coefficients even though the scale matrix does not.
      z <- b * rep(exp(-theta), times = blk[["n_levels"]])
      dim(z) <- c(d, blk[["n_levels"]])
      student_lpdf_core(RTMB::colSums(z * z), sum(theta), blk[["dist_nu"]],
                        d, blk[["n_levels"]])
    },
    vcov = function(theta, blk) {
      V <- diag(exp(theta)^2, nrow = blk[["dim"]])
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    start = function(dim) numeric(dim)
  ),
  homdiag = list(
    npar = function(dim) 1L,
    sd_idx = function(dim) 1L,
    from_natural = function(sds, C, blk) homogeneous_sd(sds,
                                                        blk[["covstruct"]]),
    nll = function(b, theta, blk) {
      sum(RTMB::dnorm(b, 0, exp(theta), log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- diag(exp(theta)^2, nrow = blk[["dim"]])
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
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
      ar1_lpdf(b, exp(theta[1]), theta[2], blk[["dim"]])
    },
    vcov = function(theta, blk) {
      d <- blk[["dim"]]
      rho <- theta[2] / sqrt(1 + theta[2]^2)
      V <- exp(theta[1])^2 * rho^abs(outer(seq_len(d), seq_len(d), "-"))
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
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
      d <- blk[["dim"]]
      sdv <- exp(theta[seq_len(d)])
      a <- 1 / (d - 1)
      rho <- -a + (1 + a) / (1 + exp(-theta[d + 1L]))
      C <- diag(d) * (1 - rho) + rho
      Sigma <- C * (RTMB::matrix(sdv, ncol = 1) %*% RTMB::matrix(sdv, nrow = 1))
      dim(b) <- c(d, length(b) %/% d)
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      d <- blk[["dim"]]
      sdv <- exp(theta[seq_len(d)])
      a <- 1 / (d - 1)
      rho <- -a + (1 + a) / (1 + exp(-theta[d + 1L]))
      C <- diag(d) * (1 - rho) + rho
      V <- C * (sdv %o% sdv)
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    start = function(dim) numeric(dim + 1L)
  )
)

# Continuous-position AR / exponential covariance over num_factor()
# levels: Sigma_ij = sd^2 * exp(-rate * |t_i - t_j|). theta = (log sd,
# log rate); blk$aux_D holds the distance matrix.
covstruct_registry[["ou"]] <- list(
  npar = function(dim) 2L,
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    Sigma <- exp(2 * theta[1]) * exp(-exp(theta[2]) * blk[["aux_D"]])
    dim(b) <- c(blk[["dim"]], length(b) %/% blk[["dim"]])
    sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- exp(theta[1])^2 * exp(-exp(theta[2]) * blk[["aux_D"]])
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(dim) c(0, 0)
)

# Heterogeneous Toeplitz (glmmTMB parameter count 2 * dim - 1): d
# log-sds plus banded correlations rho_k = phi_k / sqrt(1 + phi_k^2).
# PD is not guaranteed for every parameter value (same as glmmTMB); the
# optimizer stays in the feasible region.
covstruct_registry[["toep"]] <- list(
  npar = function(dim) 2L * dim - 1L,
  sd_idx = function(dim) seq_len(dim),
  nll = function(b, theta, blk) {
    "[<-" <- RTMB::ADoverload("[<-")
    d <- blk[["dim"]]
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
    d <- blk[["dim"]]
    sdv <- exp(theta[seq_len(d)])
    phi <- theta[d + seq_len(d - 1L)]
    cvec <- c(1, phi / sqrt(1 + phi^2))
    C <- matrix(cvec[abs(outer(seq_len(d), seq_len(d), "-")) + 1L], d, d)
    V <- C * (sdv %o% sdv)
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(dim) numeric(2L * dim - 1L)
)

# Heterogeneous AR(1): d log-sds plus one phi (rho = phi/sqrt(1+phi^2)).
covstruct_registry[["hetar1"]] <- list(
  npar = function(dim) dim + 1L,
  sd_idx = function(dim) seq_len(dim),
  nll = function(b, theta, blk) {
    d <- blk[["dim"]]
    ar1_lpdf(b, exp(theta[seq_len(d)]), theta[d + 1L], d)
  },
  vcov = function(theta, blk) {
    d <- blk[["dim"]]
    sdv <- exp(theta[seq_len(d)])
    rho <- theta[d + 1L] / sqrt(1 + theta[d + 1L]^2)
    V <- rho^abs(outer(seq_len(d), seq_len(d), "-")) * (sdv %o% sdv)
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(dim) numeric(dim + 1L)
)

# Homogeneous compound symmetry: one log-sd plus one correlation on
# (-1/(d-1), 1).
covstruct_registry[["homcs"]] <- list(
  npar = function(dim) 2L,
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    d <- blk[["dim"]]
    a <- 1 / (d - 1)
    rho <- -a + (1 + a) / (1 + exp(-theta[2]))
    C <- diag(d) * (1 - rho) + rho
    Sigma <- exp(2 * theta[1]) * C
    dim(b) <- c(d, length(b) %/% d)
    sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
  },
  vcov = function(theta, blk) {
    d <- blk[["dim"]]
    a <- 1 / (d - 1)
    rho <- -a + (1 + a) / (1 + exp(-theta[2]))
    V <- exp(theta[1])^2 * (diag(d) * (1 - rho) + rho)
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(dim) c(0, 0)
)

# Homogeneous Toeplitz: one log-sd plus d-1 banded correlations.
covstruct_registry[["homtoep"]] <- list(
  npar = function(dim) dim,
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    "[<-" <- RTMB::ADoverload("[<-")
    d <- blk[["dim"]]
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
    d <- blk[["dim"]]
    phi <- theta[1L + seq_len(d - 1L)]
    cvec <- c(1, phi / sqrt(1 + phi^2))
    C <- matrix(cvec[abs(outer(seq_len(d), seq_len(d), "-")) + 1L], d, d)
    V <- exp(theta[1])^2 * C
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(dim) numeric(dim)
)

#' Spatial structures over num_factor() coordinates (blk$aux_D holds the
#' distance matrix), glmmTMB parameterizations: `theta = (log sd,
#' log range[, log shape])`.
#'
#' @noRd
spatial_entry <- function(corr_fn, npar_k) {
  list(
    npar = function(dim) npar_k,
    sd_idx = function(dim) 1L,
    nll = function(b, theta, blk) {
      Sigma <- exp(2 * theta[1]) * corr_fn(blk[["aux_D"]], theta)
      dim(b) <- c(blk[["dim"]], length(b) %/% blk[["dim"]])
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- exp(theta[1])^2 * corr_fn(blk[["aux_D"]], theta)
      V <- as.matrix(V)
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    start = function(dim) numeric(npar_k)
  )
}

covstruct_registry[["exp"]] <- spatial_entry(
  function(D, theta) exp(-D / exp(theta[2])), 2L
)

covstruct_registry[["gau"]] <- spatial_entry(
  function(D, theta) exp(-(D / exp(theta[2]))^2), 2L
)

# Matern correlation 2^(1-nu)/Gamma(nu) (d/range)^nu K_nu(d/range); the
# zero-distance diagonal is masked in (data-only mask, branch-free).
# Bounded internal transforms - smoothness in (0.1, 5), range floored at
# 5% of the median distance - keep every intermediate representable, so
# the optimizer cannot step into 0 * Inf territory (which is where both
# we and glmmTMB otherwise die). Compare fits at the covariance level,
# not the theta level.
covstruct_registry[["mat"]] <- spatial_entry(
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
covstruct_registry[["gr_cov"]] <- list(
  npar = function(dim) dim + dim * (dim - 1L) / 2L,
  sd_idx = function(dim) seq_len(dim),
  # the natural values describe the WITHIN-level covariance; the
  # across-level structure is the fixed matrix A
  from_natural = function(sds, C, blk) {
    if (blk[["dim"]] == 1L) return(log(sds[1L]))
    c(log(sds), us_theta_cor(C))
  },
  nll = function(b, theta, blk) {
    if (blk[["dim"]] == 1L) {
      Sigma <- exp(2 * theta[1]) * blk[["aux_A"]]
      return(sum(RTMB::dmvnorm(b, 0, Sigma, log = TRUE)))
    }
    S <- us_sigma(theta, blk[["dim"]])
    D <- blk[["dim"]] * blk[["n_levels"]]
    K <- RTMB::matrix(as.vector(blk[["aux_A"]])[blk[["aux_kron"]]$ia] *
                        as.vector(S)[blk[["aux_kron"]]$is], D, D)
    sum(RTMB::dmvnorm(b, 0, K, log = TRUE))
  },
  vcov = function(theta, blk) {
    d <- blk[["dim"]]
    V <- if (d == 1L) {
      matrix(exp(theta[1])^2, 1, 1)
    } else {
      sdv <- exp(theta[seq_len(d)])
      us_chol_cor(theta[-seq_len(d)], d) * (sdv %o% sdv)
    }
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(dim) numeric(dim + dim * (dim - 1L) / 2L)
)

# Known-precision random effects: level-major b ~ N(0, Q^-1 (x) Sigma)
# with Q a fixed (sparse) precision over the grouping levels and Sigma
# the unstructured d x d within-level covariance - gr(g, prec = Q). The
# sparse GMRF density keeps large structures (big phylogenies, spatial
# graphs) tractable where the dense gr(cov=) path is not.
#
# Correlated slopes use the precision-side Kronecker identity: the
# inverse of A (x) Sigma is A^-1 (x) Sigma^-1, so the block precision is
# Q (x) Sigma^-1, which stays as sparse as Q itself (one dense d x d
# block per stored entry of Q). It is assembled on the tape as an
# AD-weighted sum of the fixed sparse matrices Q (x) E_ab (blk$aux_Qk),
# so nothing but the d(d+1)/2 weights depends on the parameters.
covstruct_registry[["gr_prec"]] <- list(
  npar = function(dim) dim + dim * (dim - 1L) / 2L,
  sd_idx = function(dim) seq_len(dim),
  # the natural values describe the WITHIN-level covariance; the
  # across-level structure is the fixed precision Q
  from_natural = function(sds, C, blk) {
    if (blk[["dim"]] == 1L) return(log(sds[1L]))
    c(log(sds), us_theta_cor(C))
  },
  nll = function(b, theta, blk) {
    if (blk[["dim"]] == 1L) {
      Qs <- exp(-2 * theta[1]) * blk[["aux_Q"]]
      return(sum(RTMB::dgmrf(b, 0, Qs, log = TRUE)))
    }
    # RTMB::solve, not base's: the S4 advector method is not imported
    Sinv <- RTMB::solve(us_sigma(theta, blk[["dim"]]))
    # no zero to start from: an AD-weighted sparse matrix has no
    # additive identity, so the first term seeds the sum
    Qs <- NULL
    for (e in blk[["aux_Qk"]]) {
      term <- Sinv[e[["a"]], e[["b"]]] * e[["M"]]
      Qs <- if (is.null(Qs)) term else Qs + term
    }
    sum(RTMB::dgmrf(b, 0, Qs, log = TRUE))
  },
  vcov = function(theta, blk) {
    d <- blk[["dim"]]
    V <- if (d == 1L) {
      matrix(exp(theta[1])^2, 1, 1)
    } else {
      sdv <- exp(theta[seq_len(d)])
      us_chol_cor(theta[-seq_len(d)], d) * (sdv %o% sdv)
    }
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(dim) numeric(dim + dim * (dim - 1L) / 2L)
)

#' Index maps for the level-major dense Kronecker covariance
#' `A (x) Sigma`: entry `(i, j)` of the `d * n_levels` square is
#' `A[level(i), level(j)] * Sigma[coef(i), coef(j)]`, so the matrix is
#' one elementwise product of two gathers from `as.vector(A)` and
#' `as.vector(Sigma)`. The gather indices are data constants, so they
#' are built once here; on the tape only `Sigma` moves, and indexing a
#' vector keeps the advector class that `kronecker()` would strip.
#'
#' Used both for a single `gr(cov = )` term with correlated slopes and
#' for an `|ID|`-merged group, where `d` is the TOTAL merged dimension.
#'
#' @noRd
kron_cov_index <- function(d, n_levels) {
  r <- seq_len(d * n_levels)
  l1 <- (r - 1L) %/% d + 1L
  c1 <- (r - 1L) %% d + 1L
  list(ia = as.vector(outer(l1, l1,
                            function(a, b) (b - 1L) * n_levels + a)),
       is = as.vector(outer(c1, c1, function(a, b) (b - 1L) * d + a)))
}

#' Fixed sparse pieces of the level-major Kronecker precision
#' Q (x) Sigma^-1: one matrix per entry (a, b) of the d x d within-level
#' precision, so the tape only ever multiplies them by a scalar. The
#' (a, b) and (b, a) halves share a matrix because Sigma^-1 is symmetric.
#'
#' @noRd
kron_prec_parts <- function(Q, d) {
  out <- list()
  for (a in seq_len(d)) {
    for (b in seq.int(a, d)) {
      E <- matrix(0, d, d)
      E[a, b] <- 1
      E[b, a] <- 1
      M <- methods::as(methods::as(Matrix::kronecker(Q, E),
                                   "generalMatrix"), "CsparseMatrix")
      out[[length(out) + 1L]] <- list(a = a, b = b, M = M)
    }
  }
  out
}

# ------------------------------------------------------- CAR / ICAR / BYM2
#
# Spatial conditional-autoregressive fields over a fixed adjacency
# matrix - brms's car(M, gr = g, type = ). The field lives on the
# grouping levels (one intercept per location), so the block is
# dim = 1 with n_levels locations, exactly like gr(prec = Q); what
# car() adds is that the precision is BUILT from the graph and carries
# its own hyperparameters.
#
# Parameterization follows brms: theta = log(sdcar) for the intrinsic
# forms and (log sdcar, logit rho) for the two that mix, with
# tau = 1 / sdcar^2 the precision multiplier.
#
#   escar          Q = tau (D - rho W), proper for rho in (0, 1)
#   icar / esicar  Q = tau (D - W), intrinsic, sum-to-zero constrained
#   bym2           sd^2 [(1 - rho) I + (rho / scale) K^-1], the
#                  Riebler et al. scaled mixture brms implements
#
# THE NORMALIZING CONSTANT is analytic in every case, so no
# normalize-trick / on-tape log-determinant is needed and the density
# costs one sparse matrix-vector product:
#
#   escar   log|Q| = n log tau + sum_i log d_i + sum_i log(1 - rho e_i)
#           with e_i the eigenvalues of D^-1/2 W D^-1/2 (fixed data).
#   icar    the graph Laplacian L = D - W is rank n - c with c the
#           number of connected components, so log|Q|* = (n - c) log
#           tau + log|L|*. The constrained density below turns that
#           into an exact n log tau (see the constraint note).
#
# THE SUM-TO-ZERO CONSTRAINT. An intrinsic CAR is improper: L annihilates
# the indicator of each connected component, so shifting the field
# inside a component and the intercept the other way leaves the
# likelihood untouched and the ML problem is rank deficient, not merely
# ill conditioned. We adopt brms's remedy - a soft sum-to-zero
# constraint whose precision rides on tau, as it does in brms's
# non-centered zcar parameterization - so the whole block precision is
#
#   Q(tau) = tau * K,   K = L + sum_j kappa_j s_j s_j',
#   kappa_j = 1 / (con_sd n_j)^2,
#
# with s_j the indicator of component j (n_j levels). K is fixed data,
# so log|Q| = n log tau + log|K| is exact and constant-free, and the
# density is a proper Gaussian - which is what makes ranef(), predict()
# and simulate() well defined on the block. The component sums are then
# pinned at an sd of con_sd n_j sdcar rather than exactly zero, so the
# fit approaches the hard-constrained (brms esicar) likelihood as
# con_sd -> 0; `esicar` selects the same density.
#
# con_sd defaults to brms's 1e-3, so the same call is the same model
# here and there. Tightening it walks the fit onto the hard-constrained
# (esicar) likelihood quadratically, and the walk is worth knowing:
# measured on a 4 x 4 lattice against a hard sum-to-zero reference,
# 1e-3 is off by 4.7e-4 in the log-likelihood (3.6e-5 relative in
# sdcar), 1e-4 by 4.7e-6 (3.6e-7), 1e-5 by 4.5e-8 (3.7e-9), 1e-6 by
# 9.2e-10, and 1e-7 loses to roundoff (1.1e-4). The bias at the default
# is four orders below the parameter's own standard error, and the
# tighter settings cost optimizer robustness - the constraint direction
# carries a factor con_sd^-2 of the block Hessian, and over 25 lattice
# refits nlminb reported false convergence 0 times at 1e-3, once at
# 1e-4 and 6 times at 1e-5 - so the loose default is the better trade.
#
# The price of any sum-to-zero constraint is a dense rank-c update
# inside the block's Laplace Hessian, which caps the practical field
# size in the low thousands.
car_con_sd_default <- 1e-3
car_types <- c("escar", "esicar", "icar", "bym2")

#' brms's validate_car_matrix plus the level matching data_ac() does: a
#' symmetric binary adjacency over exactly the locations the data show,
#' in the block's level order.
#'
#' @noRd
car_adjacency <- function(M, locs) {
  if (length(dim(M)) != 2L || nrow(M) != ncol(M)) {
    stop("car(): M must be a square adjacency matrix", call. = FALSE)
  }
  Md <- as.matrix(M)
  # NA first: every check below compares, and a comparison against NA is
  # NA, which reaches the user as "missing value where TRUE/FALSE
  # needed" from isSymmetric() - naming neither the matrix nor the cell
  if (anyNA(Md)) {
    bad <- which(is.na(Md), arr.ind = TRUE)
    stop("car(): M has ", nrow(bad), " missing entry/entries (the first ",
         "at row ", bad[1L, 1L], ", column ", bad[1L, 2L], "); an ",
         "adjacency matrix needs a 0 for every non-neighbor pair",
         call. = FALSE)
  }
  # Locations are matched by NAME, so one set of names is enough - but
  # when both are present they have to agree, or the row a location gets
  # and the column it gets are two different places.
  rn <- rownames(Md)
  cn <- colnames(Md)
  if (is.null(rn) && is.null(cn)) {
    stop("car(): M needs dimnames naming the locations", call. = FALSE)
  }
  if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) {
    stop("car(): M's rownames and colnames must name the same locations ",
         "in the same order", call. = FALSE)
  }
  nms <- rn %||% cn
  dup <- unique(nms[duplicated(nms)])
  if (length(dup)) {
    stop("car(): M names location(s) ",
         paste(utils::head(dup, 5), collapse = ", "),
         " more than once", call. = FALSE)
  }
  miss <- setdiff(locs, nms)
  if (length(miss)) {
    stop("car(): M has no row for location(s) ",
         paste(utils::head(miss, 5), collapse = ", "),
         if (length(miss) > 5) " ..." else "", call. = FALSE)
  }
  # positional subsetting, so a matrix carrying only rownames (or only
  # colnames) is read exactly like one carrying both
  pos <- match(locs, nms)
  W <- methods::as(Matrix::Matrix(Md[pos, pos, drop = FALSE],
                                  sparse = TRUE), "generalMatrix")
  if (!Matrix::isSymmetric(W, check.attributes = FALSE)) {
    stop("car(): M must be symmetric", call. = FALSE)
  }
  if (any(Matrix::diag(W) != 0)) {
    stop("car(): M must have a zero diagonal (no location neighbors ",
         "itself)", call. = FALSE)
  }
  W <- methods::as(W, "CsparseMatrix")
  # brms's validate_car_matrix reads an adjacency as non-negative
  # weights. Binarizing a negative entry would turn a stated repulsion
  # into a neighbor, so refuse it rather than reinterpret it.
  if (any(W@x < 0)) {
    stop("car(): M has ", sum(W@x < 0), " negative entry/entries (the ",
         "smallest is ", format(min(W@x), digits = 4), "); an adjacency ",
         "matrix holds non-negative weights", call. = FALSE)
  }
  if (any(W@x != 1)) {
    message("car(): converting all non-zero values in M to 1.")
    W@x[W@x != 1] <- 1
  }
  dimnames(W) <- list(locs, locs)
  W
}

#' The three finite-element matrices, under either the INLA (M0/M1/M2)
#' or the fmesher (c0/g1/g2) spelling. The mesh size is whatever the
#' matrices say it is - see spde_node_index() for why the data cannot be
#' asked instead.
#'
#' @noRd
spde_matrices <- function(fem) {
  nms <- names(fem) %||% character(0)
  key <- if (all(c("M0", "M1", "M2") %in% nms)) c("M0", "M1", "M2")
         else if (all(c("c0", "g1", "g2") %in% nms)) c("c0", "g1", "g2")
         else NULL
  if (is.null(key)) {
    stop("spde(): fem must be a list holding M0, M1, M2 (INLA) or ",
         "c0, g1, g2 (fmesher::fm_fem)", call. = FALSE)
  }
  M0 <- fem[[key[1L]]]
  if (length(dim(M0)) != 2L || nrow(M0) != ncol(M0)) {
    stop("spde(): ", key[1L], " must be a square matrix, one row per ",
         "mesh node", call. = FALSE)
  }
  n <- nrow(M0)
  out <- list()
  for (i in seq_along(key)) {
    Mi <- fem[[key[i]]]
    if (length(dim(Mi)) != 2L || nrow(Mi) != n || ncol(Mi) != n) {
      stop("spde(): ", key[i], " is ",
           paste(dim(Mi) %||% length(Mi), collapse = " x "), " but ",
           key[1L], " has ", n, " mesh nodes; the three finite-element ",
           "matrices must be square and the same size", call. = FALSE)
    }
    out[[c("M0", "M1", "M2")[i]]] <-
      methods::as(methods::as(Matrix::Matrix(Mi, sparse = TRUE),
                              "generalMatrix"), "CsparseMatrix")
  }
  out
}

#' Node indices for an spde() grouping variable.
#'
#' The finite-element matrices are indexed by mesh ROW NUMBER and carry
#' no dimnames to match labels against, so nothing reconciles a level
#' ordering with the mesh. Deriving the ordering from the grouping
#' variable - factor levels, or sort(unique(as.character(gv))) - reads
#' integer node ids lexicographically ("1", "10", "11", "2", ...), which
#' permutes the field against its own precision and converges silently to
#' a different model. The contract is therefore explicit: gr holds mesh
#' row indices, in 1..nrow(M0). Anything else is refused rather than
#' guessed at. Unobserved nodes are fine (they keep their prior); gaps
#' are not an error.
#'
#' @noRd
spde_node_index <- function(gv, n, gr_expr) {
  idx <- spde_node_numeric(gv)
  if (is.null(idx) || any(idx != trunc(idx))) {
    stop("spde(): mesh nodes are indexed by their ROW NUMBER in the ",
         "finite-element matrices, which carry no location names, so ",
         "gr must hold whole-number node indices in 1..", n, "; '",
         deparse1(gr_expr), "' holds labels that are not whole numbers. ",
         "Map each observation onto its mesh row first and pass those ",
         "indices as gr", call. = FALSE)
  }
  idx <- as.integer(idx)
  bad <- unique(idx[idx < 1L | idx > n])
  if (length(bad)) {
    stop("spde(): gr = ", deparse1(gr_expr), " holds node index/indices ",
         paste(utils::head(sort(bad), 5), collapse = ", "),
         ", outside the mesh's 1..", n, " rows", call. = FALSE)
  }
  idx
}

#' The numeric node ids behind an spde() grouping variable, or NULL when
#' the values are not numbers. A factor is read through its LEVELS, never
#' its integer codes: the codes are a level ordering, which is exactly
#' the thing that must not decide mesh rows.
#'
#' @noRd
spde_node_numeric <- function(gv) {
  v <- if (is.factor(gv)) {
    lv <- suppressWarnings(as.numeric(levels(gv)))
    if (anyNA(lv)) return(NULL)
    lv[as.integer(gv)]
  } else if (is.numeric(gv)) {
    as.numeric(gv)
  } else {
    w <- suppressWarnings(as.numeric(as.character(gv)))
    if (anyNA(w)) return(NULL)
    w
  }
  if (anyNA(v)) NULL else v
}

#' The block level labels an spde() grouping variable maps onto, used by
#' both frame assembly and newdata prediction so the two spellings of a
#' node ("3", 3L, factor "3") always land on the same column.
#'
#' @noRd
spde_node_labels <- function(gv) {
  v <- spde_node_numeric(gv)
  if (is.null(v) || any(v != trunc(v))) return(as.character(gv))
  as.character(as.integer(v))
}

#' Number of `theta` entries a CAR block needs. The two types that mix an
#' independent part with the spatial one carry a `rho`; the intrinsic
#' forms do not.
#'
#' @noRd
car_npar <- function(type) if (type %in% c("escar", "bym2")) 2L else 1L

#' Starting `theta` for a CAR block: unit `sdcar`, and `rho` at the middle
#' of its range where the type has one.
#'
#' @noRd
car_start <- function(type) numeric(car_npar(type))

# brms bounds both mixing parameters on (0, 1)
#' Maps an unbounded `theta` entry onto a CAR mixing parameter in (0, 1).
#'
#' @noRd
car_rho <- function(x) 1 / (1 + exp(-x))

#' Connected components of a symmetric adjacency matrix, by breadth-first
#' sweep over the sparse column pattern (isolated levels are their own
#' component, which is what the rank correction counts).
#'
#' @noRd
car_components <- function(W) {
  n <- nrow(W)
  Wt <- methods::as(W, "TsparseMatrix")
  nbr <- split(c(Wt@j + 1L, Wt@i + 1L),
               factor(c(Wt@i + 1L, Wt@j + 1L), levels = seq_len(n)))
  comp <- integer(n)
  k <- 0L
  for (s in seq_len(n)) {
    if (comp[s]) next
    k <- k + 1L
    queue <- s
    comp[s] <- k
    while (length(queue)) {
      v <- queue[1L]
      queue <- queue[-1L]
      nb <- nbr[[v]]
      new <- nb[comp[nb] == 0L]
      if (length(new)) {
        comp[new] <- k
        queue <- c(queue, unique(new))
      }
    }
  }
  comp
}

#' brms's scaling factor for BYM2 (brms:::.car_scale): the geometric mean
#' of the marginal variances of a unit-scale ICAR field under the
#' sum-to-zero constraint. Reproduced here (perturbation included) so the
#' rho of a bym2 fit means what it means in brms.
#'
#' @noRd
car_scale_factor <- function(W) {
  n <- nrow(W)
  Q <- Matrix::Diagonal(n, Matrix::rowSums(W)) - W
  Q <- Q + Matrix::Diagonal(n) * max(Matrix::diag(Q)) *
    sqrt(.Machine$double.eps)
  Sigma <- Matrix::solve(Q)
  A <- matrix(1, 1, n)
  Wc <- Sigma %*% t(A)
  Sigma <- Sigma - Wc %*% solve(A %*% Wc) %*% Matrix::t(Wc)
  exp(mean(log(Matrix::diag(Sigma))))
}

#' Everything the density needs, precomputed once from the adjacency
#' matrix over the grouping levels. `type` decides which fields matter;
#' all of them are fixed data.
#'
#' @noRd
car_aux <- function(W, type, con_sd = car_con_sd_default) {
  n <- nrow(W)
  deg <- as.numeric(Matrix::rowSums(W))
  L <- methods::as(Matrix::Diagonal(n, deg) - W, "generalMatrix")
  aux <- list(type = type, n = n, W = W, L = L, deg = deg)
  if (type == "escar") {
    if (any(deg == 0)) {
      stop("car(type = \"escar\"): every location needs at least one ",
           "neighbor; location(s) ",
           paste(rownames(W)[deg == 0], collapse = ", "),
           " have none. Use type = \"icar\" instead", call. = FALSE)
    }
    isq <- diag(1 / sqrt(deg), nrow = n)
    ev <- eigen(isq %*% as.matrix(W) %*% isq, symmetric = TRUE,
                only.values = TRUE)$values
    # The normalized adjacency's spectrum lies in [-1, 1], but LAPACK
    # returns the extreme eigenvalues a few ulp outside it. The density
    # takes log(1 - rho * eig) with rho = plogis(theta2), which reaches
    # 1 to double precision by theta2 = 36: a stray eig = 1 + 2e-16
    # makes the argument negative and the whole objective NaN, and the
    # optimizer only reports "NA/NaN function evaluation". Clamping
    # costs nothing elsewhere - the determinant is unchanged to twelve
    # digits - and keeps the boundary evaluable.
    aux$eigW <- pmin(pmax(ev, -1), 1)
    aux$ldet_deg <- sum(log(deg))
    return(aux)
  }
  comp <- car_components(W)
  nj <- as.numeric(table(comp))
  aux$con_sd <- con_sd
  aux$kappa0 <- 1 / (con_sd * nj)^2
  aux$Sgrp <- Matrix::sparseMatrix(i = comp, j = seq_len(n), x = 1,
                                   dims = c(length(nj), n))
  K <- L + Matrix::t(aux$Sgrp) %*% Matrix::Diagonal(length(nj),
                                                    aux$kappa0) %*%
    aux$Sgrp
  aux$K <- methods::as(methods::as(K, "generalMatrix"), "CsparseMatrix")
  aux$ldet_K <- as.numeric(Matrix::determinant(aux$K,
                                               logarithm = TRUE)$modulus)
  aux$n_comp <- length(nj)
  if (type == "bym2") {
    # the scaled mixture needs the covariance of the unit-scale ICAR
    aux$Kinv <- as.matrix(Matrix::solve(aux$K))
    aux$scale <- car_scale_factor(W)
  }
  aux
}

#' Numeric covariance of the whole field at a theta (draws, VarCorr
#' details); off the tape, so a dense solve is fine.
#'
#' @noRd
car_cov <- function(theta, blk) {
  a <- blk[["aux_car"]]
  s2 <- exp(2 * theta[1])
  if (a[["type"]] == "escar") {
    rho <- car_rho(theta[2])
    Q <- Matrix::Diagonal(a[["n"]], a[["deg"]]) - rho * a[["W"]]
    return(s2 * as.matrix(Matrix::solve(Q)))
  }
  if (a[["type"]] == "bym2") {
    rho <- car_rho(theta[2])
    return(s2 * ((1 - rho) * diag(a[["n"]]) +
                   (rho / a[["scale"]]) * a[["Kinv"]]))
  }
  s2 * as.matrix(Matrix::solve(a[["K"]]))
}

covstruct_registry[["car"]] <- list(
  npar = function(dim) {
    stop("car npar needs the type; handled at the frame call site",
         call. = FALSE)
  },
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    a <- blk[["aux_car"]]
    n <- a[["n"]]
    if (a[["type"]] == "bym2") {
      # the mixture has no sparse precision, so the dense marginal
      # covariance is the honest form; one MVN observation
      rho <- car_rho(theta[2])
      Sigma <- RTMB::matrix(
        (1 - rho) * exp(2 * theta[1]) * as.vector(diag(n)) +
          (rho * exp(2 * theta[1]) / a[["scale"]]) *
          as.vector(a[["Kinv"]]), n, n)
      return(sum(RTMB::dmvnorm(b, 0, Sigma, log = TRUE)))
    }
    tau <- exp(-2 * theta[1])
    if (a[["type"]] == "escar") {
      rho <- car_rho(theta[2])
      quad <- sum(a[["deg"]] * b^2) -
        rho * sum(b * as.vector(a[["W"]] %*% b))
      ldet <- n * log(tau) + a[["ldet_deg"]] +
        sum(log(1 - rho * a[["eigW"]]))
    } else {
      quad <- sum(b * as.vector(a[["K"]] %*% b))
      ldet <- n * log(tau) + a[["ldet_K"]]
    }
    0.5 * (ldet - n * log(2 * pi)) - 0.5 * tau * quad
  },
  vcov = function(theta, blk) {
    # sdcar, brms's reported scale; the field covariance is car_cov()
    matrix(exp(theta[1])^2, 1, 1,
           dimnames = list("sd(car)", "sd(car)"))
  },
  start = function(dim) {
    stop("car start needs the type; handled at the frame call site",
         call. = FALSE)
  }
)

# ------------------------------------------------------------------ SPDE
#
# Gaussian Matern field through the SPDE / finite-element
# representation (Lindgren-Rue-Lindstrom): on a fixed mesh the field's
# precision is
#
#   Q(kappa, tau) = tau^2 (kappa^4 M0 + 2 kappa^2 M1 + M2)
#
# with M0/M1/M2 (fmesher's c0/g1/g2, INLA's M0/M1/M2) the mesh's finite
# element matrices - fixed sparse data. Only the three weights depend on
# the parameters, so the whole precision is an AD-weighted sum of fixed
# sparse matrices and the log-determinant comes from RTMB's sparse GMRF
# density; theta = (log tau, log kappa). alpha = 2 (nu = 1 in the plane)
# is the only order the three-matrix form covers, which is the standard
# choice and what fm_fem() returns.
#' Number of `theta` entries for an spde block: `log tau` and `log kappa`.
#'
#' @noRd
spde_npar <- function() 2L

#' Starting `theta` for an spde block: unit `tau` and unit `kappa`.
#'
#' @noRd
spde_start <- function() c(0, 0)

# Range and marginal sd follow the alpha = 2, d = 2 identities
# (Lindgren et al. 2011 eq. 2): range = sqrt(8 nu)/kappa with nu = 1,
# sigma^2 = 1 / (4 pi kappa^2 tau^2).
#' Marginal range of an spde field at a `theta`, on the natural scale.
#'
#' @noRd
spde_range <- function(theta) sqrt(8) / exp(theta[2])

#' Marginal standard deviation of an spde field at a `theta`, on the
#' natural scale.
#'
#' @noRd
spde_sd <- function(theta) {
  1 / (exp(theta[1]) * exp(theta[2]) * sqrt(4 * pi))
}

covstruct_registry[["spde"]] <- list(
  npar = function(dim) {
    stop("spde npar is fixed at 2; handled at the frame call site",
         call. = FALSE)
  },
  sd_idx = function(dim) integer(0),
  nll = function(b, theta, blk) {
    a <- blk[["aux_spde"]]
    kap2 <- exp(2 * theta[2])
    Q <- exp(2 * theta[1]) * (kap2 * kap2 * a[["M0"]] +
                                (2 * kap2) * a[["M1"]] + a[["M2"]])
    sum(RTMB::dgmrf(b, 0, Q, log = TRUE))
  },
  vcov = function(theta, blk) {
    matrix(spde_sd(theta)^2, 1, 1,
           dimnames = list("sd(spde)", "sd(spde)"))
  },
  start = function(dim) {
    stop("spde start is handled at the frame call site", call. = FALSE)
  }
)

#' Numeric precision of an spde block at a theta (draws, diagnostics).
#'
#' @noRd
spde_prec <- function(theta, blk) {
  a <- blk[["aux_spde"]]
  kap2 <- exp(2 * theta[2])
  exp(2 * theta[1]) * (kap2 * kap2 * a[["M0"]] + (2 * kap2) * a[["M1"]] +
                         a[["M2"]])
}

# Exact Gaussian process over observed positions (gp(...) without k=):
# anisotropic squared-exponential kernel
# sd^2 exp(-sum_j d_j^2 / (2 rho_j^2)), assembled from per-dimension
# squared-difference matrices (blk$aux_D2) so the tape sees data
# matrices and advector lengthscales; iso = TRUE shares one rho.
# theta = (log sd, log rho) iso / (log sd, log rho_1..log rho_D)
# otherwise; parameter count depends on the dimension count, so npar
# and start are handled at the frame call sites (gp_npar/gp_start),
# like rr.
# the 1e-6 nugget keeps the notoriously ill-conditioned SE kernel
# Cholesky-factorizable as the range grows (standard GP practice)
#' Number of `theta` entries for an exact gp block: one `log sd` plus one
#' lengthscale, shared when `iso` and one per input dimension otherwise.
#'
#' @noRd
gp_npar <- function(D, iso) if (isTRUE(iso)) 2L else 1L + as.integer(D)

#' Starting `theta` for an exact gp block: unit sd and unit lengthscales.
#'
#' @noRd
gp_start <- function(D, iso) numeric(gp_npar(D, iso))

#' Correlation matrix of an exact gp block at a `theta`, assembled on the
#' tape from the per-dimension squared-difference matrices. The nugget on
#' the diagonal keeps the factorization stable at long lengthscales.
#'
#' @noRd
gp_corr <- function(theta, blk) {
  Q <- 0
  for (j in seq_along(blk[["aux_D2"]])) {
    rho <- if (isTRUE(blk[["gp_iso"]])) exp(theta[2]) else exp(theta[1 + j])
    Q <- Q + blk[["aux_D2"]][[j]] / (2 * rho^2)
  }
  exp(-Q) + diag(1e-6, nrow(blk[["aux_D2"]][[1]]))
}

#' Numeric `K(X*, X)` of a fitted exact-gp block at new coordinates Xnew
#' (n_new x D) against the block's positions. The nugget rides on
#' coincident points so kriging at observed positions collapses to
#' exact interpolation (`K*` row = K row -> indicator weights).
#'
#' @noRd
gp_cross_cov <- function(theta, blk, Xnew, pos) {
  Q <- 0
  same <- 1
  for (j in seq_len(ncol(pos))) {
    dj <- outer(Xnew[, j], pos[, j], "-")
    rho <- if (isTRUE(blk[["gp_iso"]])) exp(theta[2]) else exp(theta[1 + j])
    Q <- Q + dj^2 / (2 * rho^2)
    same <- same * (dj == 0)
  }
  exp(2 * theta[1]) * (exp(-Q) + 1e-6 * same)
}

covstruct_registry[["gp"]] <- list(
  npar = function(dim) {
    stop("gp npar needs the dimension count; handled at the frame ",
         "call site", call. = FALSE)
  },
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    Sigma <- exp(2 * theta[1]) * gp_corr(theta, blk)
    dim(b) <- c(blk[["dim"]], length(b) %/% blk[["dim"]])
    sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- as.matrix(exp(theta[1])^2 * gp_corr(theta, blk))
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(dim) {
    stop("gp start needs the dimension count; handled at the frame ",
         "call site", call. = FALSE)
  }
)

#' Hilbert-space GP approximation (Riutort-Mayol et al.): tensor-product
#' sine basis in Z, independent coefficients whose prior SDs follow the
#' D-dim SE-kernel spectral density at the multi-index frequencies
#' (blk$aux_omega, M x D):
#' `S(w) = sd^2 (2 pi)^{D/2} prod_j rho_j exp(-0.5 sum_j rho_j^2 w_j^2)`.
#'
#' @noRd
hsgp_sds <- function(theta, omega, iso = TRUE) {
  D <- ncol(omega)
  log_rho <- if (isTRUE(iso)) rep(theta[2], D) else theta[1 + seq_len(D)]
  logS <- 2 * theta[1] + (D / 2) * log(2 * pi) + sum(log_rho) -
    0.5 * as.vector(omega^2 %*% exp(2 * log_rho))
  exp(logS / 2)
}

#' Tensor-product sine basis at centered coordinates xc (n x D) for the
#' multi-index frequencies omega (M x D) and boundaries L (length D).
#'
#' @noRd
hsgp_basis <- function(xc, omega, L) {
  Phi <- 1
  for (j in seq_len(ncol(omega))) {
    Phi <- Phi * sin(outer(xc[, j] + L[j], omega[, j])) / sqrt(L[j])
  }
  Phi
}

#' Largest pairwise Euclidean distance among the rows of P: brms's input
#' scale for gp() terms (brms:::.data_gp does
#' `dmax <- sqrt(max(diff_quad(Xgp)))`). The per-pair arithmetic here is
#' brms's, so the answer is bit-identical, but the sweep is chunked and
#' the 1-D case is closed-form so a large n never materializes the n x n
#' distance matrix that brms builds.
#'
#' @noRd
gp_max_dist <- function(P) {
  D <- ncol(P)
  if (D == 1L) {
    return(sqrt((max(P[, 1]) - min(P[, 1]))^2))
  }
  n <- nrow(P)
  step <- max(1L, as.integer(ceiling(1e6 / n)))
  best <- 0
  for (s in seq.int(1L, n, by = step)) {
    ii <- seq.int(s, min(s + step - 1L, n))
    q <- 0
    for (j in seq_len(D)) {
      q <- q + outer(P[ii, j], P[, j], function(a, b) (a - b)^2)
    }
    best <- max(best, max(q))
  }
  sqrt(best)
}

#' brms's boundary rule (brms:::choose_L): one range taken over the whole
#' centered matrix - not per column - so a shared, isotropy-preserving
#' box, floored at 1 so a domain narrower than the unit scale still gets
#' a boundary of at least c.
#'
#' @noRd
gp_choose_L <- function(xc, cvec) {
  cvec * max(1, max(xc) - min(xc))
}

#' Start values for the Hilbert-space form. Because the inputs are
#' rescaled to unit maximum pairwise distance, the domain geometry is
#' fixed and a constant start is scale-free. rho = 1 would push the
#' spectral density at the top frequencies to underflow (exp(-rho^2 w^2)
#' with w up to `k*pi/(2c)`); rho = 0.1 keeps every basis SD representable
#' while still starting well inside the smooth region.
#'
#' @noRd
hsgp_start <- function(D, iso) {
  c(0, rep(log(0.1), gp_npar(D, iso) - 1L))
}

covstruct_registry[["hsgp"]] <- list(
  npar = function(dim) {
    stop("hsgp npar needs the dimension count; handled at the frame ",
         "call site", call. = FALSE)
  },
  sd_idx = function(dim) 1L,
  nll = function(b, theta, blk) {
    sum(RTMB::dnorm(b, 0, hsgp_sds(theta, blk[["aux_omega"]], blk[["gp_iso"]]),
                    log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- diag(as.numeric(hsgp_sds(theta, blk[["aux_omega"]],
                                  blk[["gp_iso"]]))^2,
              nrow = blk[["dim"]])
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(dim) {
    stop("hsgp start needs the dimension count; handled at the frame ",
         "call site", call. = FALSE)
  }
)

# Known, fully fixed within-level covariance (glmmTMB equalto): b ~
# N(0, V) with V supplied by the user - zero parameters (meta-analysis
# sampling covariances and similar).
covstruct_registry[["equalto"]] <- list(
  npar = function(dim) 0L,
  sd_idx = function(dim) integer(0),
  nll = function(b, theta, blk) {
    dim(b) <- c(blk[["dim"]], length(b) %/% blk[["dim"]])
    sum(RTMB::dmvnorm(t(b), 0, blk[["aux_A"]], log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- blk[["aux_A"]]
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(dim) numeric(0)
)

#' Reduced-rank / factor-analytic (glmmTMB rr): per-level coefficients
#' c = Lambda f with f iid standard normal factors (the block's b segment)
#' and Lambda the D x rank loadings matrix, columnwise lower-triangular
#' (theta). Covariance = Lambda Lambda', rank-deficient by design. The
#' expansion from factors to coefficients happens in expand_b(); npar and
#' start are handled at the frame call sites (they need the rank).
#'
#' @noRd
rr_npar <- function(dim, rank) {
  as.integer(dim * rank - rank * (rank - 1L) / 2L)
}

#' start at `Lambda = [I_rank; 0]`: identified, unit factor scale
#'
#' @noRd
rr_start <- function(dim, rank) {
  th <- numeric(rr_npar(dim, rank))
  pos <- 0L
  for (j in seq_len(rank)) {
    th[pos + 1L] <- 1
    pos <- pos + (dim - j + 1L)
  }
  th
}

#' The columnwise lower-triangular loadings matrix that a reduced-rank
#' theta segment encodes (AD-safe).
#'
#' @noRd
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

covstruct_registry[["rr"]] <- list(
  npar = function(dim) {
    stop("rr npar needs the rank; handled at the frame call site",
         call. = FALSE)
  },
  sd_idx = function(dim) integer(0),
  nll = function(b, theta, blk) {
    sum(RTMB::dnorm(b, 0, 1, log = TRUE))
  },
  vcov = function(theta, blk) {
    L <- rr_loadings(theta, blk[["dim"]], blk[["rank"]])
    V <- as.matrix(L %*% t(L))
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(dim) {
    stop("rr start needs the rank; handled at the frame call site",
         call. = FALSE)
  }
)

#' Coefficient-space vector the Z matrices multiply: identical to b
#' except for rr blocks, whose factors expand through the loadings.
#' AD-safe (`RTMB::matrix` + `[<-` overload) and numeric-safe.
#'
#' @noRd
expand_b <- function(frame, b, theta) {
  if (!isTRUE(frame[["has_rr"]])) return(b)
  "[<-" <- RTMB::ADoverload("[<-")
  cvec <- rep(b[1] * 0, frame[["n_c"]])   # keeps the advector class if taped
  for (bk in frame[["re_blocks"]]) {
    if (bk[["covstruct"]] == "rr") {
      L <- rr_loadings(theta[bk[["theta_idx"]]], bk[["dim"]], bk[["rank"]])
      Fm <- RTMB::matrix(b[bk[["b_idx"]]], bk[["rank"]], bk[["n_levels"]])
      cvec[bk[["c_idx"]]] <- as.vector(L %*% Fm)
    } else {
      cvec[bk[["c_idx"]]] <- b[bk[["b_idx"]]]
    }
  }
  cvec
}

# ------------------------------------ correlation parameters ---------
#
# WHICH THETA POSITIONS HOLD A CORRELATION, and under which map. The LKJ
# prior (R/priors.R) is a density on a whole correlation matrix, so it
# needs a block's correlation parameters as a GROUP, and it needs the map
# from those unbounded parameters onto correlation space. The map differs
# by structure:
#
#   "chol"  the row-normalized unit lower-triangular Cholesky factor of a
#           d x d correlation matrix (us, us_t, gr_cov, gr_prec): every
#           entry of the matrix is free, one parameter per strictly-lower
#           position, in the column-major order us_chol_L() fills.
#   "ar1"   one parameter, rho = t / sqrt(1 + t^2) on (-1, 1).
#   "cs"    one parameter, rho = -a + (1 + a) / (1 + exp(-t)) on (-a, 1)
#           with a = 1 / (d - 1). A compound-symmetric matrix is positive
#           definite only above -1/(d - 1), so its correlation is bounded
#           away from -1 and the map says so.
#
# A structure absent from the table either has no correlation parameter
# at all or has one no LKJ density fits; the second kind is named in
# `lkj_refusals` with its reason, so that a prior addressed to it is
# refused by name instead of silently doing nothing.

#' The correlation segment of an `us`-style theta: the strictly-lower
#' Cholesky entries, which follow the `dim` log standard deviations.
#'
#' @noRd
cor_spec_chol <- function(d) {
  if (d < 2L) return(NULL)
  list(kind = "chol", d = d, idx = d + seq_len(d * (d - 1L) / 2L))
}

covstruct_registry[["us"]]$cor_spec <- cor_spec_chol
covstruct_registry[["us_t"]]$cor_spec <- cor_spec_chol
covstruct_registry[["gr_cov"]]$cor_spec <- cor_spec_chol
covstruct_registry[["gr_prec"]]$cor_spec <- cor_spec_chol
covstruct_registry[["cs"]]$cor_spec <- function(d) {
  if (d < 2L) return(NULL)
  list(kind = "cs", d = d, idx = d + 1L)
}
covstruct_registry[["homcs"]]$cor_spec <- function(d) {
  if (d < 2L) return(NULL)
  list(kind = "cs", d = d, idx = 2L)
}
covstruct_registry[["ar1"]]$cor_spec <- function(d) {
  if (d < 2L) return(NULL)
  list(kind = "ar1", d = d, idx = 2L)
}
covstruct_registry[["hetar1"]]$cor_spec <- function(d) {
  if (d < 2L) return(NULL)
  list(kind = "ar1", d = d, idx = d + 1L)
}

# Structures whose correlations exist but carry no LKJ density.
lkj_refusals <- c(
  toep = paste("a banded parameterization that is not positive definite",
               "everywhere, so there is no correlation matrix over the",
               "whole parameter space to put a density on"),
  homtoep = paste("a banded parameterization that is not positive",
                  "definite everywhere, so there is no correlation",
                  "matrix over the whole parameter space to put a",
                  "density on")
)

#' A block's correlation parameters (`kind`, `d`, and the `idx`
#' positions WITHIN the block's theta segment), or NULL when it has
#' none.
#'
#' @noRd
block_cor_spec <- function(bk) {
  f <- covstruct_registry[[bk[["covstruct"]]]][["cor_spec"]]
  if (is.null(f)) return(NULL)
  f(bk[["dim"]])
}

#' How many of a block's theta positions are correlations.
#'
#' @noRd
block_n_cor <- function(bk) length(block_cor_spec(bk)[["idx"]] %||% integer(0))

# The three readers below exist so that `covstruct_registry` and
# `lkj_refusals` can stay internal while an extension package still asks
# the questions the non-centering plan and the sampling default priors
# ask (dev/draws-extraction.md). Exporting the registry would have made
# every field of every structure's entry into stable API; these are the
# three questions actually asked, and nothing else.

#' Whether any Cholesky factor is registered for a block's structure,
#' which is the first half of non-centering eligibility.
#'
#' @noRd
covstruct_has_chol <- function(bk) {
  reg <- covstruct_registry[[bk[["covstruct"]]]]
  !is.null(reg[["chol_sd"]]) || !is.null(reg[["chol_L"]])
}

#' The positions WITHIN a block's theta segment that are standard
#' deviations.
#'
#' @noRd
block_sd_idx <- function(bk) {
  covstruct_registry[[bk[["covstruct"]]]]$sd_idx(bk[["dim"]])
}

#' What correlation prior a block can carry: `"none"` (it has no
#' correlation parameters), `"lkj"` (an LKJ density fits them), or
#' `"unsupported"` (it has them, but its parameterization has no
#' correlation matrix over the whole of it for a density to be about).
#'
#' Three-valued because both answers are needed and they are not
#' complements: a caller choosing a default wants `"lkj"`, and a caller
#' listing the slots left flat wants `"unsupported"` alone - a block
#' with no correlations at all is not a gap.
#'
#' @noRd
block_cor_prior <- function(bk) {
  # the refusal is tested FIRST: a structure on the refusal list need
  # not register a cor_spec at all (toep does not), and asking for the
  # spec first would report it as having no correlation rather than as
  # having one nothing can prior
  if (bk[["covstruct"]] %in% names(lkj_refusals)) return("unsupported")
  if (is.null(block_cor_spec(bk))) return("none")
  "lkj"
}

# --------------------------------- non-centered parameterization -----
#
# `frm_sample(reparameterize = TRUE)` samples z ~ N(0, I) in place of a
# block's own coefficients and computes b = L(theta) z on the tape. The
# centered joint posterior of (b, theta) is a funnel (the width of b's
# prior is itself being sampled), and NUTS cannot adapt one step size to
# both ends of it. The non-centered pair (z, theta) is a product of
# independent-ish pieces, which is why brms writes every block that way.
#
# The ML fit is untouched by this. The Laplace approximation integrates
# b out, and that integral is invariant under a linear change of the
# integrated variable, so the fitted objective, its tape, and every
# number the fit reports are the same either way.
#
# A structure joins the lane by declaring ONE factor accessor:
#
#   chol_sd(theta, blk)  the per-level standard deviations, when the
#                        factor is diagonal (length `blk$dim`, recycled
#                        over the levels, level-major like `b`), or
#   chol_L(theta, blk)   the per-level lower-triangular `L` with
#                        `L L' = Sigma`, the within-level covariance,
#
# plus, for a structure whose LEVELS are correlated, the fixed
# across-level factor
#
#   chol_A(blk)          lower-triangular `LA` with `LA LA' = A`.
#
# The whole block is then `B = L Z LA'` on the d x n_levels layout,
# which is `vec(B) = (LA (x) L) vec(Z)`: the Cholesky factor of the
# Kronecker covariance, assembled from its two small factors instead of
# factorized as one big one.
#
# One accessor serves both directions on purpose: the same map builds b
# on the tape and back-transforms each posterior draw afterwards, and a
# second implementation of it would be a second model.
#
# A structure with NO accessor samples centered, and `frm_sample()`
# names it. The absentees and their reasons:
#
#   us_t, diag_t  a Student-t latent is a SCALE MIXTURE of gaussians;
#                 non-centering it needs the mixing variable sampled
#                 too, which is a different construction, not a factor.
#   car, spde,
#   gr_prec       the density is a sparse PRECISION. The factor of its
#                 inverse is dense, so a solve per leapfrog step would
#                 cost more than the funnel does.
#   gp            a dense kernel over the observed positions: a full
#                 n x n factorization per gradient, on the tape.
#   ou, exp, gau,
#   mat           the same, over a spatial field's positions.
#   toep, homtoep the banded parameterization does not guarantee a
#                 positive definite matrix (glmmTMB's does not either),
#                 so no factor exists over the whole parameter space and
#                 `b = L z` would not be a bijection there.
#   rr            already non-centered: its `b` IS a standard normal
#                 factor vector, expanded through the loadings by
#                 `expand_b()`. Nothing to transform, nothing to report.

#' Lower Cholesky factor of the equicorrelation matrix
#' `(1 - rho) I + rho J`.
#'
#' Every row below the diagonal repeats one value per column, so the
#' factor follows an O(d) recursion rather than the general O(d^3) one:
#' with `t_j` the squared norm of row j's entries left of the diagonal,
#' `L[j, j] = sqrt(1 - t_j)` and the shared below-diagonal entry is
#' `(rho - t_j) / L[j, j]`. AD-safe.
#'
#' @noRd
cs_chol_cor <- function(rho, d) {
  "[<-" <- RTMB::ADoverload("[<-")
  L <- RTMB::matrix(rep(rho * 0, d * d), d, d)
  L[1L, 1L] <- 1
  if (d > 1L) L[seq.int(2L, d), 1L] <- rho
  tj <- rho * rho
  for (j in seq_len(d - 1L) + 1L) {
    dj <- sqrt(1 - tj)
    L[j, j] <- dj
    if (j < d) {
      cj <- (rho - tj) / dj
      L[seq.int(j + 1L, d), j] <- cj
      tj <- tj + cj * cj
    }
  }
  L
}

#' Lower Cholesky factor of the AR(1) correlation matrix `rho^|i-j|`.
#'
#' Analytic, so the factor costs no factorization: row `i` holds
#' `rho^(i-1)` in column 1 and `rho^(i-j) sqrt(1 - rho^2)` in columns
#' `2..i`. Built from data index matrices and one advector power vector,
#' so the tape sees `d^2` multiplications and no branch. AD-safe.
#'
#' @noRd
ar1_chol_cor <- function(rho, d) {
  "[<-" <- RTMB::ADoverload("[<-")
  # sequential products, as the ar1 density itself builds them: `pow`
  # would route a negative rho through exp/log
  pw <- rep(rho, d)
  pw[1] <- 1
  for (k in seq_len(d - 1L) + 1L) pw[k] <- pw[k - 1L] * rho
  ii <- row(diag(d))
  jj <- col(diag(d))
  keep <- as.vector(ii >= jj) * 1
  first <- as.vector(jj == 1L) * 1
  idx <- as.vector(pmax(ii - jj, 0L)) + 1L
  s <- sqrt(1 - rho * rho)
  RTMB::matrix(pw[idx] * keep * (first + s * (keep - first)), d, d)
}

#' The `us` factor: `Sigma = D C D` with `C = Lr Lr'`, so
#' `(D Lr)(D Lr)' = Sigma` and the recycling multiplies row `i` by
#' `sd_i`. Shared with `gr_cov`, whose within-level covariance is the
#' same unstructured one.
#'
#' @noRd
us_chol_scaled <- function(theta, d) {
  if (d == 1L) return(RTMB::matrix(exp(theta[1]), 1L, 1L))
  us_chol_L(theta[-seq_len(d)], d) * exp(theta[seq_len(d)])
}

covstruct_registry[["us"]]$chol_L <- function(theta, blk) {
  us_chol_scaled(theta, blk[["dim"]])
}
covstruct_registry[["diag"]]$chol_sd <- function(theta, blk) exp(theta)
covstruct_registry[["homdiag"]]$chol_sd <- function(theta, blk) {
  rep(exp(theta[1]), blk[["dim"]])
}
covstruct_registry[["hsgp"]]$chol_sd <- function(theta, blk) {
  hsgp_sds(theta, blk[["aux_omega"]], blk[["gp_iso"]])
}
covstruct_registry[["equalto"]]$chol_L <- function(theta, blk) {
  t(chol(blk[["aux_A"]]))
}
covstruct_registry[["cs"]]$chol_L <- function(theta, blk) {
  d <- blk[["dim"]]
  a <- 1 / (d - 1)
  rho <- -a + (1 + a) / (1 + exp(-theta[d + 1L]))
  cs_chol_cor(rho, d) * exp(theta[seq_len(d)])
}
covstruct_registry[["homcs"]]$chol_L <- function(theta, blk) {
  d <- blk[["dim"]]
  a <- 1 / (d - 1)
  rho <- -a + (1 + a) / (1 + exp(-theta[2]))
  cs_chol_cor(rho, d) * exp(theta[1])
}
covstruct_registry[["ar1"]]$chol_L <- function(theta, blk) {
  rho <- theta[2] / sqrt(1 + theta[2]^2)
  ar1_chol_cor(rho, blk[["dim"]]) * exp(theta[1])
}
covstruct_registry[["hetar1"]]$chol_L <- function(theta, blk) {
  d <- blk[["dim"]]
  rho <- theta[d + 1L] / sqrt(1 + theta[d + 1L]^2)
  ar1_chol_cor(rho, d) * exp(theta[seq_len(d)])
}
# gr(cov = A): the covariance is the Kronecker product A (x) Sigma, and
# a Kronecker product's factor is the Kronecker product of the factors.
# A is fixed data, so its factor is a constant; only the small
# within-level one moves with theta.
covstruct_registry[["gr_cov"]]$chol_L <- function(theta, blk) {
  us_chol_scaled(theta, blk[["dim"]])
}
covstruct_registry[["gr_cov"]]$chol_A <- function(blk) {
  t(chol(as.matrix(blk[["aux_A"]])))
}

# WHICH BLOCKS ARE NON-CENTERED, and why it is not all of them.
#
# The funnel is made by the SCALE: `b | sd ~ N(0, sd^2 C)` narrows as
# `sd` shrinks, and pulling the factor out of the density removes it.
# A CORRELATION parameter makes no funnel (it is bounded, and its
# effect on the width of `b` is bounded with it), so non-centering one
# buys no geometry of its own. Both kinds of parameter still have to be
# SAFE to hand a non-centered chain, and that is a question about the
# prior, not about the factor: remove the funnel and the chain is free
# to walk any flat tail the prior leaves open.
#
# Before 0.39 a correlated block stayed centered for exactly that
# reason. frmtmb parameterizes a correlation by an unbounded
# row-normalized Cholesky `theta`, and FLAT on that theta is
# `(1 - rho^2)^-3/2` on the correlation: improper, with all its mass at
# |rho| = 1. On sleepstudy `(Days | Subject)` the PROFILE log-likelihood
# is flat in that theta beyond |theta| ~ 100 and only 4.4 nats below the
# peak (dev/benchmarks.md), so the posterior really did have infinite
# mass out there; a non-centered chain walked straight down it to
# theta = 2e6 with a bulk-ESS of 1. The LKJ prior (R/priors.R) closes
# that tail, the formula route now applies `lkj(1)` by default, and the
# rule that was written to survive this change did:
#
#   a block is non-centered only when every parameter it has is either a
#   standard deviation or a correlation with a registered factor, AND
#   every one of those parameters carries a prior (ncp_plan()).
#
# The second half is the whole of the safety argument. A block with a
# flat prior anywhere in its theta is still centered, correlated or not.

#' The non-centering a block gets: `"full"`, or `NA` for a block that
#' stays centered.
#'
#' @noRd
ncp_mode <- function(bk) {
  if (is_student_block(bk)) return(NA_character_)
  cs <- bk[["covstruct"]]
  reg <- covstruct_registry[[cs]]
  if (is.null(reg[["chol_sd"]]) && is.null(reg[["chol_L"]])) {
    return(NA_character_)
  }
  n_sd <- length(reg$sd_idx(bk[["dim"]]))
  n_cor <- block_n_cor(bk)
  n_par <- switch(
    cs,
    # a fixed covariance has no parameters at all, so its factor is a
    # constant and there is nothing to expose
    equalto = 0L,
    # hsgp's non-scale parameters are LENGTHSCALES, not correlations:
    # they enter the block's standard deviations, so the factor is
    # diagonal and the whole of it comes out
    hsgp = n_sd,
    tryCatch(reg$npar(bk[["dim"]]), error = function(e) NA_integer_)
  )
  # every parameter the block has must be one the factor consumes and
  # the prior lane can reach; a structure that grows a third kind stays
  # out of the lane until it says which kind it is
  if (is.na(n_par) || n_par != n_sd + n_cor) return(NA_character_)
  "full"
}

#' Whether a block can be sampled non-centered at all.
#'
#' @noRd
ncp_eligible <- function(bk) !is.na(ncp_mode(bk))

#' `b = S(theta) z` for one block, level-major in and level-major out.
#' Runs on the tape (advector `z` and `theta`) and off it (the per-draw
#' back-transform), which is the point of the shared accessor: the same
#' map has to build `b` in both places or they describe two models.
#'
#' @noRd
ncp_scale_b <- function(bk, z, theta) {
  reg <- covstruct_registry[[bk[["covstruct"]]]]
  if (!is.null(reg[["chol_sd"]])) {
    return(z * rep(reg[["chol_sd"]](theta, bk), times = bk[["n_levels"]]))
  }
  L <- reg[["chol_L"]](theta, bk)
  dim(z) <- c(bk[["dim"]], bk[["n_levels"]])
  B <- L %*% z
  if (!is.null(reg[["chol_A"]])) B <- B %*% t(reg[["chol_A"]](bk))
  as.vector(B)
}

#' `z = S(theta)^-1 b`, the inverse map. Numeric only: it exists to put
#' the sampler's starting z where the ML mode's b is, so that a
#' non-centered chain starts at the same point in the model as a
#' centered one.
#'
#' @noRd
ncp_unscale_b <- function(bk, b, theta) {
  reg <- covstruct_registry[[bk[["covstruct"]]]]
  if (!is.null(reg[["chol_sd"]])) {
    return(b / rep(reg[["chol_sd"]](theta, bk), times = bk[["n_levels"]]))
  }
  L <- reg[["chol_L"]](theta, bk)
  B <- matrix(b, bk[["dim"]], bk[["n_levels"]])
  Z <- solve(L, B)
  if (!is.null(reg[["chol_A"]])) {
    Z <- t(solve(reg[["chol_A"]](bk), t(Z)))
  }
  as.vector(Z)
}

#' The block's dense per-level factor at a numeric `theta`, whichever
#' accessor the structure declares. For the tests that recompute
#' `b = L z` by hand, and for reporting.
#'
#' @noRd
ncp_block_chol <- function(bk, theta) {
  reg <- covstruct_registry[[bk[["covstruct"]]]]
  if (!is.null(reg[["chol_sd"]])) {
    return(diag(as.numeric(reg[["chol_sd"]](theta, bk)),
                nrow = bk[["dim"]]))
  }
  as.matrix(reg[["chol_L"]](theta, bk))
}

# Smooth wiggly blocks are iid-Gaussian with one variance (the inverse
# smoothing parameter); reuse the homdiag machinery under its own name so
# blocks stay self-describing.
covstruct_registry[["smooth"]] <- covstruct_registry[["homdiag"]]
