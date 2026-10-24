# Registry of random-effect covariance structures. One entry per
# structure, and ONE entry shape: the built-in structures below and a
# structure another package registers with frmtmb_register_covstruct()
# (R/covstruct-seam.R) are read by the same code through cs_entry(). An
# entry supplies:
#
#   npar(blk)           number of parameters in `theta` for the block
#   start(blk)          starting values for theta
#   sd_idx(blk)         positions in theta that are log standard
#                       deviations
#   nll(b, theta, blk)  AD log-DENSITY (the name is historical) of the
#                       block's random effects, summed, vectorized over
#                       levels (never loop over levels)
#   vcov(theta, blk)    numeric within-level covariance for VarCorr()
#
# and optionally from_natural(), cor_spec(blk), chol_sd(), chol_L(),
# chol_A() (each documented where it is used), plus:
#
#   natural(theta, blk) the rows confint_varcorr() reports, for a block
#                       that is not a brms grouping-factor block (a
#                       smooth, a gp(), a car() or spde() field): its
#                       theta is not a set of sd_<group>__<coef> and
#                       cor_ parameters. See cs_natural_rows().
#   draw(theta, blk)    one draw of the block's whole `b`, for a block
#                       whose levels are correlated; the default draws
#                       each level from vcov()
#   levels_coupled      TRUE when a fixed matrix correlates the LEVELS
#                       (gr(cov =), gr(prec =), car(), spde()): no
#                       unseen-level law, one unit for the sandwich, no
#                       per-group integral
#   expand              list(apply, jacobian, condsd) for a block whose
#                       coefficients are a fixed linear map of `b`
#                       (esicar); active on a block with `expands = TRUE`
#   term                the predictor special that builds the block, for
#                       a structure that is not written as a bar term:
#                       list(args, parse, build, levels_of,
#                       levels_complete, usage)
#
# Parameterizations follow glmmTMB's covstruct vignette so that fitted
# covariances are directly comparable. `b` for a block is level-major: the
# `dim` coefficients of each level are contiguous (mkReTrms ordering), so
# matrix(b, nrow = dim) has one level per column.

#' An entry function of the block that reads only its dimension.
#'
#' Most structures' parameter count and layout depend on the block's
#' size alone; a field or a reduced-rank block reads more of the block,
#' which is why the entry functions take the block and not its size.
#'
#' @noRd
by_dim <- function(f) {
  force(f)
  function(blk) f(blk[["dim"]])
}

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
    frm_stop("The requested correlation matrix is not positive definite",
             call. = FALSE)
  })
  # column-major recycling divides row i by diag(Lc)[i]
  L <- Lc / diag(Lc)
  L[lower.tri(L)]
}

#' Unstructured d x d covariance from its theta segment (AD-safe).
#' `floored = TRUE` reads the standard deviations through `sd_floored()`,
#' as the random-effect blocks do; a family's residual covariance keeps
#' `exp()`.
#'
#' @noRd
us_sigma <- function(theta, d, floored = FALSE) {
  if (d == 1L) return(RTMB::matrix(exp(2 * theta[1]), 1, 1))
  sdv <- if (floored) {
    sd_floored(theta[seq_len(d)])
  } else {
    exp(theta[seq_len(d)])
  }
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
    frm_stop("A '", what, "' block has one shared standard deviation; got ",
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
#' independent univariate ones. So is a `||` term, which is one `diag()`
#' block. brms builds one term the same way.
#'
#' All the Student-t terms of one grouping factor also share one `u_j`
#' per level, as in brms. An example is
#' `(1 | gr(g, dist = "student")) + (0 + x | gr(g, dist = "student"))`.
#' A term in `mu` and a term in `sigma` also share it. At each level,
#' the coefficients of these terms are one multivariate t, with a
#' block-diagonal scale matrix. So the example is the model
#' `(1 + x || gr(g, dist = "student"))`.
#'
#' A `gr(g, by = f, dist = "student")` term has one block per level of
#' `f`. Each block holds different levels of `g`, so alone the blocks
#' share nothing. Beside a plain Student-t term on `g`, each level of
#' `g` shares one `u_j` with the by-term's block that holds it.
#'
#' A term on another grouping factor has its own `u_j`. A gaussian term
#' on `g` has none. The Student-t terms of one factor take one
#' `dist_nu`, because brms has one `df` per factor.
#' `simulate(re_formula = )` keeps or redraws the terms together at the
#' levels they share. `quadrature = TRUE` refuses terms that share a
#' level, because it integrates one random effect at a time.
#'
#' **`Sigma` is the SCALE, not the covariance.** The variance is
#' `Sigma * nu / (nu - 2)`, so a standard deviation is
#' `scale * sqrt(nu / (nu - 2))`. `VarCorr()` reports the scale in its
#' `sd` rows, as brms's `sd_` is; `confint()`, `variables()`
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

#' Mark the Student-t blocks that share one latent scale.
#'
#' brms 2.23.0 gives every Student-t term of a grouping factor ONE
#' mixing variable per level of the factor: `(1 | gr(g, dist =
#' "student")) + (0 + x | gr(g, dist = "student"))` writes `dfm_1` into
#' both terms. So does a term in `mu` beside one in `sigma`, and a
#' `gr(g, by = f)` term beside a plain one, where `dfm_1` runs over
#' every level of `g` (dev/prilabel-log/brms-ref.txt,
#' dev/prilabel-rev-log/misc.txt). A term on another factor gets its
#' own `dfm_2`. At each level of the factor the coefficients of every
#' block that holds the level are then ONE multivariate t, with a
#' block-diagonal scale matrix.
#'
#' The scale is shared per LEVEL. frmtmb splits a `gr(g, by = f)` term
#' into one block per by-level, each over its own levels of `g`, so the
#' blocks of one factor need not hold the same levels. Each block of a
#' set gets `t_share` (the indices of every block of the set, in block
#' order), `t_level` (its levels' positions in the union of the set's
#' levels) and `t_n` (the size of that union). A set whose blocks hold
#' disjoint levels, one by-term alone, shares nothing and is not
#' marked: it is the model of independent blocks. brms has one `df` per
#' factor, so a marked set must agree on `dist_nu`.
#'
#' @noRd
student_share_sets <- function(re_blocks) {
  st <- which(vapply(re_blocks, is_student_block, NA))
  if (length(st) < 2L) return(re_blocks)
  key <- vapply(re_blocks[st], function(bk) bk[["group_name"]] %||% "", "")
  for (g in unique(key[duplicated(key)])) {
    ids <- st[key == g]
    lvs <- lapply(re_blocks[ids], function(bk) as.character(bk[["levels"]]))
    u <- unique(unlist(lvs))
    held <- tabulate(unlist(lapply(lvs, match, u)), nbins = length(u))
    if (max(held) < 2L) next
    lab <- vapply(re_blocks[ids], function(bk) bk[["term_label"]], "")
    nu <- vapply(re_blocks[ids], function(bk) bk[["dist_nu"]], 0)
    if (length(unique(nu)) > 1L) {
      frm_stop("The Student-t terms of grouping factor '", g, "' (",
               paste(lab, collapse = "; "), ") share one latent scale, as ",
               "in brms, which has one df per grouping factor, so they ",
               "take one dist_nu. Got ",
               paste(unique(nu), collapse = " and "),
               ". Write the same dist_nu in each", call. = FALSE)
    }
    for (k in seq_along(ids)) {
      i <- ids[k]
      re_blocks[[i]][["t_share"]] <- ids
      re_blocks[[i]][["t_level"]] <- match(lvs[[k]], u)
      re_blocks[[i]][["t_n"]] <- length(u)
    }
  }
  re_blocks
}

#' The pieces of one Student-t block's density that a shared latent
#' scale adds up: per level, the quadratic form `q_j = b_j' S^-1 b_j`,
#' and the log-determinant of the scale's factor, exactly as the block's
#' own `nll` computes them (scalar, `us_t` and `diag_t` alike).
#'
#' @noRd
student_t_parts <- function(b, theta, bk) {
  d <- bk[["dim"]]
  n <- bk[["n_levels"]]
  if (d == 1L) {
    z <- b * exp(-theta[1L])
    return(list(q = z * z, ldet = theta[1L]))
  }
  if (identical(bk[["covstruct"]], "diag_t")) {
    z <- b * rep(exp(-theta), times = n)
    dim(z) <- c(d, n)
    return(list(q = RTMB::colSums(z * z), ldet = sum(theta)))
  }
  W <- us_chol_L(theta[-seq_len(d)], d) * sd_floored(theta[seq_len(d)])
  dim(b) <- c(d, n)
  z <- RTMB::solve(W, b)
  list(q = RTMB::colSums(z * z), ldet = sum(log(RTMB::diag(W))))
}

#' Log-density of the blocks of one shared Student-t set: at each level
#' `j` of the factor, one multivariate t over the coefficients of every
#' block that holds `j`, of dimension `D_j`, the sum of those blocks'
#' dimensions. Its block-diagonal scale makes the quadratic forms and
#' the log-determinants add, so the density is
#'
#'   sum_j { lgamma((nu + D_j) / 2) - lgamma(nu / 2)
#'          - D_j / 2 log(nu pi) - (nu + D_j) / 2 log1p(q_j / nu) }
#'   - sum_k n_k ldet_k,
#'
#' with `q_j` the sum of the blocks' quadratic forms at `j` and `n_k`
#' block `k`'s number of levels. `b` and `theta` are the model's whole
#' vectors.
#'
#' @noRd
student_shared_lpdf <- function(bks, b, theta) {
  "[<-" <- RTMB::ADoverload("[<-")
  nu <- bks[[1L]][["dist_nu"]]
  n_u <- bks[[1L]][["t_n"]]
  D <- numeric(n_u)
  q <- numeric(n_u)
  ldet <- 0
  for (bk in bks) {
    p <- student_t_parts(b[bk[["b_idx"]]], theta[bk[["theta_idx"]]], bk)
    at <- bk[["t_level"]]
    # scattered into the union's levels; positions the block does not
    # hold stay 0, which adds nothing to the level's quadratic form
    qk <- numeric(n_u)
    qk[at] <- p$q
    q <- q + qk
    D[at] <- D[at] + bk[["dim"]]
    ldet <- ldet + length(at) * p$ldet
  }
  sum(lgamma_shift_diff(nu / 2, D / 2) - D / 2 * (log(nu) + log(pi))) -
    ldet - sum((nu + D) / 2 * log1p(q / nu))
}

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
#' The head term needs the same care and for the same reason.
#' `lgamma((nu + d)/2) - lgamma(nu/2)` written out cancels: measured
#' against a 300-bit reference on `d = 4` over 30 levels it is wrong by
#' 6.7e-05 at `nu = 1e10` and by 2.6 at `nu = 1e14`, which puts a false
#' stationary point in the gaussian limit this comment sends users to.
#' `lgamma_shift_diff()` carries the numbers.
#'
#' @noRd
student_lpdf_core <- function(q, ldet_W, nu, d, n_levels) {
  n_levels * (lgamma_shift_diff(nu / 2, d / 2) -
                d / 2 * (log(nu) + log(pi)) - ldet_W) -
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

covstruct_core <- list(
  us = list(
    npar = by_dim(function(dim) dim + dim * (dim - 1L) / 2L),
    sd_idx = by_dim(function(dim) seq_len(dim)),
    from_natural = function(sds, C, blk) {
      if (blk[["dim"]] == 1L) return(log(sds[1L]))
      c(log(sds), us_theta_cor(C))
    },
    nll = function(b, theta, blk) {
      d <- blk[["dim"]]
      if (d == 1L) {
        return(sum(RTMB::dnorm(b, 0, exp(theta), log = TRUE)))
      }
      sdv <- sd_floored(theta[seq_len(d)])
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
    start = by_dim(function(dim) numeric(dim + dim * (dim - 1L) / 2L))
  ),
  diag = list(
    npar = by_dim(function(dim) dim),
    sd_idx = by_dim(function(dim) seq_len(dim)),
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
    start = by_dim(function(dim) numeric(dim))
  ),
  # Student-t latents, brms `gr(g, dist = "student")`. Same theta as the
  # gaussian `us` block, plus a FIXED `blk$dist_nu`. `theta` holds the
  # SCALE, not the standard deviation. See student_lpdf().
  us_t = list(
    npar = by_dim(function(dim) dim + dim * (dim - 1L) / 2L),
    sd_idx = by_dim(function(dim) seq_len(dim)),
    from_natural = function(sds, C, blk) {
      if (blk[["dim"]] == 1L) return(log(sds[1L]))
      c(log(sds), us_theta_cor(C))
    },
    nll = function(b, theta, blk) {
      d <- blk[["dim"]]
      if (d == 1L) return(student_lpdf_scalar(b, theta[1], blk[["dist_nu"]]))
      W <- us_chol_L(theta[-seq_len(d)], d) * sd_floored(theta[seq_len(d)])
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
    start = by_dim(function(dim) numeric(dim + dim * (dim - 1L) / 2L))
  ),
  diag_t = list(
    npar = by_dim(function(dim) dim),
    sd_idx = by_dim(function(dim) seq_len(dim)),
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
    start = by_dim(function(dim) numeric(dim))
  ),
  homdiag = list(
    npar = by_dim(function(dim) 1L),
    sd_idx = by_dim(function(dim) 1L),
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
    start = by_dim(function(dim) 0)
  ),
  # Homogeneous AR(1) over the term's factor levels, glmmTMB
  # parameterization: theta = (log sd, phi) with rho = phi/sqrt(1+phi^2).
  ar1 = list(
    npar = by_dim(function(dim) 2L),
    sd_idx = by_dim(function(dim) 1L),
    nll = function(b, theta, blk) {
      ar1_lpdf(b, sd_floored(theta[1]), theta[2], blk[["dim"]])
    },
    vcov = function(theta, blk) {
      d <- blk[["dim"]]
      rho <- theta[2] / sqrt(1 + theta[2]^2)
      V <- exp(theta[1])^2 * rho^abs(outer(seq_len(d), seq_len(d), "-"))
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    start = by_dim(function(dim) c(0, 0))
  ),
  # Heterogeneous compound symmetry, glmmTMB-style parameter count (d+1):
  # d log-sds plus one correlation mapped onto (-1/(d-1), 1).
  cs = list(
    npar = by_dim(function(dim) dim + 1L),
    sd_idx = by_dim(function(dim) seq_len(dim)),
    nll = function(b, theta, blk) {
      d <- blk[["dim"]]
      sdv <- sd_floored(theta[seq_len(d)])
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
    start = by_dim(function(dim) numeric(dim + 1L))
  )
)

#' The variance `sd^2` of a block with one scale, from its log sd, with
#' a floor where `exp()` underflows.
#'
#' Below a log sd of about -372, `exp(2 * log_sd)` is exactly 0, the
#' block's covariance is the zero matrix, and RTMB's `dmvnorm()` of a
#' nonzero `b` against it returned `+Inf`, a log density no point can
#' have. The fit never goes there, but a sampler's first long leapfrog
#' step does: `frm_sample(fit)` on an exact `y ~ gp(x)` fit stepped to
#' log sd -1138, read the `+Inf` as a point of infinite density, and
#' stayed there with every later transition divergent. With the floor
#' the zero matrix is `1e-300` times the correlation, and a `b` that is
#' not exactly zero has a log density near `-|b|^2 * 1e300`, or `-Inf`
#' where that overflows, which a sampler rejects.
#'
#' The floor is outside every cancelling sum (rtmb-pitfalls item 11),
#' and it is below half an ulp of any variance above 9e-285, so the sum
#' is the variance itself, bit for bit, wherever a fit can be: the
#' objective of every fit is unchanged.
#'
#' @noRd
sd2_floored <- function(log_sd) exp(2 * log_sd) + 1e-300

#' The standard deviations of a block with several scales, from their
#' log sds, with a floor where `exp()` underflows.
#'
#' `sd2_floored()` serves a block whose covariance is one variance
#' times a correlation. A block with a scale per coefficient (`us` and
#' `us_t` with two or more coefficients, `cs`, `toep`, `hetar1`,
#' `gr(cov =)` and `gr(prec =)` with correlated slopes) builds its
#' covariance from the standard deviations themselves, and so does
#' `ar1` through its standardized values. Where one log sd fell below
#' about -745, the covariance had a zero row: `dmvnorm()` gave `+Inf`
#' (`us`, `cs`, `toep`, `gr(cov =)`), the AR(1) density `NaN`
#' (`Inf - Inf`), and the solves of `us_t` and `gr(prec =)` stopped
#' (dev/moparity-findings.md, item 5).
#'
#' `1e-150` is the square root of `sd2_floored()`'s floor, so a floored
#' variance is the same `1e-300`. It is below half an ulp of any sd of
#' 1.1e-134 or more (log sd -308.45), so there the sum is the sd itself,
#' bit for bit, and so are the objective and its gradient.
#'
#' @noRd
sd_floored <- function(log_sd) exp(log_sd) + 1e-150

# Continuous-position AR / exponential covariance over num_factor()
# levels: Sigma_ij = sd^2 * exp(-rate * |t_i - t_j|). theta = (log sd,
# log rate); blk$aux_D holds the distance matrix.
covstruct_core[["ou"]] <- list(
  npar = by_dim(function(dim) 2L),
  sd_idx = by_dim(function(dim) 1L),
  nll = function(b, theta, blk) {
    Sigma <- sd2_floored(theta[1]) * exp(-exp(theta[2]) * blk[["aux_D"]])
    dim(b) <- c(blk[["dim"]], length(b) %/% blk[["dim"]])
    sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- exp(theta[1])^2 * exp(-exp(theta[2]) * blk[["aux_D"]])
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = by_dim(function(dim) c(0, 0))
)

# Heterogeneous Toeplitz (glmmTMB parameter count 2 * dim - 1): d
# log-sds plus banded correlations rho_k = phi_k / sqrt(1 + phi_k^2).
# PD is not guaranteed for every parameter value (same as glmmTMB); the
# optimizer stays in the feasible region.
covstruct_core[["toep"]] <- list(
  npar = by_dim(function(dim) 2L * dim - 1L),
  sd_idx = by_dim(function(dim) seq_len(dim)),
  nll = function(b, theta, blk) {
    "[<-" <- RTMB::ADoverload("[<-")
    d <- blk[["dim"]]
    sdv <- sd_floored(theta[seq_len(d)])
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
  start = by_dim(function(dim) numeric(2L * dim - 1L))
)

# Heterogeneous AR(1): d log-sds plus one phi (rho = phi/sqrt(1+phi^2)).
covstruct_core[["hetar1"]] <- list(
  npar = by_dim(function(dim) dim + 1L),
  sd_idx = by_dim(function(dim) seq_len(dim)),
  nll = function(b, theta, blk) {
    d <- blk[["dim"]]
    ar1_lpdf(b, sd_floored(theta[seq_len(d)]), theta[d + 1L], d)
  },
  vcov = function(theta, blk) {
    d <- blk[["dim"]]
    sdv <- exp(theta[seq_len(d)])
    rho <- theta[d + 1L] / sqrt(1 + theta[d + 1L]^2)
    V <- rho^abs(outer(seq_len(d), seq_len(d), "-")) * (sdv %o% sdv)
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = by_dim(function(dim) numeric(dim + 1L))
)

# Homogeneous compound symmetry: one log-sd plus one correlation on
# (-1/(d-1), 1).
covstruct_core[["homcs"]] <- list(
  npar = by_dim(function(dim) 2L),
  sd_idx = by_dim(function(dim) 1L),
  nll = function(b, theta, blk) {
    d <- blk[["dim"]]
    a <- 1 / (d - 1)
    rho <- -a + (1 + a) / (1 + exp(-theta[2]))
    C <- diag(d) * (1 - rho) + rho
    Sigma <- sd2_floored(theta[1]) * C
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
  start = by_dim(function(dim) c(0, 0))
)

# Homogeneous Toeplitz: one log-sd plus d-1 banded correlations.
covstruct_core[["homtoep"]] <- list(
  npar = by_dim(function(dim) dim),
  sd_idx = by_dim(function(dim) 1L),
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
    Sigma <- sd2_floored(theta[1]) * C
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
  start = by_dim(function(dim) numeric(dim))
)

#' Spatial structures over num_factor() coordinates (blk$aux_D holds the
#' distance matrix), glmmTMB parameterizations: `theta = (log sd,
#' log range[, log shape])`.
#'
#' @noRd
spatial_entry <- function(corr_fn, npar_k) {
  list(
    npar = by_dim(function(dim) npar_k),
    sd_idx = by_dim(function(dim) 1L),
    nll = function(b, theta, blk) {
      Sigma <- sd2_floored(theta[1]) * corr_fn(blk[["aux_D"]], theta)
      dim(b) <- c(blk[["dim"]], length(b) %/% blk[["dim"]])
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- exp(theta[1])^2 * corr_fn(blk[["aux_D"]], theta)
      V <- as.matrix(V)
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    start = by_dim(function(dim) numeric(npar_k))
  )
}

covstruct_core[["exp"]] <- spatial_entry(
  function(D, theta) exp(-D / exp(theta[2])), 2L
)

covstruct_core[["gau"]] <- spatial_entry(
  function(D, theta) exp(-(D / exp(theta[2]))^2), 2L
)

# Matern correlation 2^(1-nu)/Gamma(nu) (d/range)^nu K_nu(d/range); the
# zero-distance diagonal is masked in (data-only mask, branch-free).
# Bounded internal transforms - smoothness in (0.1, 5), range floored at
# 5% of the median distance - keep every intermediate representable, so
# the optimizer cannot step into 0 * Inf territory (which is where both
# we and glmmTMB otherwise die). Compare fits at the covariance level,
# not the theta level.
covstruct_core[["mat"]] <- spatial_entry(
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
covstruct_core[["gr_cov"]] <- list(
  npar = by_dim(function(dim) dim + dim * (dim - 1L) / 2L),
  sd_idx = by_dim(function(dim) seq_len(dim)),
  # the natural values describe the WITHIN-level covariance; the
  # across-level structure is the fixed matrix A
  from_natural = function(sds, C, blk) {
    if (blk[["dim"]] == 1L) return(log(sds[1L]))
    c(log(sds), us_theta_cor(C))
  },
  nll = function(b, theta, blk) {
    if (blk[["dim"]] == 1L) {
      Sigma <- sd2_floored(theta[1]) * blk[["aux_A"]]
      return(sum(RTMB::dmvnorm(b, 0, Sigma, log = TRUE)))
    }
    S <- us_sigma(theta, blk[["dim"]], floored = TRUE)
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
  start = by_dim(function(dim) numeric(dim + dim * (dim - 1L) / 2L))
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
covstruct_core[["gr_prec"]] <- list(
  npar = by_dim(function(dim) dim + dim * (dim - 1L) / 2L),
  sd_idx = by_dim(function(dim) seq_len(dim)),
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
    Sinv <- RTMB::solve(us_sigma(theta, blk[["dim"]], floored = TRUE))
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
  start = by_dim(function(dim) numeric(dim + dim * (dim - 1L) / 2L))
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

# Exact Gaussian process over observed positions (gp(...) without k=):
# anisotropic squared-exponential kernel
# sd^2 exp(-sum_j d_j^2 / (2 rho_j^2)), assembled from per-dimension
# squared-difference matrices (blk$aux_D2) so the tape sees data
# matrices and advector lengthscales; iso = TRUE shares one rho.
# theta = (log sd, log rho) iso / (log sd, log rho_1..log rho_D)
# otherwise; parameter count depends on the dimension count, so npar
# and start are handled at the frame call sites (gp_npar/gp_start),
# like rr.
# the nugget (gp_nugget) keeps the notoriously ill-conditioned SE kernel
# Cholesky-factorizable as the range grows (standard GP practice)
#' The exact gp()'s nugget, on the CORRELATION: white noise of
#' `gp_nugget * sd^2` at each position, part of the fitted model. One
#' constant, read by the density, the kriging and the sampler's draw of
#' the field, so a change to it changes all three together; brms's
#' Stan program uses 1e-12 absolute (dev/brms-likelihood-tests.md).
#'
#' @noRd
gp_nugget <- 1e-6

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
  exp(-Q) + diag(gp_nugget, nrow(blk[["aux_D2"]][[1]]))
}

#' The rows confint_varcorr() reports for a gp() or hsgp() block: the
#' marginal sd and one range per input dimension (one when iso), each
#' on the log scale.
#'
#' hsgp estimates the lengthscale on brms's rescaled inputs, but the
#' reported range belongs in data units. The scale factor is a data
#' constant, so the shift on the log scale is exact and the se rides
#' through unchanged. The exact gp keeps the raw scale (dmax NULL).
#'
#' @noRd
gp_natural <- function(theta, blk) {
  log_dmax <- log(blk[["gp_dmax"]] %||% 1)
  nr <- length(theta) - 1L
  rng <- if (nr == 1L) "range(gp)" else {
    paste0("range(gp, ", blk[["gp_vars"]][seq_len(nr)], ")")
  }
  list(term = c("sd(gp)", rng), type = c("sd", rep("range", nr)),
       est = c(theta[1L], theta[1L + seq_len(nr)] + log_dmax),
       grad = diag(length(theta)))
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
  exp(2 * theta[1]) * (exp(-Q) + gp_nugget * same)
}

covstruct_core[["gp"]] <- list(
  # the parameter count follows the input dimension count, not the
  # block dimension (the positions)
  npar = function(blk) gp_npar(blk[["gp_D"]], blk[["gp_iso"]]),
  sd_idx = by_dim(function(dim) 1L),
  nll = function(b, theta, blk) {
    Sigma <- sd2_floored(theta[1]) * gp_corr(theta, blk)
    dim(b) <- c(blk[["dim"]], length(b) %/% blk[["dim"]])
    sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- as.matrix(exp(theta[1])^2 * gp_corr(theta, blk))
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(blk) gp_start(blk[["gp_D"]], blk[["gp_iso"]]),
  natural = function(theta, blk) gp_natural(theta, blk)
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

covstruct_core[["hsgp"]] <- list(
  npar = function(blk) gp_npar(blk[["gp_D"]], blk[["gp_iso"]]),
  sd_idx = by_dim(function(dim) 1L),
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
  start = function(blk) hsgp_start(blk[["gp_D"]], blk[["gp_iso"]]),
  natural = function(theta, blk) gp_natural(theta, blk)
)

# Known, fully fixed within-level covariance (glmmTMB equalto): b ~
# N(0, V) with V supplied by the user - zero parameters (meta-analysis
# sampling covariances and similar).
covstruct_core[["equalto"]] <- list(
  npar = by_dim(function(dim) 0L),
  sd_idx = by_dim(function(dim) integer(0)),
  nll = function(b, theta, blk) {
    dim(b) <- c(blk[["dim"]], length(b) %/% blk[["dim"]])
    sum(RTMB::dmvnorm(t(b), 0, blk[["aux_A"]], log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- blk[["aux_A"]]
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = by_dim(function(dim) numeric(0))
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

covstruct_core[["rr"]] <- list(
  npar = function(blk) rr_npar(blk[["dim"]], blk[["rank"]]),
  sd_idx = by_dim(function(dim) integer(0)),
  nll = function(b, theta, blk) {
    sum(RTMB::dnorm(b, 0, 1, log = TRUE))
  },
  vcov = function(theta, blk) {
    L <- rr_loadings(theta, blk[["dim"]], blk[["rank"]])
    V <- as.matrix(L %*% t(L))
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  start = function(blk) rr_start(blk[["dim"]], blk[["rank"]])
)

#' Whether a block's coefficients are a fixed linear map of its
#' parameters, applied by its structure's `expand` slot: the
#' hard-constrained intrinsic CAR, whose coefficients are the CENTERED
#' parameters. One predicate, because the frame flag and expand_b() have
#' to agree about which blocks need the pass. The block says so
#' (`expands`, set where the term is built) and the structure says how;
#' a block that says so under a structure with no `expand` is refused
#' rather than evaluated as the identity, which would be another model.
#'
#' @noRd
block_expands <- function(bk) {
  if (!isTRUE(bk[["expands"]])) return(FALSE)
  if (is.null(cs_entry(bk[["covstruct"]])[["expand"]])) {
    frm_stop("internal: the '", bk[["covstruct"]], "' block ",
             bk[["term_label"]] %||% "", " maps its parameters, and its ",
             "structure registers no expand slot", call. = FALSE)
  }
  TRUE
}

#' Whether a block's coefficients are its parameters POSITIONWISE, so
#' that `cvec[c_idx[i]]` is a function of `b[b_idx[i]]` alone.
#'
#' This is the property every caller that perturbs, bounds or batches
#' `b` one entry at a time depends on, and `expand_b()` is where it can
#' fail: an `rr` block's coefficients are the loadings times that
#' level's factors, and an esicar block's are `b` with each connected
#' component's mean removed, so one `b` entry reaches every coefficient
#' of its component. Everything else copies `b` across unchanged.
#'
#' Stated as one predicate beside `expand_b()` rather than tested where
#' it is needed, because the two have to agree and a caller that
#' spelled the test itself got it wrong: the finite-difference route
#' used `length(c_idx) != length(b_idx)` as a proxy for "not
#' positionwise", which catches `rr` below full rank and misses both
#' `rr` AT full rank and esicar, where the lengths are equal and the
#' map is not the identity. That was a silent wrong `Est.Error` on an
#' esicar fit (`dev/reunc-findings.md`, punch round 2).
#'
#' @noRd
block_b_positionwise <- function(bk) {
  !identical(bk[["covstruct"]], "rr") && !block_expands(bk)
}

#' Whether a frame's `b` and its coefficient vector are different
#' objects, so that `expand_b()` has to run. DERIVED from the blocks,
#' with the frame's own flag only as a fast path.
#'
#' Deriving it rather than trusting `frame$has_expand` is the point.
#' `[[` does not partially match, so a frame serialized before that
#' field existed returns NULL, and a gate that took NULL for FALSE
#' would evaluate an esicar block WITHOUT centering while its density
#' still carried the constrained normalizer: two different models, no
#' error, and only off the mode, which is exactly the regime
#' `imp_frozen_proposal()` and `cluster_scores_at()` work in. The
#' blocks cannot go stale the way a cached boolean can, because the
#' same build sets `expands` and chooses the density.
#'
#' @noRd
frame_needs_expand <- function(frame) {
  if (isTRUE(frame[["has_expand"]]) || isTRUE(frame[["has_rr"]])) {
    return(TRUE)
  }
  bks <- frame[["re_blocks"]] %||% list()
  # the same disjunction as block_b_positionwise(), and it is called
  # rather than repeated so the two cannot drift
  any(vapply(bks, function(bk) !block_b_positionwise(bk), logical(1)))
}

#' Coefficient-space vector the Z matrices multiply: identical to b
#' except for rr blocks, whose factors expand through the loadings, and
#' esicar blocks, whose field is b with each component's mean removed.
#' AD-safe (`RTMB::matrix` + `[<-` overload) and numeric-safe.
#'
#' Both exceptions have to happen HERE rather than in the block density,
#' because the vector that reaches the linear predictor is the
#' coefficient vector and not the parameter vector, and for esicar the
#' whole point is that the predictor never sees the component means.
#'
#' @noRd
expand_b <- function(frame, b, theta) {
  if (!frame_needs_expand(frame)) return(b)
  "[<-" <- RTMB::ADoverload("[<-")
  cvec <- rep(b[1] * 0, frame[["n_c"]])   # keeps the advector class if taped
  for (bk in frame[["re_blocks"]]) {
    if (bk[["covstruct"]] == "rr") {
      L <- rr_loadings(theta[bk[["theta_idx"]]], bk[["dim"]], bk[["rank"]])
      Fm <- RTMB::matrix(b[bk[["b_idx"]]], bk[["rank"]], bk[["n_levels"]])
      cvec[bk[["c_idx"]]] <- as.vector(L %*% Fm)
    } else if (block_expands(bk)) {
      ex <- cs_entry(bk[["covstruct"]])[["expand"]]
      cvec[bk[["c_idx"]]] <- ex[["apply"]](b[bk[["b_idx"]]], bk)
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

covstruct_core[["us"]]$cor_spec <- by_dim(cor_spec_chol)
covstruct_core[["us_t"]]$cor_spec <- by_dim(cor_spec_chol)
covstruct_core[["gr_cov"]]$cor_spec <- by_dim(cor_spec_chol)
covstruct_core[["gr_prec"]]$cor_spec <- by_dim(cor_spec_chol)
covstruct_core[["cs"]]$cor_spec <- by_dim(function(d) {
  if (d < 2L) return(NULL)
  list(kind = "cs", d = d, idx = d + 1L)
})
covstruct_core[["homcs"]]$cor_spec <- by_dim(function(d) {
  if (d < 2L) return(NULL)
  list(kind = "cs", d = d, idx = 2L)
})
covstruct_core[["ar1"]]$cor_spec <- by_dim(function(d) {
  if (d < 2L) return(NULL)
  list(kind = "ar1", d = d, idx = 2L)
})
covstruct_core[["hetar1"]]$cor_spec <- by_dim(function(d) {
  if (d < 2L) return(NULL)
  list(kind = "ar1", d = d, idx = d + 1L)
})

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
  f <- cs_entry(bk[["covstruct"]])[["cor_spec"]]
  if (is.null(f)) return(NULL)
  f(bk)
}

# The structures that are diagonal by construction; block_has_cor().
uncorrelated_structures <- c("diag", "diag_t", "homdiag")

#' Does a block correlate its columns?
#'
#' Not when its structure is diagonal by construction (`diag`,
#' `diag_t`, `homdiag`): every correlation is then zero. brms's
#' `cor = FALSE` block, which is what a `||` term is, has no `cor_`
#' parameter, so such a block names, reports and summarizes none
#' either. The test is the structure and not its parameter count:
#' `equalto()` estimates nothing and still has a correlation, a fixed
#' one.
#'
#' @noRd
block_has_cor <- function(bk) {
  bk[["dim"]] > 1L && !bk[["covstruct"]] %in% uncorrelated_structures
}

#' How many of a block's theta positions are correlations.
#'
#' @noRd
block_n_cor <- function(bk) length(block_cor_spec(bk)[["idx"]] %||% integer(0))

# The three readers below exist so that the registry and
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
  reg <- cs_entry(bk[["covstruct"]])
  !is.null(reg[["chol_sd"]]) || !is.null(reg[["chol_L"]])
}

#' The positions WITHIN a block's theta segment that are standard
#' deviations.
#'
#' @noRd
block_sd_idx <- function(bk) {
  cs_entry(bk[["covstruct"]])$sd_idx(bk)
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
#   ou, exp, gau,
#   mat           a dense covariance over a spatial field's positions: a
#                 full n x n factorization per gradient, on the tape.
#                 An exact gp() pays that cost (see its chol_L below):
#                 its centered field did not mix.
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

covstruct_core[["us"]]$chol_L <- function(theta, blk) {
  us_chol_scaled(theta, blk[["dim"]])
}
covstruct_core[["diag"]]$chol_sd <- function(theta, blk) exp(theta)
covstruct_core[["homdiag"]]$chol_sd <- function(theta, blk) {
  rep(exp(theta[1]), blk[["dim"]])
}
covstruct_core[["hsgp"]]$chol_sd <- function(theta, blk) {
  hsgp_sds(theta, blk[["aux_omega"]], blk[["gp_iso"]])
}
covstruct_core[["equalto"]]$chol_L <- function(theta, blk) {
  t(chol(blk[["aux_A"]]))
}
covstruct_core[["cs"]]$chol_L <- function(theta, blk) {
  d <- blk[["dim"]]
  a <- 1 / (d - 1)
  rho <- -a + (1 + a) / (1 + exp(-theta[d + 1L]))
  cs_chol_cor(rho, d) * exp(theta[seq_len(d)])
}
covstruct_core[["homcs"]]$chol_L <- function(theta, blk) {
  d <- blk[["dim"]]
  a <- 1 / (d - 1)
  rho <- -a + (1 + a) / (1 + exp(-theta[2]))
  cs_chol_cor(rho, d) * exp(theta[1])
}
covstruct_core[["ar1"]]$chol_L <- function(theta, blk) {
  rho <- theta[2] / sqrt(1 + theta[2]^2)
  ar1_chol_cor(rho, blk[["dim"]]) * exp(theta[1])
}
covstruct_core[["hetar1"]]$chol_L <- function(theta, blk) {
  d <- blk[["dim"]]
  rho <- theta[d + 1L] / sqrt(1 + theta[d + 1L]^2)
  ar1_chol_cor(rho, d) * exp(theta[seq_len(d)])
}
# gr(cov = A): the covariance is the Kronecker product A (x) Sigma, and
# a Kronecker product's factor is the Kronecker product of the factors.
# A is fixed data, so its factor is a constant; only the small
# within-level one moves with theta.
covstruct_core[["gr_cov"]]$chol_L <- function(theta, blk) {
  us_chol_scaled(theta, blk[["dim"]])
}
covstruct_core[["gr_cov"]]$chol_A <- function(blk) {
  t(chol(as.matrix(blk[["aux_A"]])))
}

# An exact gp(): the field over its distinct positions has covariance
# sd^2 K(lscale), K the kernel with gp_nugget on its diagonal, so its
# factor is sd times the Cholesky factor of K, taken on the tape. brms
# samples the same field as sdgp * L zgp with standard normal zgp. The
# factorization is n^3 / 3 per gradient for n positions, which is what
# the centered field avoided; the centered field is a funnel the sampler
# does not mix in (on lane surface's y ~ gp(x), 4 chains of 1000: 1286
# of 2000 transitions at the maximum tree depth, R-hat 1.52, bulk ESS
# 7; dev/moparity-log/gp-cmp.txt, dev/moparity-findings.md, item 7).
covstruct_core[["gp"]]$chol_L <- function(theta, blk) {
  t(chol(gp_corr(theta, blk))) * exp(theta[1])
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
  reg <- cs_entry(cs)
  if (is.null(reg[["chol_sd"]]) && is.null(reg[["chol_L"]])) {
    return(NA_character_)
  }
  n_sd <- length(reg$sd_idx(bk))
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
    # an exact gp()'s length scales enter its factor too, through the
    # kernel, and class lscale reaches them, as it reaches hsgp's,
    # so the whole of it comes out
    gp = n_sd,
    tryCatch(reg$npar(bk), error = function(e) NA_integer_)
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
  reg <- cs_entry(bk[["covstruct"]])
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
  reg <- cs_entry(bk[["covstruct"]])
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
  reg <- cs_entry(bk[["covstruct"]])
  if (!is.null(reg[["chol_sd"]])) {
    return(diag(as.numeric(reg[["chol_sd"]](theta, bk)),
                nrow = bk[["dim"]]))
  }
  as.matrix(reg[["chol_L"]](theta, bk))
}

# Smooth wiggly blocks are iid-Gaussian with one variance (the inverse
# smoothing parameter); reuse the homdiag machinery under its own name so
# blocks stay self-describing.
covstruct_core[["smooth"]] <- covstruct_core[["homdiag"]]
# brms has no sd_<group>__<coef> for a smooth: its one parameter is the
# smoothing sd, which confint_varcorr() reports under its own name
covstruct_core[["smooth"]]$natural <- function(theta, blk) {
  list(term = "sd(wiggle)", type = "sd", est = theta[1L],
       grad = matrix(1, 1L, 1L))
}

# gr(cov = A) and gr(prec = Q) correlate the LEVELS through the fixed
# matrix, so a draw of the block is one draw of the whole field. The
# within-level covariance is the registry vcov(); the across-level
# structure is A (x) Sigma, or the precision Q (x) Sigma^-1.
covstruct_core[["gr_cov"]]$levels_coupled <- TRUE
covstruct_core[["gr_cov"]]$draw <- function(theta, blk) {
  S <- cs_entry("gr_cov")$vcov(theta, blk)
  K <- kronecker(blk[["aux_A"]], S)
  drop(crossprod(chol(K), stats::rnorm(nrow(K))))
}
covstruct_core[["gr_prec"]]$levels_coupled <- TRUE
covstruct_core[["gr_prec"]]$draw <- function(theta, blk) {
  # x = U^-1 z with U'U = Q has covariance Q^-1; for correlated
  # slopes Q is the level-major Kronecker precision Q (x) Sigma^-1
  Qb <- if (blk[["dim"]] == 1L) {
    exp(-2 * theta[1L]) * blk[["aux_Q"]]
  } else {
    S <- cs_entry("gr_prec")$vcov(theta, blk)
    Matrix::kronecker(blk[["aux_Q"]], methods::as(solve(unname(S)),
                                                 "generalMatrix"))
  }
  U <- Matrix::chol(Qb)
  as.vector(Matrix::solve(U, stats::rnorm(length(blk[["b_idx"]]))))
}

# ------------------------------------------------------------ the registry
#
# The built-in entries above are the registry's first state; a package
# adds its own with frmtmb_register_covstruct() (R/covstruct-seam.R).
# An environment, so that a registration from another package's
# .onLoad() reaches a sealed namespace, keyed like the addition-term
# registry and saved and restored the same way by the test helper.

frmtmb_covstruct_registry <- new.env(parent = emptyenv())
frmtmb_covstruct_registry$reg <- covstruct_core

#' Every registered structure, built-in first, in registration order.
#'
#' @noRd
cs_entries <- function() frmtmb_covstruct_registry$reg

#' The names of every registered structure.
#'
#' @noRd
cs_names <- function() names(frmtmb_covstruct_registry$reg)

#' One structure's entry. A name nothing registered is refused here,
#' and when an extension supplies the structure the refusal names it: a
#' fit or a frame built with that package loaded and read without it
#' would otherwise fail deep inside a method on a NULL entry.
#'
#' @noRd
cs_entry <- function(name) {
  e <- frmtmb_covstruct_registry$reg[[name]]
  if (!is.null(e)) return(e)
  # loading the extension that owns the name registers it
  home <- extension_special_home(name)
  if (!is.null(home)) {
    require_extension(home, paste0("The ", name, "() structure"))
    e <- frmtmb_covstruct_registry$reg[[name]]
    if (!is.null(e)) return(e)
  }
  frm_stop("The ", name, "() structure is not registered in this session. ",
           "A fit that uses it needs the package that registered it ",
           "loaded, or the registration repeated",
           class = "frmtmb_extension_missing", call. = FALSE)
}

#' Whether a block is a field rather than a brms grouping-factor block:
#' its structure reports its own rows (`natural`), so its theta is not a
#' set of `sd_<group>__<coef>` and `cor_` parameters, brms counts no
#' grouping factor for it, and VarCorr() leaves it out.
#'
#' @noRd
cs_is_field <- function(bk) {
  !is.null(cs_entry(bk[["covstruct"]])[["natural"]])
}

#' Whether a fixed matrix correlates the block's LEVELS, so that an
#' unseen level has no marginal law, the block is one unit of the
#' random-effect prior, and the integral does not split over groups.
#'
#' @noRd
cs_levels_coupled <- function(bk) {
  isTRUE(cs_entry(bk[["covstruct"]])[["levels_coupled"]])
}

#' The natural-scale rows of a field block at `theta`: `term`, `type`
#' (`"sd"`, `"range"`, `"prop"`, `"cor"` or `"raw"`, which decides the
#' back-transform), `est` on the transformed scale (log, logit, Fisher-z
#' or none), and `grad`, the jacobian of `est` in theta, one row per
#' term. Checked here, once, so that a registered structure's slip
#' reaches its author as a message and not as a wrong interval.
#'
#' @noRd
cs_natural_rows <- function(bk, theta) {
  f <- cs_entry(bk[["covstruct"]])[["natural"]]
  if (is.null(f)) return(NULL)
  nr <- f(theta, bk)
  k <- length(nr[["term"]])
  g <- nr[["grad"]]
  if (!is.character(nr[["term"]]) || !k ||
      length(nr[["type"]]) != k || length(nr[["est"]]) != k ||
      !all(nr[["type"]] %in% c("sd", "range", "prop", "cor", "raw")) ||
      !is.matrix(g) || nrow(g) != k || ncol(g) != length(theta)) {
    frm_stop("The natural = slot of the '", bk[["covstruct"]], "' structure ",
             "must return list(term, type, est, grad): one entry per row, ",
             "type one of sd, range, prop, cor, raw, and grad a matrix with ",
             "one row per term and one column per theta entry (",
             length(theta), ")", call. = FALSE)
  }
  nr
}

#' The standard error of one natural row from the covariance of its
#' block's theta. A row that reads one theta entry takes that entry's
#' own variance, so a NaN elsewhere in the block cannot reach it.
#'
#' @noRd
cs_natural_se <- function(g, V) {
  kk <- which(g != 0)
  if (!length(kk)) return(0)
  if (length(kk) == 1L) return(abs(g[kk]) * sqrt(V[kk, kk]))
  sqrt(max(drop(g[kk] %*% V[kk, kk, drop = FALSE] %*% g[kk]), 0))
}

#' The level labels a field block's grouping values map onto, the same
#' spelling at frame assembly and on newdata. A structure whose levels
#' are not the values' own strings (spde()'s mesh row numbers) says so
#' through its term's `levels_of`.
#'
#' @noRd
cs_levels_of <- function(bk, gv) {
  f <- cs_entry(bk[["covstruct"]])[["term"]][["levels_of"]]
  if (is.null(f)) as.character(gv) else as.character(f(gv))
}

#' Whether a field block's levels are its whole index set (every mesh
#' node of an spde() block, observed or not), so that a value outside
#' them is not an unseen level of a distribution but a position the
#' structure does not have.
#'
#' @noRd
cs_levels_complete <- function(bk) {
  isTRUE(cs_entry(bk[["covstruct"]])[["term"]][["levels_complete"]])
}
