# Cluster-robust (sandwich) covariance for frmtmb fits.
#
# Why per-CLUSTER and never per-OBSERVATION: the objective is a
# MARGINAL likelihood. An observation's contribution to it is not
# defined once a random effect is integrated out, so `sandwich::estfun`
# has nothing to return and the package does not provide one. A
# CLUSTER's contribution is defined, as long as the marginal likelihood
# factors over clusters - which it does exactly when every random
# effect (and every grouped likelihood term) is nested in, or equal to,
# the clustering factor. All of the guards below exist to establish
# that condition before any number is computed.
#
# How the scores are computed (see dev/sandwich/probe-scores.R for the
# two candidates and their timings). Give the objective one extra
# parameter per cluster, `clw`, multiplying that cluster's per-row data
# weights. Evaluating the gradient at clw = e_g leaves an objective
# that IS cluster g's marginal negative log-likelihood: the clusters
# whose weight is zero contribute their random effects' Gaussian prior
# integrated over the whole space, which is exactly 1, so their
# Laplace term is exactly 0 and carries exactly zero gradient. So
# `obj$gr(c(theta_hat, e_g))` read on the theta block IS the cluster
# score, to machine precision, from ONE tape and G gradient
# evaluations.

#' Small-sample factors, in the clubSandwich spelling and definitions
#' (Pustejovsky & Tipton 2018, table 1). clubSandwich applies the
#' square root to the estimating matrices; these are the factors on the
#' variance.
#'
#' @noRd
cr_adjust <- function(type, G, N, p) {
  switch(type,
         CR0 = 1,
         CR1 = G / (G - 1),
         CR1p = G / (G - p),
         CR1S = G * (N - 1) / ((G - 1) * (N - p)),
         frm_stop("unreachable", call. = FALSE))
}

#' Turns the user's `cluster` into a factor of length `n_obs`, in the
#' row order of the fitted model frame.
#'
#' @noRd
resolve_cluster <- function(fit, cluster) {
  n <- fit$frame[["n_obs"]]
  if (inherits(cluster, "formula")) {
    if (length(cluster) != 2L) {
      frm_stop("`cluster` must be a one-sided formula such as ~ g",
               call. = FALSE)
    }
    if (length(all.vars(cluster)) > 1L) {
      frm_stop("`cluster` names more than one variable. Multiway ",
               "clustering is not implemented; pass one factor, e.g. ",
               "cluster = ~ interaction(firm, year), for a single crossed ",
               "clustering factor", call. = FALSE)
    }
    ex <- cluster[[2L]]
    env <- environment(cluster) %||% parent.frame()
    v <- tryCatch(eval(ex, fit$frame[["data_frame"]], env),
                  error = function(e) NULL)
    if (is.null(v) || length(v) != n) {
      # a clustering variable that is in the data but not in the model
      # formula is not in the stored frame; go back to the call's data
      # and apply the fit's own na.action so the rows still line up
      d0 <- tryCatch(eval(fit$call$data, env), error = function(e) NULL)
      v2 <- if (is.null(d0)) NULL else {
        tryCatch(eval(ex, d0, env), error = function(e) NULL)
      }
      if (!is.null(v2) && !is.null(fit$frame[["na_action"]])) {
        v2 <- v2[-unclass(fit$frame[["na_action"]])]
      }
      v <- v2 %||% v
    }
    if (is.null(v)) {
      frm_stop("could not find `", deparse1(ex), "` in the model frame or ",
               "in the data the model was fitted to; pass the clustering ",
               "factor itself as `cluster`", call. = FALSE)
    }
  } else if (is.character(cluster) && length(cluster) == 1L &&
             cluster %in% names(fit$frame[["data_frame"]])) {
    v <- fit$frame[["data_frame"]][[cluster]]
  } else {
    v <- cluster
  }
  if (length(v) != n) {
    frm_stop("`cluster` has length ", length(v), " but the model was ",
             "fitted to ", n, " rows", call. = FALSE)
  }
  if (anyNA(v)) {
    frm_stop("`cluster` cannot contain missing values", call. = FALSE)
  }
  droplevels(as.factor(v))
}

#' Independent units of the random-effect prior, as column sets of the
#' coefficient space. For nearly every covariance structure the levels
#' of a block are independent given `theta`, so a unit is one level (a
#' run of `dim` consecutive coefficient columns). A structure that ties
#' the levels together through a fixed matrix over them
#' (`cs_levels_coupled()`: `gr(cov=)`, `gr(prec=)`, `car()` and the
#' SPDE) makes the whole block one unit.
#' Smooths, `gp()` and `hsgp()` already come through with `n_levels =
#' 1`, which lands in the same place without a special case.
#'
#' @noRd
re_prior_units <- function(frame) {
  u <- integer(frame[["n_c"]] %||% 0L)
  if (!length(u)) return(u)
  id <- 0L
  for (bk in frame[["re_blocks"]]) {
    cols <- bk[["c_idx"]]
    if (cs_levels_coupled(bk)) {
      id <- id + 1L
      u[cols] <- id
      next
    }
    d <- bk[["dim"]]
    nlev <- length(cols) %/% d
    u[cols] <- id + rep(seq_len(nlev), each = d)
    id <- id + nlev
  }
  u
}

#' Row / coefficient-column incidence of the random effects, unioned
#' over the linear predictors (an `id`-shared block enters more than
#' one of them).
#'
#' @noRd
re_incidence <- function(frame) {
  ii <- integer(0)
  jj <- integer(0)
  for (lp in frame[["linpreds"]]) {
    if (is.null(lp[["Z"]])) next
    tz <- methods::as(methods::as(lp[["Z"]], "generalMatrix"), "TsparseMatrix")
    keep <- tz@x != 0
    ii <- c(ii, tz@i[keep] + 1L)
    jj <- c(jj, tz@j[keep] + 1L)
  }
  list(i = ii, j = jj)
}

#' Every guard that has to hold before a cluster score means anything.
#' Each refusal names the offending structure and the way out, because
#' "the likelihood does not factor" is not something a user can debug
#' from a covariance matrix full of plausible numbers.
#'
#' @noRd
cluster_guard <- function(fit, cl) {
  frame <- fit$frame

  if (isTRUE(fit$REML)) {
    frm_stop("vcov_cluster() needs a maximum-likelihood fit. The ",
             "restricted likelihood integrates the fixed effects out ",
             "jointly, so it does not factor over clusters and a ",
             "per-cluster score is not defined. Refit with REML = FALSE",
             call. = FALSE)
  }
  if (isTRUE(fit$control$profile)) {
    frm_stop("vcov_cluster() cannot use a fit made with ",
             "frmtmb_control(profile = TRUE): the profiled fixed effects ",
             "are integrated jointly, so the objective does not factor ",
             "over clusters. Refit with profile = FALSE", call. = FALSE)
  }
  if (isTRUE(fit$quadrature)) {
    frm_stop("vcov_cluster() cannot use a fit made with quadrature = ",
             "TRUE: the Gauss-Kronrod rule integrates a cluster with no ",
             "data only to quadrature accuracy, so the per-cluster scores ",
             "would not be exact. Refit with quadrature = FALSE (Laplace)",
             call. = FALSE)
  }
  if (!is.null(fit$importance)) {
    frm_stop("vcov_cluster() cannot use a fit made with an importance ",
             "correction: that objective is a sum over GROUPS of ",
             "reweighted integrals, not the per-row sum the cluster ",
             "scores are read off, so a cluster carrying no rows of its ",
             "own still moves it. Refit with importance = 0 (Laplace)",
             call. = FALSE)
  }
  if (!is.null(fit$prior)) {
    frm_stop("vcov_cluster() cannot use a fit made with priors: the ",
             "optimized objective is penalized, the penalty belongs to no ",
             "cluster, and the sandwich is not a covariance for a ",
             "penalized estimator. Refit without priors, or use ",
             "frm_bootstrap()", call. = FALSE)
  }
  if (length(frame[["autocor"]] %||% list())) {
    frm_stop("vcov_cluster() does not support the residual correlation ",
             "term ", frame[["autocor"]][[1L]]$label,
             ": its density is a joint ",
             "one over each group rather than a product over rows (cov = ",
             "TRUE), or each row's mean reads the residuals of earlier ",
             "rows (cov = FALSE), so a cluster's weight does not reach ",
             "all of its contribution. Use frm_bootstrap()",
             call. = FALSE)
  }
  # A structured family says for itself whether a per-cluster score
  # exists. A group-level mixture's does, once every group sits inside
  # one cluster, which is checked at the end of this function; a joint
  # density over a whole sequence has none to offer.
  for (resp_ in fit$spec$responses) {
    structure_gate(fam_structure(resp_$family), "cluster_robust",
                   structure_generic(resp_$family, "vcov_cluster()"))
  }
  if (isTRUE(fit$spec$rescor)) {
    frm_stop("vcov_cluster() does not support rescor = TRUE: the ",
             "responses enter one joint density per row, which the ",
             "cluster weights do not reach. Use frm_bootstrap()",
             call. = FALSE)
  }
  if (length(frame[["subset_rows"]] %||% list())) {
    frm_stop("vcov_cluster() does not support a multivariate fit whose ",
             "responses use subset(): each response has its own rows, and ",
             "`cluster` is given for the rows of the whole data. Use ",
             "frm_bootstrap()", call. = FALSE)
  }
  if (length(frame[["mi_map"]] %||% list()) || !is.null(frame[["me"]])) {
    frm_stop("vcov_cluster() does not support mi() / me() fits: the latent ",
             "values are parameters of the outer problem and their ",
             "contribution belongs to no cluster. Use frm_bootstrap()",
             call. = FALSE)
  }
  if (nlevels(cl) < 2L) {
    frm_stop("`cluster` has ", nlevels(cl), " level(s); a cluster-robust ",
             "covariance needs at least 2", call. = FALSE)
  }

  ci <- as.integer(cl)
  # random effects: every prior-independent unit must live in ONE
  # cluster, or the marginal likelihood does not factor
  u <- re_prior_units(frame)
  if (length(u)) {
    inc <- re_incidence(frame)
    if (length(inc$i)) {
      uu <- u[inc$j]
      pairs <- unique(cbind(uu, ci[inc$i]))
      bad <- pairs[duplicated(pairs[, 1L]), 1L]
      if (length(bad)) {
        blk <- frame[["re_blocks"]][[which(vapply(
          frame[["re_blocks"]],
          function(bk) any(u[bk[["c_idx"]]] %in% bad), TRUE))[1L]]]
        frm_stop("the random effect ", blk[["term_label"]], " crosses ",
                 "`cluster`: one of its levels loads on rows in more than ",
                 "one cluster, so the marginal likelihood does not factor ",
                 "over clusters and a per-cluster score is not defined. ",
                 "Cluster at a factor that every random effect is nested ",
                 "in, or use frm_bootstrap()", call. = FALSE)
      }
    }
  }
  # group-level (latent-class) mixtures: same condition on the
  # mixture's own grouping. A family that does not declare
  # cluster_robust was refused above, so a block still standing here
  # carries a group incidence this can read.
  for (r in names(frame[["blocks"]] %||% list())) {
    mg <- frame_block_of(frame, r)
    tg <- methods::as(methods::as(mg$Gt, "generalMatrix"), "TsparseMatrix")
    pairs <- unique(cbind(tg@i + 1L, ci[tg@j + 1L]))
    if (anyDuplicated(pairs[, 1L]) > 0L) {
      frm_stop("the group-level mixture grouping for response ", r,
               " crosses `cluster`: one mixture group spans more than one ",
               "cluster, so the likelihood does not factor. Cluster at ",
               "that grouping factor or coarser", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Per-cluster score matrix
#'
#' The gradient of each cluster's contribution to the marginal
#' log-likelihood, evaluated at the fitted estimates. This is the
#' quantity `sandwich::estfun()` would return if the marginalized
#' objective had per-observation contributions, which it does not - see
#' [vcov_cluster()] for the argument and the guards.
#'
#' Rows sum to the gradient of the full objective, which is (near) zero
#' at a converged optimum; that identity is the cheapest check that the
#' clustering factor really does split the likelihood.
#'
#' @param object A `frmtmb_fit`.
#' @param cluster The clustering factor: a one-sided formula such as
#'   `~ g` (looked up in the model frame, then in the data the model
#'   was fitted to), a variable name, or a vector with one entry per
#'   fitted row.
#' @return A matrix with one row per cluster level and one column per
#'   estimated outer parameter, named as in `vcov(object, full =
#'   TRUE)`. Signs follow `sandwich::estfun()`: these are scores of the
#'   log-likelihood, not of the objective.
#' @seealso [vcov_cluster()], which sandwiches these between two copies
#'   of the inverse observed information.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200), g = factor(rep(1:25, 8)))
#' dd$y <- rnorm(200, 1 + 0.5 * dd$x + rnorm(25, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = FALSE)
#'
#' S <- cluster_scores(fit, ~ g)
#' dim(S)
#' # the scores add up to the gradient at the optimum
#' max(abs(colSums(S)))
#' @export
cluster_scores <- function(object, cluster) {
  if (!inherits(object, "frmtmb_fit")) {
    frm_stop("cluster_scores() needs a model fitted by frm(), not ",
             arg_desc(object), call. = FALSE)
  }
  cl <- resolve_cluster(object, cluster)
  cluster_guard(object, cl)
  cluster_scores_at(object, cl)
}

#' The score machinery itself, on an already-validated clustering.
#'
#' @noRd
cluster_scores_at <- function(fit, cl) {
  frame <- fit$frame
  gs <- levels(cl)
  # the hook is set on a PRIVATE copy of the frame; the frame stored on
  # the fit never carries it, so no other method can see the extra
  # parameter
  frame[["cluster_w"]] <- lapply(
    stats::setNames(nm = names(frame[["aterm_values"]])),
    function(r) as.integer(cl))
  # the full-shape parameter list: `env$parameters` holds a mapped
  # component in the map's reduced shape, so a partly mapped one (the
  # ordinal families' disc, held at 1 beside a free coefficient) failed
  # MakeADFun with "A map factor length must equal parameter length"
  tpl <- fit$obj$env$parList(fit$opt$par)
  tpl[["clw"]] <- rep(1, length(gs))
  ri <- fit$obj$env$random
  random <- if (length(ri)) unique(names(fit$obj$env$par)[ri])
  obj <- RTMB::MakeADFun(build_objective(frame), tpl, random = random,
                         map = frame[["map"]], silent = TRUE)
  keep <- names(obj$par) != "clw"
  if (!identical(names(obj$par)[keep], names(fit$opt$par))) {
    frm_stop("the cluster-weighted objective did not reproduce the ",
             "fitted parameter vector; this is a bug in frmtmb",
             call. = FALSE)
  }
  p <- obj$par
  p[keep] <- fit$opt$par
  S <- matrix(NA_real_, length(gs), sum(keep),
              dimnames = list(gs, outer_par_names(fit)))
  for (i in seq_along(gs)) {
    p[!keep] <- as.numeric(seq_along(gs) == i)
    # obj$gr is the gradient of the NEGATIVE log-likelihood; estfun
    # convention is the log-likelihood
    S[i, ] <- -obj$gr(p)[keep]
  }
  S
}

#' Cluster-robust (sandwich) covariance
#'
#' The misspecification-robust covariance of the estimates when
#' observations are correlated within clusters in a way the model does
#' not describe. With `B` the inverse observed information of the
#' marginal likelihood at the optimum (that is, `vcov(object, full =
#' TRUE)`) and `s_g` the score of cluster `g`
#' ([cluster_scores()]),
#'
#' \deqn{V = a_{G} \, B \left(\sum_g s_g s_g'\right) B}
#'
#' the Liang-Zeger / White cluster sandwich for an M-estimator. The
#' small-sample factor \eqn{a_G} is the clubSandwich one, with `G`
#' clusters, `N` rows and `p` estimated parameters:
#' `"CR0"` is 1, `"CR1"` is `G / (G - 1)`, `"CR1p"` is `G / (G - p)`,
#' and `"CR1S"` is `G * (N - 1) / ((G - 1) * (N - p))`, the Stata
#' `vce(cluster)` factor.
#'
#' # What it is a covariance for
#'
#' The sandwich is taken over the WHOLE outer parameter vector, so the
#' returned fixed-effect block already carries the cost of estimating
#' the covariance parameters. `clubSandwich::vcovCR()` on an `lmerMod`
#' instead conditions on the variance parameters - it sandwiches only
#' the mean-model bread - so the two agree exactly only when the
#' fixed-effect / covariance-parameter block of the observed
#' information vanishes. To reproduce that conditional form, build it
#' from [cluster_scores()] and the fixed-effect block of `solve(vcov(
#' object, full = TRUE))`.
#'
#' `"CR2"` (Bell-McCaffrey) and `"CR3"` are refused rather than
#' approximated. Both are defined through the hat matrix of a linear
#' (or GLS) model; a Laplace-marginal likelihood with a nonlinear link
#' has no such representation, and nothing in the literature defines
#' the adjustment for one. Use `clubSandwich` on a matched linear model
#' if `"CR2"` is what the analysis needs.
#'
#' # When it is refused
#'
#' The marginal likelihood factors over clusters only when every random
#' effect is nested in (or equal to) the clustering factor, and when
#' every likelihood term is a product over rows. Fits that fail either
#' condition are refused with the reason: a random effect whose level
#' spans two clusters (including crossed effects, `mm()` pooled levels,
#' a global smooth, `gp()`, `car()` and the SPDE), a group-level
#' mixture whose groups span clusters, an `autocor()` residual, a
#' family whose structure does not declare `cluster_robust`,
#' `rescor = TRUE`, `mi()`/`me()`, `REML = TRUE`, `profile = TRUE`,
#' `quadrature = TRUE`, and any fit made with priors. [frm_bootstrap()]
#' is the fallback in every one of those cases.
#'
#' Few clusters make the estimator badly biased whatever the
#' correction; the usual advice is at least 30-50, and `"CR1"` with a
#' `t(G - 1)` reference distribution below that. [confint()] and
#' [hypothesis()] use the `t(G - 1)` reference automatically when they
#' are handed a matrix from this function.
#'
#' @references
#' Liang, K.-Y. and Zeger, S. L. (1986). Longitudinal data analysis
#' using generalized linear models. *Biometrika* 73, 13-22.
#'
#' White, H. (1982). Maximum likelihood estimation of misspecified
#' models. *Econometrica* 50, 1-25.
#'
#' Cameron, A. C. and Miller, D. L. (2015). A practitioner's guide to
#' cluster-robust inference. *Journal of Human Resources* 50, 317-372.
#'
#' Pustejovsky, J. E. and Tipton, E. (2018). Small-sample methods for
#' cluster-robust variance estimation and hypothesis testing in fixed
#' effects models. *Journal of Business & Economic Statistics* 36,
#' 672-683.
#'
#' @param object A `frmtmb_fit`.
#' @param cluster The clustering factor, as in [cluster_scores()].
#' @param type Small-sample correction: `"CR0"` (none), `"CR1"`,
#'   `"CR1p"` or `"CR1S"`, spelled and defined as in `clubSandwich`.
#' @param full If `TRUE`, return the whole outer parameter block
#'   (covariance parameters included), named as in `vcov(full = TRUE)`;
#'   otherwise the fixed-effect block, as [vcov.frmtmb_fit()] returns.
#' @return A covariance matrix, with attributes `type`, `nclusters` and
#'   `df` (the `G - 1` reference degrees of freedom). When the fixed
#'   design is ill-conditioned, the model's covariance is held as a
#'   factor (see [frm_lp_basis()]). Then the sandwich is formed through
#'   that factor, and the matrix also has the attribute `root`, a matrix
#'   `F` with `V = F F'`. [hypothesis()] reads `root` when it is given
#'   the matrix as `vcov`. A product formed from the whole matrix loses
#'   its digits on such a fit, and the function warns once per fit.
#' @seealso [cluster_scores()] for the pieces, [vcov.frmtmb_fit()]
#'   which accepts `cluster` and forwards here, and [frm_bootstrap()]
#'   for the cases this refuses.
#'
#' @srrstats {RE4.6} Cluster-robust covariance of the model parameters,
#'   as an alternative to the model-based `vcov()`: the same inverse
#'   observed information sandwiching the outer product of the
#'   per-cluster scores.
#'
#' @examples
#' set.seed(1)
#' G <- 40
#' dd <- data.frame(g = factor(rep(seq_len(G), each = 6)),
#'                  x = rnorm(G * 6))
#' # cluster-specific error scale the model does not describe
#' dd$y <- 1 + 0.5 * dd$x +
#'   rnorm(G * 6, 0, rep(runif(G, 0.3, 2.5), each = 6))
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
#'            REML = FALSE)
#'
#' sqrt(diag(vcov(fit)))                          # model-based
#' sqrt(diag(vcov_cluster(fit, ~ g, "CR1")))      # cluster-robust
#'
#' # feeds straight into the inference methods, which need the whole
#' # outer parameter vector and pick up the t(G - 1) reference from it
#' confint(fit, parm = "x",
#'         vcov = vcov_cluster(fit, ~ g, "CR1", full = TRUE))
#' @export
vcov_cluster <- function(object, cluster, type = c("CR0", "CR1",
                                                   "CR1p", "CR1S"),
                         full = FALSE) {
  if (!inherits(object, "frmtmb_fit")) {
    frm_stop("vcov_cluster() needs a model fitted by frm(), not ",
             arg_desc(object), call. = FALSE)
  }
  check_flag(full, "full")
  if (is.character(type) && length(type) == 1L &&
      type %in% c("CR2", "CR3", "CR4")) {
    frm_stop("type = '", type, "' is not defined for a marginal-likelihood ",
             "fit. The Bell-McCaffrey family is built from the hat matrix ",
             "of a linear (or GLS) model, which a Laplace-marginal ",
             "likelihood with a nonlinear link does not have, and frmtmb ",
             "will not ship an adjustment with no derivation behind it. ",
             "Use type = 'CR1', or clubSandwich::vcovCR() on a matched ",
             "linear model", call. = FALSE)
  }
  type <- frm_match_arg(type)
  cl <- resolve_cluster(object, cluster)
  cs <- cluster_sandwich_full(object, cl, type)
  V <- cs$V
  Fr <- cs$Fr
  G <- cs$G
  nm <- rownames(V)
  keep <- seq_along(nm)
  if (!full) {
    # by COMPONENT, not by name: vcov() does the same, and two linear
    # predictors can carry the same coefficient name
    cmp <- outer_par_map(object)$comp
    keep <- c(which(cmp == "beta"), which(cmp == "betad"))
    V <- V[keep, keep, drop = FALSE]
    dimnames(V) <- list(estimated_coef_names(object),
                        estimated_coef_names(object))
  }
  if (!is.null(Fr)) {
    Fr <- Fr[keep, , drop = FALSE]
    rownames(Fr) <- rownames(V)
    attr(V, "root") <- Fr
  }
  attr(V, "type") <- type
  attr(V, "nclusters") <- G
  attr(V, "df") <- G - 1
  if (!is.null(Fr)) warn_cov_products(object, "vcov_cluster()")
  V
}

#' The cluster sandwich over the whole outer parameter vector: `V`, with
#' `NaN` in the rows of a parameter the fit lost (as `vcov(full = TRUE)`
#' shows it), and `Vp`, the same sandwich before that mask.
#'
#' The bread is the covariance conditional on the lost directions
#' (`cov_fixed_prop`) where the fit lost a parameter. Formed from the
#' covariance that `vcov()` shows, one `NaN` row spread to every entry
#' of `B S'S B`: on a fit with a `mo()` category no row is in, a
#' one-value `mo()` term or a variance at its boundary, every entry was
#' `NaN`. The cluster scores are 0 along an exactly flat direction, so
#' for any quantity that does not move along one, the sandwich from
#' this bread is the sandwich of the determined parameters.
#'
#' @noRd
cluster_sandwich_full <- function(object, cl, type) {
  cluster_guard(object, cl)
  S <- cluster_scores_at(object, cl)
  # vcov_estimated(), not vcov(): vcov()'s warning is about the matrix
  # it returns, and this one is formed below
  B <- vcov_estimated(object, full = TRUE)
  nm <- outer_par_labels(object)
  if (!identical(dim(B), c(length(nm), length(nm)))) {
    frm_stop("the model-based covariance and the score matrix do not ",
             "line up; this is a bug in frmtmb", call. = FALSE)
  }
  lost <- !is.finite(diag(B))
  Bp <- B
  if (any(lost)) {
    sdr <- tryCatch(sdr_of(object), error = keep_missing(NULL))
    Vprop <- sdr$cov_fixed_prop
    if (!is.null(Vprop) && all(dim(Vprop) == dim(B))) {
      Bp <- as.matrix(Vprop)
    }
  }
  G <- nrow(S)
  N <- object$frame[["n_obs"]]
  p <- length(object$opt$par)
  if (type %in% c("CR1p", "CR1S") && G <= p) {
    frm_warning("type = '", type, "' needs more clusters (", G,
                ") than estimated parameters (", p,
                "); the small-sample factor is not positive", call. = FALSE)
  }
  cf <- cr_adjust(type, G, N, p)
  Fr <- cluster_root(object, S, cf, nrow(B))
  if (is.null(Fr)) {
    Vp <- Bp %*% crossprod(S) %*% Bp
    Vp <- cf * Vp
    Vp <- (Vp + t(Vp)) / 2                 # kill the asymmetric round-off
  } else {
    Vp <- tcrossprod(Fr)
    Fr[lost, ] <- NaN
  }
  dimnames(Vp) <- list(nm, nm)
  # a parameter the check lost has no row in the bread
  V <- Vp
  V[lost, ] <- NaN
  V[, lost] <- NaN
  list(V = V, Vp = Vp, Fr = Fr, G = G, lost = lost)
}

#' A factor `F` of the cluster-robust covariance, `V = F F'`, on a fit
#' whose covariance is held as a factor `root` (se_reparam_cov()):
#' `F = sqrt(cf) root K'` with `K = S root`, so the sandwich
#' `B S'S B` is never formed from the bread `B = root root'`. On a
#' degree-8 raw polynomial with `(1 | g)` the sandwich from the whole
#' bread was off by 2.3 to 12 on the correlation scale, and even the
#' exact bread rounded to double gave 3.5 to 10; through the factor the
#' error was 6.4e-6 to 1.5e-5 (review of lane illse, item 1). With more
#' clusters than parameters `F` is reduced to a square factor by a QR
#' of `F'`. NULL on any other fit, or when `cf` is not positive.
#'
#' @noRd
cluster_root <- function(object, S, cf, np) {
  root <- tryCatch(sdr_of(object)$cov_fixed_root, error = keep_missing(NULL))
  if (is.null(root) || nrow(root) != np || ncol(S) != np ||
        !is.finite(cf) || cf <= 0) {
    return(NULL)
  }
  K <- S %*% root
  Fr <- sqrt(cf) * root %*% t(K)
  if (ncol(Fr) > np) {
    q <- qr(t(Fr))
    Fr <- t(qr.R(q)[, order(q$pivot), drop = FALSE])
  }
  if (all(is.finite(Fr))) Fr
}

#' Resolves a user-supplied `vcov` argument (a matrix, or a function of
#' the fit) against the outer parameter names, and reports the
#' reference degrees of freedom that came with it.
#'
#' @noRd
resolve_vcov_arg <- function(fit, V, what) {
  if (is.function(V)) V <- V(fit)
  nm <- outer_par_names(fit)
  if (!is.matrix(V) || nrow(V) != ncol(V)) {
    frm_stop(what, "(vcov = ) needs a square matrix (or a function ",
             "returning one)", call. = FALSE)
  }
  if (nrow(V) != length(nm)) {
    hint <- if (nrow(V) == length(estimated_coef_names(fit))) {
      " That is the fixed-effect block; pass full = TRUE."
    } else ""
    frm_stop(what, "(vcov = ) needs a ", length(nm), " x ", length(nm),
             " matrix over the whole outer parameter vector, as ",
             "vcov(object, full = TRUE) and vcov_cluster(full = TRUE) ",
             "return; got ", nrow(V), " x ", ncol(V), ".", hint,
             call. = FALSE)
  }
  # a factor that came with the matrix (vcov_cluster() on a fit whose
  # covariance is held as one), for the routes that form products
  root <- attr(V, "root")
  if (!is.matrix(root) || nrow(root) != nrow(V)) root <- NULL
  list(V = V, df = attr(V, "df"), root = root)
}
