# Custom-family checking, and the downstream-package glue: emmeans,
# marginaleffects, lme4::getME.
#
# The sampling surface that used to sit between them left for the
# frmtmb.sample package (dev/draws-extraction.md). The prior and bound
# machinery the FIT route reaches - frm(), par_template() and
# frm_simulate() - was already in R/priors.R rather than here, which is
# what let that cut miss this file's prior code entirely.

#' Check a custom family's log-density for AD safety
#'
#' Tapes the family's `lpdf` on test values and compares the AD gradient
#' against central finite differences. A mismatch usually means the lpdf
#' uses operations the tape cannot see (base `matrix()`/`c()` on
#' advectors, branching on parameter values, `min`/`max`, clamping).
#'
#' @param family A `frmtmb_family` (from [frmtmb_family()] /
#'   [custom_family()]).
#' @param y A response vector of test data.
#' @param dpars Named list of numeric test values, one entry per dpar
#'   (each of length 1 or `length(y)`).
#' @param aterms Named list of addition-term values (e.g. `trials`).
#' @param tol Maximum relative gradient error.
#' @return Invisibly `TRUE`; signals an error on failure.
#' @examples
#' set.seed(1)
#' y <- rpois(50, 3)
#'
#' # a hand-written poisson: check it before fitting anything with it
#' ok <- custom_family(
#'   "my_poisson", dpars = "mu", links = list(mu = "log"),
#'   lpdf = function(y, dpars, aterms) {
#'     y * log(dpars$mu) - dpars$mu - lgamma(y + 1)
#'   },
#'   type = "discrete"
#' )
#' check_custom_family(ok, y = y, dpars = list(mu = rep(2.5, 50)))
#'
#' # base matrix() strips the advector class, so the tape sees constants
#' # and the gradient is silently wrong. The check catches it.
#' bad <- custom_family(
#'   "bad", dpars = "mu", links = list(mu = "log"),
#'   lpdf = function(y, dpars, aterms) {
#'     m <- matrix(dpars$mu, ncol = 1)
#'     y * log(m[, 1]) - m[, 1] - lgamma(y + 1)
#'   },
#'   type = "discrete"
#' )
#' try(check_custom_family(bad, y = y, dpars = list(mu = rep(2.5, 50))))
#' @export
check_custom_family <- function(family, y, dpars, aterms = list(),
                                tol = 1e-4) {
  stopifnot(inherits(family, "frmtmb_family"))
  check_positive(tol, "tol")
  if (!setequal(names(dpars), family$dpars)) {
    stop("`dpars` must supply test values for exactly: ",
         paste(family$dpars, collapse = ", "), call. = FALSE)
  }
  f <- function(p) -sum(family$lpdf(y, p, aterms))
  v0 <- f(lapply(dpars, as.numeric))
  if (!is.finite(v0)) {
    stop("lpdf is not finite at the test values", call. = FALSE)
  }
  obj <- tryCatch(
    RTMB::MakeADFun(f, dpars, silent = TRUE),
    error = function(e) {
      stop("Failed to tape the lpdf: ", conditionMessage(e),
           ". Typical cause: base matrix()/c() stripping the advector ",
           "class, or branching on parameter values", call. = FALSE)
    }
  )
  if (abs(obj$fn(obj$par) - v0) > 1e-8 * max(1, abs(v0))) {
    stop("Taped lpdf disagrees with its plain-numeric value (",
         format(obj$fn(obj$par)), " vs ", format(v0), "): the lpdf uses ",
         "operations that behave differently on the AD tape (base ",
         "matrix()/c() on advectors are the usual culprits)",
         call. = FALSE)
  }
  g <- as.vector(obj$gr(obj$par))
  p0 <- obj$par
  h <- 1e-6 * pmax(abs(p0), 1)
  fd <- vapply(seq_along(p0), function(i) {
    pp <- p0; pp[i] <- pp[i] + h[i]
    pm <- p0; pm[i] <- pm[i] - h[i]
    (obj$fn(pp) - obj$fn(pm)) / (2 * h[i])
  }, numeric(1))
  rel <- abs(g - fd) / pmax(abs(fd), 1)
  if (any(rel > tol)) {
    stop("AD gradient disagrees with finite differences (max relative ",
         "error ", format(max(rel), digits = 3), ")", call. = FALSE)
  }
  invisible(TRUE)
}

#' emmeans support: registered in .onLoad, only for the parametric fixed
#' part of the mu linear predictor of univariate, non-nl fits.
#'
#' An ordinal fit therefore lands on the LATENT scale, which is
#' emmeans's own `mode = "latent"` convention for `clm`-like models: the
#' K-1 thresholds are not coefficients of this design, so the basis has
#' no intercept and the marginal means are the latent predictor without
#' a threshold offset. Contrasts are unaffected by that offset;
#' category probabilities are a different question, answered by
#' `predict(type = "response")` and `conditional_effects()`.
#'
#' @noRd
emm_mu_linpred <- function(object) {
  if (length(object$spec$responses) > 1) {
    stop("emmeans support is univariate-only for now", call. = FALSE)
  }
  rspec <- object$spec$responses[[1]]
  lp <- object$frame$linpreds[[linpred_key(rspec$resp_name, "mu")]]
  if (is.null(lp) || !is.null(lp$nl_body)) {
    stop("emmeans support needs a linear mu predictor", call. = FALSE)
  }
  lp
}

# marginaleffects support: the four extension generics plus the
# class-whitelist option set in .onLoad. Predictions under set_coef are
# conditional on the estimated random-effect modes (the glmmTMB/lmer
# convention for delta-method slopes).

#' @exportS3Method marginaleffects::get_coef
get_coef.frmtmb_fit <- function(model, ...) {
  est <- model$estimates
  bd <- est$betad
  if (length(model$frame$betad_fixed_idx)) {
    bd <- bd[-model$frame$betad_fixed_idx]
  }
  stats::setNames(c(est$beta, bd), estimated_coef_names(model))
}

#' @exportS3Method marginaleffects::set_coef
set_coef.frmtmb_fit <- function(model, coefs, ...) {
  tpl <- model$frame$par_template
  nb <- length(tpl$beta)
  model$estimates$beta[] <- coefs[seq_len(nb)]
  if (!is.null(tpl$betad)) {
    keep <- setdiff(seq_along(tpl$betad), model$frame$betad_fixed_idx)
    model$estimates$betad[keep] <- coefs[nb + seq_along(keep)]
  }
  model
}

#' @exportS3Method marginaleffects::get_vcov
get_vcov.frmtmb_fit <- function(model, ...) {
  vcov(model)
}

#' @exportS3Method marginaleffects::get_predict
get_predict.frmtmb_fit <- function(model, newdata, type = "response",
                                   ...) {
  type <- if (identical(type, "link")) "link" else "response"
  p <- predict(model, newdata = newdata, type = type)
  if (is.matrix(p)) {
    # A categorical outcome predicts a DISTRIBUTION per row (an ordinal
    # family's K category probabilities, a multinomial's D cell means),
    # so one newdata row is several predictions. marginaleffects keys
    # those with a `group` column and repeats the rowid; without it the
    # flattened matrix was handed back as n * K unrelated rows numbered
    # 1..nK, which silently misaligns every downstream contrast.
    g <- colnames(p) %||% as.character(seq_len(ncol(p)))
    return(data.frame(
      rowid = rep(seq_len(nrow(p)), times = ncol(p)),
      group = rep(g, each = nrow(p)),
      estimate = as.numeric(p)
    ))
  }
  data.frame(rowid = seq_along(p), estimate = as.numeric(p))
}

#' @exportS3Method emmeans::recover_data
recover_data.frmtmb_fit <- function(object, ..., data = NULL) {
  lp <- emm_mu_linpred(object)
  emmeans::recover_data(object$call,
                        stats::delete.response(lp$terms),
                        na.action = NULL,
                        data = data %||% model.frame(object), ...)
}

#' @exportS3Method emmeans::emm_basis
emm_basis.frmtmb_fit <- function(object, trms, xlev, grid, ...) {
  lp <- emm_mu_linpred(object)
  m <- stats::model.frame(trms, grid, na.action = stats::na.pass,
                          xlev = xlev)
  X <- stats::model.matrix(trms, m, contrasts.arg = lp$contrasts)
  # The fitted design is not always model.matrix()'s: an ordinal family
  # drops the intercept (the K-1 thresholds take its place), and a
  # rank-deficient fit drops aliased columns. Select the fitted columns
  # by name, or the basis is not conformable with bhat and emmeans fails
  # with "Non-conformable elements in reference grid".
  pn <- lp$param_colnames[seq_len(lp$n_param_cols)]
  if (!identical(colnames(X), pn)) {
    keep <- match(pn, colnames(X))
    if (anyNA(keep)) {
      stop("emmeans support cannot rebuild the fitted design: column(s) ",
           paste(pn[is.na(keep)], collapse = ", "),
           " are missing from the reference grid", call. = FALSE)
    }
    X <- X[, keep, drop = FALSE]
  }
  idx <- lp$idx[seq_len(lp$n_param_cols)]
  bhat <- object$estimates[[lp$par]][idx]
  V <- vcov(object)[idx, idx, drop = FALSE]
  list(X = X, bhat = bhat, nbasis = matrix(NA), V = V,
       dffun = function(k, dfargs) Inf, dfargs = list(), misc = list())
}

# --- lme4::getME ------------------------------------------------------
# Only the pieces that mean the same thing here. The vocabulary stays
# small on purpose: a name that would have to be faked (Lambdat, u, the
# lme4 sparse-Cholesky machinery) is worse than a name that errors,
# because downstream code cannot tell a wrong answer from a right one.

frmtmb_getME_vocab <- c("X", "Z", "Zt", "beta", "fixef", "b", "theta",
                        "lower", "sigma", "flist", "n_rtrms",
                        "n_rfacs")

#' Blocks that carry a real grouping factor. Smooths and Gaussian
#' processes are stored as random-effect blocks too, but their "levels"
#' are basis functions, so they have no factor to report.
#'
#' @noRd
getME_group_blocks <- function(object) {
  Filter(function(bk) !bk$covstruct %in% c("smooth", "gp", "hsgp") &&
           !is.null(bk$components[[1L]]$bar),
         object$frame$re_blocks)
}

#' Labels for the random-effect coefficient vector, which is also the
#' column order of every Z: level-major within a block, one entry per
#' (level, term coefficient), the way lme4 labels Zt rows.
#'
#' @noRd
re_coef_labels <- function(frame) {
  if (!length(frame$re_blocks)) return(character(0))
  unlist(lapply(frame$re_blocks, function(bk) {
    lv <- bk$levels %||% as.character(seq_len(bk$n_levels))
    paste0(rep(lv, each = bk$dim), ".", rep(bk$cnms, bk$n_levels))
  }), use.names = FALSE)
}

#' The grouping factors as they were at fit time, rebuilt from the
#' stored model frame the same way predict() rebuilds them for newdata.
#'
#' @noRd
getME_flist <- function(object) {
  out <- list()
  for (bk in getME_group_blocks(object)) {
    comp <- bk$components[[1L]]
    # a multi-membership row belongs to several levels at once, so there
    # is no per-observation grouping factor to report; lme4's flist has
    # no representation for one, and returning the first member would
    # be a wrong answer rather than a missing one
    if (!is.null(comp$mm)) next
    lp <- object$frame$linpreds[[comp$lp_key]]
    env <- object$spec$responses[[lp$resp]]$formula_env
    gv <- tryCatch(
      as.character(eval(comp$bar[[3L]], object$frame$data_frame, env)),
      error = function(e) NULL
    )
    if (length(gv) != object$frame$n_obs) {
      stop("getME(\"flist\"): cannot rebuild the grouping factor for `",
           bk$term_label, "`", call. = FALSE)
    }
    nm <- bk$group_name
    if (is.null(out[[nm]])) out[[nm]] <- factor(gv, levels = bk$levels)
  }
  out
}

#' Extract components of a fit, lme4 style
#'
#' A small [lme4::getME()] vocabulary, for downstream code written
#' against merMod objects. Registered on lme4's generic, so call it as
#' `lme4::getME(fit, "X")` or load lme4 first.
#'
#' Supported names:
#' \describe{
#'   \item{`"X"`}{Fixed-effect design matrix of the `mu` predictor.}
#'   \item{`"Z"`, `"Zt"`}{The sparse random-effect design of the `mu`
#'     predictor and its transpose. Columns (rows of `Zt`) span the
#'     whole random-effect coefficient vector, so a block belonging to
#'     another distributional parameter contributes zero columns here.}
#'   \item{`"beta"`, `"fixef"`}{The primary (`mu`-family) fixed-effect
#'     coefficients, named. Coefficients of auxiliary distributional
#'     parameters are a separate vector; use [fixef()] for all of them.}
#'   \item{`"b"`}{Conditional modes in coefficient space, aligned with
#'     the columns of `Z`. Reduced-rank (`rr()`) blocks are expanded
#'     through their loadings, so this is not the internal parameter
#'     vector.}
#'   \item{`"theta"`}{Covariance parameters on the internal
#'     (unconstrained) scale, as in `confint()`. These are not lme4's
#'     relative-covariance-factor entries.}
#'   \item{`"lower"`}{Lower bounds on `theta`. The internal
#'     parameterization is unbounded, so this is a vector of `-Inf`, not
#'     lme4's mixture of `0` and `-Inf`. Code that tests
#'     `theta == lower` to detect a singular fit will never fire; use
#'     [diagnose()] or [VarCorr()] instead.}
#'   \item{`"sigma"`}{Residual standard deviation
#'     ([sigma.frmtmb_fit()]).}
#'   \item{`"flist"`}{The grouping factors, one per distinct grouping
#'     variable. Smooth and Gaussian-process blocks are excluded: their
#'     levels are basis functions, not groups. There is no `"assign"`
#'     attribute.}
#'   \item{`"n_rtrms"`, `"n_rfacs"`}{Number of random-effect terms and
#'     of distinct grouping factors.}
#' }
#'
#' Multivariate fits have one design per response, so `"X"`, `"Z"` and
#' `"Zt"` need `resp`; the other names answer without it.
#'
#' @param object A `frmtmb_fit`.
#' @param name One or more names from the vocabulary above. A vector
#'   returns a named list.
#' @param resp Response name, for the design extractors on a
#'   multivariate fit.
#' @param ... Unused.
#' @return The requested component, or a named list when `name` names
#'   several.
#'
#' @srrstats {RE4.13} Predictor variables and their metadata are
#'   retrievable from the fitted object. `getME()` exposes the
#'   fixed-effect design `X`, the random-effect design `Z` carrying the
#'   row names of the input data, the grouping-factor structure, and the
#'   parameter vectors, using lme4's vocabulary; `model.matrix()` returns
#'   the design and `model.frame()` the stored model frame with its row
#'   names. The frame also keeps the `terms`, `xlevels`, and `contrasts`
#'   of each linear predictor, frozen at fit time and reapplied to
#'   `newdata`, plus the frozen bases of data-dependent terms such as
#'   `poly()` and `scale()`. A name outside the vocabulary errors and
#'   lists the accepted names.
#'
#' @examples
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   set.seed(1)
#'   dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#'   dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#'   fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#'   # the designs, for downstream code written against merMod objects
#'   dim(lme4::getME(fit, "X"))
#'   dim(lme4::getME(fit, "Zt"))
#'
#'   # a vector of names returns a named list
#'   str(lme4::getME(fit, c("n_rtrms", "n_rfacs", "sigma")))
#'
#'   # the conditional modes in coefficient space, aligned with Z
#'   head(lme4::getME(fit, "b"))
#'   # note: "lower" is all -Inf here, because the internal covariance
#'   # parameterization is unbounded. Use diagnose() to spot a singular
#'   # fit, not theta == lower.
#'   lme4::getME(fit, "lower")
#' }
#' @exportS3Method lme4::getME
getME.frmtmb_fit <- function(object, name, resp = NULL, ...) {
  if (missing(name) || !is.character(name) || !length(name)) {
    stop("getME() needs one or more names: ",
         paste(frmtmb_getME_vocab, collapse = ", "), call. = FALSE)
  }
  bad <- setdiff(name, frmtmb_getME_vocab)
  if (length(bad)) {
    stop("getME(): unknown name(s) ", paste(bad, collapse = ", "),
         ". Supported: ", paste(frmtmb_getME_vocab, collapse = ", "),
         call. = FALSE)
  }
  if (length(name) > 1L) {
    out <- lapply(name, function(nm) getME.frmtmb_fit(object, nm, resp))
    return(stats::setNames(out, name))
  }
  # find_linpred() raises the "disambiguate with resp =" error for a
  # multivariate fit, which is exactly the guard the designs need
  mu_lp <- function() find_linpred(object, resp, "mu")
  mu_Z <- function() {
    Z <- mu_lp()$Z
    if (is.null(Z)) {
      Z <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                x = numeric(0),
                                dims = c(object$frame$n_obs,
                                         object$frame$n_c %||% 0L))
    }
    dimnames(Z) <- list(rownames(object$frame$data_frame),
                        re_coef_labels(object$frame))
    Z
  }
  theta <- object$estimates$theta %||% numeric(0)
  names(theta) <- if (length(theta)) {
    paste0("theta_", seq_along(theta))
  }
  switch(
    name,
    X = mu_lp()$X,
    Z = mu_Z(),
    Zt = Matrix::t(mu_Z()),
    beta = ,
    fixef = object$estimates$beta,
    b = {
      bv <- coef_b(object) %||% numeric(0)
      stats::setNames(as.numeric(bv), re_coef_labels(object$frame))
    },
    theta = theta,
    lower = stats::setNames(rep(-Inf, length(theta)), names(theta)),
    sigma = sigma(object),
    flist = getME_flist(object),
    n_rtrms = length(getME_group_blocks(object)),
    n_rfacs = length(getME_flist(object))
  )
}
