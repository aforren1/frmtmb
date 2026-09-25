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
#' It differentiates with respect to the DPAR VALUES, and it supplies no
#' linear predictors, so a density that reads the linear-predictor scale
#' through [frmtmb-robust-dpars] is checked here on its fallback path
#' only. That path is the plain arithmetic, and the plain arithmetic is
#' what saturates. This function cannot reach the other path: `dpars`
#' must name exactly the family's own distributional parameters, so a
#' `.eta_<dpar>` entry is refused rather than taped. Call the lpdf
#' yourself on a list that carries one; [frmtmb-robust-dpars] shows how.
#'
#' @param family A `frmtmb_family` (from [frmtmb_family()] /
#'   [custom_family()]).
#' @param y A response vector of test data.
#' @param dpars Named list of numeric test values, one entry per dpar
#'   (each of length 1 or `length(y)`).
#' @param aterms Named list of addition-term values (e.g. `trials`).
#' @param tol Maximum relative gradient error.
#' @return Invisibly `TRUE`; signals an error on failure.
#' @seealso [frmtmb_family()] for the family this checks, and
#'   [frmtmb-robust-dpars] for the path this check cannot reach
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
  if (!inherits(family, "frmtmb_family")) {
    frm_stop("check_custom_family() needs a family made by ",
             "custom_family() or frmtmb_family(), not ", arg_desc(family),
             call. = FALSE)
  }
  check_positive(tol, "tol")
  if (!setequal(names(dpars), family[["dpars"]])) {
    frm_stop("`dpars` must supply test values for exactly: ",
             paste(family[["dpars"]], collapse = ", "), call. = FALSE)
  }
  f <- function(p) -sum(family[["lpdf"]](y, p, aterms))
  v0 <- f(lapply(dpars, as.numeric))
  if (!is.finite(v0)) {
    frm_stop("lpdf is not finite at the test values", call. = FALSE)
  }
  obj <- tryCatch(
    RTMB::MakeADFun(f, dpars, silent = TRUE),
    error = function(e) {
      frm_stop("Failed to tape the lpdf: ", conditionMessage(e),
               ". Typical cause: base matrix()/c() stripping the advector ",
               "class, or branching on parameter values", call. = FALSE)
    }
  )
  if (abs(obj$fn(obj$par) - v0) > 1e-8 * max(1, abs(v0))) {
    frm_stop("Taped lpdf disagrees with its plain-numeric value (",
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
    frm_stop("AD gradient disagrees with finite differences (max relative ",
             "error ", format(max(rel), digits = 3), ")", call. = FALSE)
  }
  invisible(TRUE)
}

# marginaleffects support: the four extension generics plus the
# class-whitelist option set in .onLoad. Predictions under set_coef are
# conditional on the estimated random-effect modes (the glmmTMB/lmer
# convention for delta-method slopes).

#' The parameter vector every delta-method seam perturbs, and its
#' names.
#'
#' It is every estimated coefficient plus the parameters an ordinal fit
#' keeps outside `beta` and `betad`: the thresholds and the `cs()`
#' coefficients, on their INTERNAL scale and under `confint()`'s names
#' (`tau_raw_1`, `bcs1_1`).
#'
#' They have to be here. marginaleffects builds its Jacobian by
#' perturbing `get_coef()` one entry at a time and re-predicting, so a
#' parameter the vector leaves out contributes nothing to the standard
#' error, and an ordinal prediction depends on the thresholds. Measured
#' on `bf(ord ~ x) + cumulative()`, `avg_slopes()` reported 0.00029 for
#' the middle category where `MASS::polr` reports 0.02016, a factor of
#' 69, on a point estimate that agreed to five decimals
#' (`dev/shapes-rev-ordse.R`).
#'
#' @noRd
interop_coef_names <- function(model) {
  nm <- estimated_coef_names(model)
  tpl <- model$frame[["par_template"]]
  for (cp in ord_extra_comps(model)) {
    v <- names(tpl[[cp]])
    if (is.null(v)) v <- paste0(cp, "_", seq_along(tpl[[cp]]))
    nm <- c(nm, v)
  }
  nm
}

#' @noRd
interop_coef_vector <- function(model) {
  est <- model$estimates
  bd <- est[["betad"]]
  if (length(model$frame[["betad_fixed_idx"]])) {
    bd <- bd[-model$frame[["betad_fixed_idx"]]]
  }
  v <- c(est[["beta"]], bd)
  for (cp in ord_extra_comps(model)) v <- c(v, as.numeric(est[[cp]]))
  stats::setNames(unname(v), interop_coef_names(model))
}

#' The covariance of `interop_coef_vector()`.
#'
#' `vcov_estimated()` when there is nothing outside the coefficients,
#' so the matrix every other fit gets is unchanged, bit for bit.
#' Otherwise the joint matrix `hypothesis()` already assembles, cut
#' down to these rows and put in this order.
#'
#' @noRd
interop_vcov <- function(model) {
  cps <- ord_extra_comps(model)
  V0 <- vcov_estimated(model)
  if (!length(cps)) return(V0)
  pc <- tryCatch(suppressWarnings(hyp_par_cov(model)),
                 error = function(e) NULL)
  nm <- interop_coef_names(model)
  if (is.null(pc)) return(V0)
  ord <- unlist(lapply(c("beta", "betad", cps),
                       function(cp) which(pc$comp == cp)))
  if (length(ord) != length(nm)) return(V0)
  V <- as.matrix(pc$V[ord, ord, drop = FALSE])
  dimnames(V) <- list(nm, nm)
  V
}

#' @exportS3Method marginaleffects::get_coef
get_coef.frmtmb_fit <- function(model, ...) {
  # Exempt: marginaleffects reaches this through its numderiv machinery
  # with `variables`, `numderiv` and `internal_call` in the dots, and
  # NEITHER static call-site count saw it, because the call is assembled
  # rather than written out. Guarding it broke avg_slopes(), slopes(),
  # avg_comparisons() and hypotheses(); dev/argspell-interop-log.txt has
  # the run. A count of zero written call sites is not evidence here.
  interop_coef_vector(model)
}

#' @exportS3Method marginaleffects::set_coef
set_coef.frmtmb_fit <- function(model, coefs, ...) {
  # Exempt: marginaleffects reaches this through its numderiv machinery
  # with `variables`, `numderiv` and `internal_call` in the dots, and
  # NEITHER static call-site count saw it, because the call is assembled
  # rather than written out. Guarding it broke avg_slopes(), slopes(),
  # avg_comparisons() and hypotheses(); dev/argspell-interop-log.txt has
  # the run. A count of zero written call sites is not evidence here.
  tpl <- model$frame[["par_template"]]
  nb <- length(tpl[["beta"]])
  model$estimates[["beta"]][] <- coefs[seq_len(nb)]
  off <- nb
  if (!is.null(tpl[["betad"]])) {
    keep <- setdiff(seq_along(tpl[["betad"]]), model$frame[["betad_fixed_idx"]])
    model$estimates[["betad"]][keep] <- coefs[nb + seq_along(keep)]
    off <- off + length(keep)
  }
  # the ordinal extras, in interop_coef_vector() order; without these
  # a perturbation of a threshold never reached the prediction
  for (cp in ord_extra_comps(model)) {
    k <- length(model$estimates[[cp]])
    model$estimates[[cp]][] <- coefs[off + seq_len(k)]
    off <- off + k
  }
  model
}

#' @exportS3Method marginaleffects::get_vcov
get_vcov.frmtmb_fit <- function(model, ...) {
  # NOT vcov(model): that is brms's population-level block since item
  # 2.6f, and it has neither the rows nor the names of get_coef(), which
  # is the vector marginaleffects perturbs. A covariance one row short
  # of the coefficient vector is the wrong matrix, not a differently
  # named one; the standard errors of avg_slopes() moved measurably
  # when the two disagreed (dev/shapes-log/interop-lane.txt).
  interop_vcov(model)
}

#' @exportS3Method marginaleffects::get_predict
get_predict.frmtmb_fit <- function(model, newdata, type = "response",
                                   ...) {
  type <- if (identical(type, "link")) "link" else "response"
  p <- frm_linpred(model, newdata = newdata, type = type)
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

# --- emmeans ----------------------------------------------------------
# Registered in .onLoad. The user-facing account is the frmtmb-emmeans
# help topic below; the notes here are about how the pieces fit.

#' emmeans support
#'
#' frmtmb registers [emmeans::emmeans()] support for `frmtmb_fit`
#' objects. The arguments that select what is averaged are brms's:
#' `dpar`, `nlpar`, `resp`, `epred` and `re_formula`, with brms's
#' defaults. Give them to `emmeans()` or `ref_grid()` directly, and
#' emmeans passes them on.
#'
#' @section What is averaged:
#' \describe{
#'   \item{default}{The linear predictor of `mu`, on the link scale,
#'     at the population level (`re_formula = NA`).}
#'   \item{`dpar = "sigma"`}{The linear predictor of that
#'     distributional parameter, on its own link scale.}
#'   \item{`nlpar = "a"`}{The linear predictor of that nonlinear
#'     parameter. It is linear in its own coefficients, so the basis is
#'     the usual design, coefficients and covariance.}
#'   \item{`resp = "y1"`}{One response of a multivariate fit.}
#'   \item{no `resp` on a multivariate fit}{Every response, stacked as
#'     the multivariate factor `rep.meas` whose levels are the response
#'     names. This is brms's layout. Average over it, condition on it
#'     (`by = "rep.meas"`) or compare its levels with `contrast()`.}
#'   \item{`epred = TRUE`}{The expected value of the response, as
#'     [fitted()] reports it, on the response scale.}
#'   \item{`re_formula = NULL`}{Keeps the group-level effects. The
#'     grouping factors join the reference grid, so the marginal mean
#'     averages over their observed levels, as in brms.}
#' }
#'
#' @section How the uncertainty is computed:
#' brms computes each marginal mean from the posterior draws. frmtmb uses
#' the Wald covariance of the estimates instead. For a linear predictor
#' with only parametric terms, the basis is emmeans's usual one: the
#' design at the reference grid, the coefficients and their
#' covariance. Otherwise the basis has brms's shape: one column per
#' grid point, the predictions at the grid as the estimates, and the
#' covariance `J V J'`, where `V` is the joint covariance of the
#' estimates and `J` is the Jacobian of the grid predictions that
#' [frm_lp_basis()] computes. That route is used for a nonlinear
#' predictor, for `epred = TRUE`, for `re_formula` other than `NA`,
#' and for a predictor with a `s()`, `t2()`, `gp()`, `mo()` or `mi()`
#' term. For a linear predictor it is exact. For a nonlinear
#' predictor or the expected response it is the delta method.
#'
#' @section Scales and averaging:
#' emmeans averages the grid estimates with linear weights. For the
#' default route that is the average of linear predictors, as for any
#' generalized linear model. With `epred = TRUE` it is the average of
#' the predicted means, which is how brms marginalizes the expected
#' response. `type = "response"` applies the inverse link of the
#' selected predictor once, after the average. With `epred = TRUE`
#' there is no link to apply, because the estimates are already on
#' the response scale, so `type = "response"` changes nothing. A
#' nonlinear `mu` is reported on its link scale, where the body's value
#' lives; `type = "response"` applies the family's inverse link to it.
#' A response transformation such as `log(y)` is detected for `mu` and
#' for `epred = TRUE` only, since it does not apply to a
#' distributional or nonlinear parameter.
#'
#' An ordinal fit is averaged on its latent scale, which is emmeans's
#' `mode = "latent"` convention for `clm`-like models. The thresholds
#' are not in the basis, so the means carry no threshold offset and
#' `type = "response"` does not transform them. Contrasts are not
#' affected by the offset.
#'
#' @section Refusals:
#' A refusal names its reason. emmeans hides an error that
#' `recover_data()` raises behind "Perhaps a 'data' or 'params'
#' argument is needed", so frmtmb returns the reason as text instead,
#' which emmeans reports as the error. Refused:
#' \itemize{
#'   \item `dpar` and `nlpar` together, and a name the model does not
#'     have (brms refuses both).
#'   \item `epred = TRUE` for an ordinal or categorical family, whose
#'     expected response is a distribution over categories. Use the
#'     latent predictor, or [frm_linpred()] with `type = "response"`.
#'   \item No `resp` on a multivariate fit whose responses have
#'     different links, since the stacked predictors would be on
#'     different scales. brms refuses this too. Name one response, or
#'     use `epred = TRUE`.
#'   \item An exact `gp()` term predicted at a position the fit did not
#'     see. Its kriging variance has no covariance between two grid
#'     points here. Use an approximate `gp(..., k = )`, or hold the
#'     covariate at an observed value with the `at` argument.
#' }
#'
#' @section Divergence from brms:
#' For a nonlinear `mu` without `nlpar`, the brms reference grid holds
#' only the covariates of the nonlinear body. frmtmb also puts in the
#' covariates of the nonlinear parameters' own formulas, so that
#' emmeans averages over them rather than holding them at an unstated
#' value.
#' @name frmtmb-emmeans
#' @return This page documents the emmeans methods and returns nothing.
#'   `emmeans()` on a `frmtmb_fit` returns an `emmGrid` object.
#' @seealso [frm_lp_basis()] for the Jacobian behind the delta method,
#'   and [fitted()] for the expected response.
#' @examples
#' \donttest{
#' if (requireNamespace("emmeans", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(f = factor(rep(c("a", "b", "c"), 40)),
#'                   x = runif(120))
#'   d$y <- 1 + as.numeric(d$f) + 2 * d$x + rnorm(120, 0, 0.3)
#'   d$z <- as.numeric(d$f) + rnorm(120)
#'
#'   # a nonlinear model: one of its parameters, then the whole mean
#'   fnl <- frm(bf(y ~ a + b * x, a ~ f, b ~ 1, nl = TRUE), data = d)
#'   emmeans::emmeans(fnl, "f", nlpar = "a")
#'   emmeans::emmeans(fnl, "f")
#'   pairs(emmeans::emmeans(fnl, "f", epred = TRUE))
#'
#'   # a multivariate model: one response, then both as rep.meas
#'   fmv <- frm(mvbf(bf(y ~ f), bf(z ~ f)), data = d)
#'   emmeans::emmeans(fmv, "f", resp = "y")
#'   emmeans::emmeans(fmv, ~ f | rep.meas)
#' }
#' }
NULL

#' Resolve brms's emmeans arguments against a fit.
#'
#' Called twice per grid, by `recover_data()` and by `emm_basis()`, so
#' both see the same predictor. `route` is `"design"` when every
#' selected predictor is linear with parametric terms only and the
#' prediction is at the population level; there the design basis is
#' exact and is what every fit got before these arguments existed.
#' Anything else takes the `"grid"` route through `frm_lp_basis()`.
#'
#' @noRd
emm_target <- function(object, resp = NULL, dpar = NULL, nlpar = NULL,
                       re_formula = NA, epred = FALSE, warn = FALSE) {
  if (identical(dpar, "mean")) {
    # brms's deprecated spelling of epred = TRUE, honored the same way
    if (warn) {
      frm_warning("dpar = 'mean' is deprecated. Please use epred = TRUE ",
                  "instead.", call. = FALSE)
    }
    epred <- TRUE
    dpar <- NULL
  }
  check_flag(epred, "epred")
  check_re_form(re_formula)
  emm_check_name(dpar, "dpar")
  emm_check_name(nlpar, "nlpar")
  if (!is.null(dpar) && !is.null(nlpar)) {
    frm_stop("'dpar' and 'nlpar' cannot be specified at the same time.",
             call. = FALSE)
  }
  all_resp <- names(object$spec$responses)
  if (is.null(resp)) {
    resp <- all_resp
  } else if (!is.character(resp) || !length(resp) || anyNA(resp) ||
             !all(resp %in% all_resp)) {
    frm_stop("Invalid argument 'resp'. Valid response variables are: ",
             paste0("'", all_resp, "'", collapse = ", "), call. = FALSE)
  }
  resp <- unique(resp)
  # re_resolve() refuses a term that matches no group-level term of the
  # fit, rather than reading a misspelled factor as "drop them all"
  use_re <- re_form_keeps(
    re_resolve(object, re_formula, "emmeans")$re_formula)
  targets <- lapply(resp, emm_target_one, object = object, dpar = dpar,
                    nlpar = nlpar, epred = epred)
  names(targets) <- resp
  design <- !epred && !use_re &&
    all(vapply(targets, function(t) emm_is_design(t$lp), NA))
  misc <- list()
  if (!epred) {
    links <- vapply(targets, `[[`, "", "link")
    if (length(unique(links)) > 1L) {
      frm_stop("emmeans without resp = stacks the responses of a ",
               "multivariate fit on one scale, and these have different ",
               "links (", paste0(resp, ": ", links, collapse = ", "),
               "). Name one response with resp =, or use epred = TRUE ",
               "for the response scale", call. = FALSE)
    }
    misc <- targets[[1L]]$misc
  }
  if (length(resp) > 1L) misc$ylevs <- list(rep.meas = resp)
  list(resp = resp, targets = targets, epred = epred,
       re_formula = re_formula, use_re = use_re,
       route = if (design) "design" else "grid", misc = misc)
}

#' @noRd
emm_check_name <- function(x, arg) {
  if (!is.null(x) && (!is.character(x) || length(x) != 1L || is.na(x))) {
    frm_stop("`", arg, "` must be a single name, not ", arg_desc(x),
             call. = FALSE)
  }
  invisible(x)
}

#' The predictor one response contributes, with brms's refusals for a
#' name the model does not have.
#'
#' @noRd
emm_target_one <- function(object, r, dpar, nlpar, epred) {
  rspec <- object$spec$responses[[r]]
  fam <- rspec$family
  where <- if (length(object$spec$responses) > 1L) {
    paste0(" for response '", r, "'")
  } else ""
  if (epred) {
    emm_check_epred_family(rspec, where)
    # the nonlinear parameters reach the mean only through a dpar's
    # body, which frm_lp_basis() differentiates as a whole
    return(list(resp = r, rspec = rspec,
                dpars = setdiff(names(rspec$dpars), rspec$nlpars)))
  }
  own_dpars <- setdiff(names(rspec$dpars), rspec$nlpars)
  nm <- if (!is.null(nlpar)) {
    if (!nlpar %in% rspec$nlpars) {
      frm_stop("Non-linear parameter '", nlpar, "' is not part of the ",
               "model", where, ". Supported parameters are: ",
               if (length(rspec$nlpars)) {
                 paste0("'", rspec$nlpars, "'", collapse = ", ")
               } else "none, since the model has no nonlinear formula",
               call. = FALSE)
    }
    nlpar
  } else if (!is.null(dpar)) {
    # a nonlinear parameter is reachable by dpar = too, the way
    # frm_linpred() and fitted() reach it
    if (!dpar %in% names(rspec$dpars)) {
      frm_stop("Distributional parameter '", dpar, "' is not part of the ",
               "model", where, ". Supported parameters are: ",
               paste0("'", own_dpars, "'", collapse = ", "), call. = FALSE)
    }
    dpar
  } else if ("mu" %in% names(rspec$dpars)) {
    "mu"
  } else {
    frm_stop("emmeans needs dpar =", where, ": family '",
             fam[["family"]], "' has no parameter named mu. Name one ",
             "of: ", paste0("'", own_dpars, "'", collapse = ", "),
             call. = FALSE)
  }
  lp <- object$frame[["linpreds"]][[linpred_key(r, nm)]]
  if (is.null(lp)) {
    frm_stop("'", nm, "'", where, " is fixed to a constant, so it has no ",
             "predictor for emmeans to average", call. = FALSE)
  }
  link <- lp[["link"]][["name"]] %||% "identity"
  # an ordinal location is latent: its inverse link maps no threshold
  # to anything, so there is no response scale to transform to
  latent <- identical(fam[["type"]], "ordinal") && nm == "mu"
  misc <- if (latent) {
    list()
  } else {
    emm_link_misc(if (nm == "mu") fam[["family"]] else "", lp[["link"]])
  }
  list(resp = r, rspec = rspec, name = nm, lp = lp,
       link = if (latent) "latent" else link, misc = misc)
}

#' The expected response is one number per grid row, or refused.
#'
#' @noRd
emm_check_epred_family <- function(rspec, where) {
  fam <- rspec$family
  if (isTRUE(fam[["type"]] %in% c("ordinal", "categorical"))) {
    frm_stop("emmeans(epred = TRUE) is refused", where, ": family '",
             fam[["family"]], "' predicts a distribution over the ",
             "categories, not one mean per row, and no delta method for ",
             "the category probabilities is built. Average the latent ",
             "predictor instead (leave epred unset, or name dpar =), or ",
             "use frm_linpred(type = \"response\") or ",
             "conditional_effects() for the probabilities", call. = FALSE)
  }
  if (!is.null(fam_structure(fam)[["fitted_mean"]])) {
    frm_stop("emmeans(epred = TRUE) is refused", where, ": the mean of ",
             "family '", fam[["family"]], "' conditions on the whole ",
             "observed response sequence, so a reference grid row has none",
             call. = FALSE)
  }
  if (is.null(fam[["post"]]$mean_fn) && !"mu" %in% names(rspec$dpars)) {
    frm_stop("emmeans(epred = TRUE) is refused", where, ": family '",
             fam[["family"]], "' declares no mean. Name a parameter with ",
             "dpar = instead", call. = FALSE)
  }
  invisible(NULL)
}

#' emmeans's link labels for one predictor. A link R's `make.link()`
#' knows goes by name, as brms passes it; any other is handed over as
#' the link object itself, which emmeans accepts, because a name it
#' does not know fails only later, at `type = "response"`.
#'
#' @noRd
emm_link_misc <- function(family, link) {
  nm <- link[["name"]] %||% "identity"
  misc <- emmeans::.std.link.labels(list(family = family, link = nm),
                                    list())
  if (!is.null(misc$tran) &&
      !nm %in% c("logit", "probit", "cauchit", "cloglog", "log", "sqrt",
                 "1/mu^2", "inverse")) {
    misc$tran <- list(linkfun = link[["linkfun"]],
                      linkinv = link[["linkinv"]],
                      mu.eta = link[["mu_eta"]],
                      valideta = function(eta) TRUE, name = nm)
  }
  misc
}

#' Whether the design basis is exact for a predictor: linear, and
#' every term parametric. A smooth's, a gp()'s, a mo()'s and an mi()'s
#' contribution is not in the fixed-effect terms, and a design built
#' from those terms alone left it out of the means without a word.
#'
#' @noRd
emm_is_design <- function(lp) {
  is.null(lp[["nl_body"]]) && !length(lp[["smooths"]]) &&
    !length(lp[["gps"]]) && !length(lp[["mo"]]) && !length(lp[["mi"]])
}

#' The predictors a grid-route prediction reads: every covariate of the
#' selected predictor and of the parameters its body reaches, the
#' grouping factors when `re_formula` keeps them, and for the
#' expected response the addition terms its mean reads (trials,
#' truncation bounds), which brms's grid holds too.
#'
#' @noRd
emm_grid_vars <- function(object, t, tg) {
  rspec <- t$rspec
  r <- t$resp
  lps <- if (tg$epred) {
    lapply(t$dpars, function(d) {
      object$frame[["linpreds"]][[linpred_key(r, d)]]
    })
  } else {
    list(t$lp)
  }
  v <- character(0)
  for (lp in Filter(Negate(is.null), lps)) {
    v <- c(v, ce_plot_vars(object, rspec, lp, r))
    for (lpk in emm_reach(object, rspec, lp, r)) {
      for (gi in lpk[["gps"]] %||% list()) {
        v <- c(v, unlist(lapply(gi$exprs, all.vars)))
      }
      if (tg$use_re) {
        gv <- emm_group_vars(rspec, lpk[["dpar"]])
        # a factor whose terms re_formula drops would only multiply the
        # grid by levels that change nothing
        if (inherits(tg$re_formula, "formula")) {
          gv <- intersect(gv, all.vars(tg$re_formula))
        }
        v <- c(v, gv)
      }
    }
  }
  if (tg$epred) {
    skip <- c("cens", "cens_y2", "se_sigma", "mi", "mi_sd", "weights")
    for (nm in setdiff(names(rspec$aterms), skip)) {
      v <- c(v, all.vars(rspec$aterms[[nm]]))
    }
  }
  # a body also names parameters and constants, which are not columns
  intersect(unique(v), names(model.frame(object)))
}

#' A predictor and every predictor its nonlinear body reaches.
#'
#' @noRd
emm_reach <- function(object, rspec, lp, r, seen = character(0)) {
  out <- list(lp)
  if (is.null(lp[["nl_body"]])) return(out)
  reach <- c(lp[["nl_pars"]] %||% rspec$nlpars, lp[["nl_dpar_refs"]])
  for (np in setdiff(reach, seen)) {
    lpn <- object$frame[["linpreds"]][[linpred_key(r, np)]]
    if (!is.null(lpn)) {
      out <- c(out, emm_reach(object, rspec, lpn, r, c(seen, np)))
    }
  }
  out
}

#' The grouping variables of one parameter's group-level terms, read
#' the way `nonpredictor_frame_vars()` reads them.
#'
#' @noRd
emm_group_vars <- function(rspec, dpar) {
  dp <- rspec$dpars[[dpar]]
  out <- character(0)
  for (rt in dp[["re"]] %||% list()) {
    out <- c(out, if (is.null(rt$mm)) all.vars(rt$bar[[3L]]) else
                    rt$mm$gvars)
  }
  for (ce in c(dp[["carterms"]] %||% list(),
               dp[["spdeterms"]] %||% list())) {
    out <- c(out, all.vars(ce$gr_expr))
  }
  out
}

#' The terms emmeans builds its reference grid from.
#'
#' @noRd
emm_terms <- function(object, tg) {
  if (tg$route == "design" && length(tg$targets) == 1L) {
    return(stats::delete.response(tg$targets[[1L]]$lp[["terms"]]))
  }
  vars <- if (tg$route == "design") {
    unlist(lapply(tg$targets, function(t) {
      all.vars(stats::delete.response(t$lp[["terms"]]))
    }))
  } else {
    unlist(lapply(tg$targets, emm_grid_vars, object = object, tg = tg))
  }
  vars <- unique(vars)
  rhs <- if (length(vars)) {
    Reduce(function(a, b) call("+", a, b), lapply(vars, as.name))
  } else 1
  f <- stats::as.formula(call("~", rhs))
  environment(f) <- tg$targets[[1L]]$rspec$formula_env %||% globalenv()
  stats::terms(f)
}

#' The call emmeans reads a response transformation from: it takes the
#' formula in the call's first argument and reads any function on its
#' left-hand side as a transformation. A transformation applies to the
#' location and to the expected response of one response, so only
#' those get a formula; brms passes a bare call("brms") everywhere
#' else, and so does this.
#'
#' The formula is `<response> ~ 1`, not the fit's own call. The fit's
#' formula carries brms's addition terms on the left, and emmeans read
#' `y | weights(w)` as a transformation named "|.weights", which
#' mislabeled the scale and made `type = "response"` a no-op with a
#' note.
#'
#' @noRd
emm_call <- function(object, tg) {
  own <- length(tg$targets) == 1L &&
    (tg$epred || identical(tg$targets[[1L]]$name, "mu"))
  if (!own) return(call("frm"))
  rspec <- tg$targets[[1L]]$rspec
  f <- stats::as.formula(call("~", rspec$resp_expr, 1),
                         env = rspec$formula_env %||% globalenv())
  as.call(list(as.name("frm"), formula = f))
}

#' @exportS3Method emmeans::recover_data
recover_data.frmtmb_fit <- function(object, ..., data = NULL, resp = NULL,
                                    dpar = NULL, nlpar = NULL,
                                    re_formula = NA, epred = FALSE) {
  # emmeans calls this inside try() and replaces whatever it raises with
  # "Perhaps a 'data' or 'params' argument is needed", which names
  # neither the reason nor the fix. A character value is emmeans's own
  # way for recover_data() to refuse: ref_grid() stops with that text.
  # So every refusal is returned, not raised.
  tryCatch({
    tg <- emm_target(object, resp = resp, dpar = dpar, nlpar = nlpar,
                     re_formula = re_formula, epred = epred)
    emmeans::recover_data(emm_call(object, tg), emm_terms(object, tg),
                          na.action = NULL,
                          data = data %||% model.frame(object), ...)
  }, error = function(e) conditionMessage(e))
}

#' @exportS3Method emmeans::emm_basis
emm_basis.frmtmb_fit <- function(object, trms, xlev, grid, ..., resp = NULL,
                                 dpar = NULL, nlpar = NULL, re_formula = NA,
                                 epred = FALSE) {
  tg <- emm_target(object, resp = resp, dpar = dpar, nlpar = nlpar,
                   re_formula = re_formula, epred = epred, warn = TRUE)
  b <- if (tg$route == "design") {
    emm_basis_design(object, tg, xlev, grid)
  } else {
    emm_basis_grid(object, tg, grid)
  }
  list(X = b$X, bhat = b$bhat, nbasis = matrix(NA), V = b$V,
       dffun = function(k, dfargs) Inf, dfargs = list(), misc = tg$misc)
}

#' The design basis: each selected predictor's fixed-effect design at the
#' grid, block diagonal over the responses so the rows are the grid once
#' per response, response slowest, which is where emmeans puts
#' `rep.meas`.
#'
#' @noRd
emm_basis_design <- function(object, tg, xlev, grid) {
  Xs <- list()
  bh <- list()
  pos <- list()
  for (t in tg$targets) {
    lp <- t$lp
    trm <- stats::delete.response(lp[["terms"]])
    xl <- xlev[intersect(names(xlev), all.vars(trm))]
    m <- stats::model.frame(trm, grid, na.action = stats::na.pass,
                            xlev = xl)
    X <- stats::model.matrix(trm, m, contrasts.arg = lp[["contrasts"]])
    # The fitted design is not always model.matrix()'s: an ordinal family
    # drops the intercept (the K-1 thresholds take its place), and a
    # rank-deficient fit drops aliased columns. Select the fitted columns
    # by name, or the basis is not conformable with bhat and emmeans
    # fails with "Non-conformable elements in reference grid".
    pn <- lp[["param_colnames"]][seq_len(lp[["n_param_cols"]])]
    if (!identical(colnames(X), pn)) {
      keep <- match(pn, colnames(X))
      if (anyNA(keep)) {
        frm_stop("emmeans support cannot rebuild the fitted design: ",
                 "column(s) ", paste(pn[is.na(keep)], collapse = ", "),
                 " are missing from the reference grid", call. = FALSE)
      }
      X <- X[, keep, drop = FALSE]
    }
    idx <- lp[["idx"]][seq_len(lp[["n_param_cols"]])]
    Xs[[t$resp]] <- X
    bh[[t$resp]] <- object$estimates[[lp[["par"]]]][idx]
    pos[[t$resp]] <- emm_vcov_pos(object, lp[["par"]], idx)
  }
  nr <- vapply(Xs, nrow, 1L)
  nc <- vapply(Xs, ncol, 1L)
  X <- matrix(0, sum(nr), sum(nc))
  for (k in seq_along(Xs)) {
    X[sum(nr[seq_len(k - 1L)]) + seq_len(nr[k]),
      sum(nc[seq_len(k - 1L)]) + seq_len(nc[k])] <- Xs[[k]]
  }
  colnames(X) <- unlist(lapply(Xs, colnames), use.names = FALSE)
  p <- unlist(pos, use.names = FALSE)
  # vcov_estimated(), not vcov(): these positions index the estimated
  # coefficient vector, and vcov() is brms's population-level block. A
  # coefficient fixed to a constant has no row there and no variance.
  Vall <- vcov_estimated(object)
  V <- matrix(0, length(p), length(p))
  ok <- !is.na(p)
  V[ok, ok] <- Vall[p[ok], p[ok], drop = FALSE]
  list(X = X, bhat = unname(unlist(bh, use.names = FALSE)), V = V)
}

#' Rows of `vcov_estimated()` for a predictor's coefficients: beta
#' first, then the ESTIMATED betad entries in their own order. `NA`
#' marks a coefficient fixed to a constant.
#'
#' @noRd
emm_vcov_pos <- function(object, par, idx) {
  tpl <- object$frame[["par_template"]]
  if (par == "beta") return(idx)
  est <- setdiff(seq_along(tpl[["betad"]]), object$frame[["betad_fixed_idx"]])
  length(tpl[["beta"]]) + match(idx, est)
}

#' The grid basis: `X` the identity over the grid (once per response),
#' `bhat` the predictions, `V = G Sigma G'` over the joint covariance.
#' Responses share coefficients only through a `|ID|` random-effect
#' block, and the joint covariance carries exactly that.
#'
#' @noRd
emm_basis_grid <- function(object, tg, grid) {
  nd <- as.data.frame(grid)
  parts <- lapply(tg$targets, function(t) {
    p <- if (tg$epred) {
      emm_epred_part(object, t, nd, tg$re_formula)
    } else {
      emm_lp_part(object, t$resp, t$name, nd, tg$re_formula)
    }
    emm_check_part(p, t$resp, length(tg$targets) > 1L)
    p
  })
  jc <- get_joint_cov(object)
  pos_all <- sort(unique(unlist(lapply(parts, `[[`, "pos"))))
  n <- vapply(parts, function(p) length(p$est), 1L)
  G <- matrix(0, sum(n), length(pos_all))
  off <- 0L
  for (p in parts) {
    rows <- off + seq_along(p$est)
    cols <- match(p$pos, pos_all)
    for (j in seq_along(cols)) {
      G[rows, cols[j]] <- G[rows, cols[j]] + p$G[, j]
    }
    off <- off + length(p$est)
  }
  V <- G %*% jc$V[pos_all, pos_all, drop = FALSE] %*% t(G)
  list(X = diag(sum(n)), bhat = unname(unlist(lapply(parts, `[[`, "est"))),
       V = (V + t(V)) / 2)
}

#' One predictor at the grid, through `frm_lp_basis()`, which tapes a
#' nonlinear body for its Jacobian.
#'
#' @noRd
emm_lp_part <- function(object, r, name, nd, re_formula) {
  b <- frm_lp_basis(object, newdata = nd, dpar = name, resp = r,
                    re_formula = re_formula)
  list(est = unname(b$eta), G = as.matrix(b$A), pos = b$coef_pos,
       extra = b$extra_var, nonest = b$nonest)
}

#' The expected response at the grid and its Jacobian, by the chain rule
#' `predict_mean_se()` uses: the mean's gradient in each dpar's linear
#' predictor (`mean_eta_grad()`) times that predictor's Jacobian from
#' `frm_lp_basis()`. Unlike `predict_mean_se()` it keeps the whole
#' Jacobian rather than the diagonal of `G V G'`, and it reaches a
#' nonlinear dpar, whose Jacobian `frm_lp_basis()` tapes.
#'
#' @noRd
emm_epred_part <- function(object, t, nd, re_formula) {
  rspec <- t$rspec
  fam <- rspec$family
  bs <- list()
  dp <- list()
  links <- list()
  for (dnm in t$dpars) {
    lp <- object$frame[["linpreds"]][[linpred_key(t$resp, dnm)]]
    if (is.null(lp)) next
    bs[[dnm]] <- frm_lp_basis(object, newdata = nd, dpar = dnm,
                              resp = t$resp, re_formula = re_formula)
    links[[dnm]] <- lp[["link"]]
    dp[[dnm]] <- lp[["link"]]$linkinv(bs[[dnm]]$eta)
  }
  av <- aterms_for_newdata(rspec, nd)
  m <- response_mean(fam, dp, av)
  n <- nrow(nd)
  if (!is.numeric(m) || !is.null(dim(m)) || length(m) != n) {
    frm_stop("emmeans(epred = TRUE) is refused: the mean of family '",
             fam[["family"]], "' is not one number per grid row",
             call. = FALSE)
  }
  pos <- unique(unlist(lapply(bs, `[[`, "coef_pos")))
  G <- matrix(0, n, length(pos))
  for (dnm in names(bs)) {
    g <- mean_eta_grad(fam, dp, av, dnm, links[[dnm]], bs[[dnm]]$eta)
    A <- as.matrix(bs[[dnm]]$A)
    cols <- match(bs[[dnm]]$coef_pos, pos)
    for (j in seq_along(cols)) {
      G[, cols[j]] <- G[, cols[j]] + g * A[, j]
    }
  }
  extra <- Reduce(`+`, lapply(bs, `[[`, "extra_var"), numeric(n))
  nonest <- Reduce(`|`, lapply(bs, `[[`, "nonest"), rep(FALSE, n))
  list(est = unname(as.numeric(m)), G = G, pos = pos, extra = extra,
       nonest = nonest)
}

#' The grid basis carries coefficient uncertainty only. A row that also
#' carries a kriging variance, or that the fit cannot estimate, is
#' refused rather than reported with a covariance that is too small or
#' an estimate that is not one.
#'
#' @noRd
emm_check_part <- function(p, r, several) {
  where <- if (several) paste0(" for response '", r, "'") else ""
  if (any(p$nonest) || !all(is.finite(p$est))) {
    frm_stop("emmeans cannot estimate every row of the reference grid",
             where, ": the fit is rank deficient there, or the prediction ",
             "is not finite. Restrict the grid with at = or by =",
             call. = FALSE)
  }
  if (any(p$extra > 0)) {
    frm_stop("emmeans cannot use this reference grid", where, ": an exact ",
             "gp() term is predicted at a position the fit did not see. ",
             "Its kriging variance is known point by point but not as a ",
             "covariance between two grid points, which the marginal means ",
             "need. Use an approximate gp(..., k = ), or hold the gp() ",
             "covariate at an observed value with at =", call. = FALSE)
  }
  invisible(NULL)
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
  Filter(function(bk) !bk[["covstruct"]] %in% c("smooth", "gp", "hsgp") &&
           !is.null(bk[["components"]][[1L]]$bar),
         object$frame[["re_blocks"]])
}

#' Labels for the random-effect coefficient vector, which is also the
#' column order of every Z: level-major within a block, one entry per
#' (level, term coefficient), the way lme4 labels Zt rows.
#'
#' @noRd
re_coef_labels <- function(frame) {
  if (!length(frame[["re_blocks"]])) return(character(0))
  unlist(lapply(frame[["re_blocks"]], function(bk) {
    lv <- bk[["levels"]] %||% as.character(seq_len(bk[["n_levels"]]))
    paste0(rep(lv,
               each = bk[["dim"]]), ".", rep(bk[["cnms"]], bk[["n_levels"]]))
  }), use.names = FALSE)
}

#' The grouping factors as they were at fit time, rebuilt from the
#' stored model frame the same way predict() rebuilds them for newdata.
#'
#' @noRd
getME_flist <- function(object) {
  out <- list()
  for (bk in getME_group_blocks(object)) {
    comp <- bk[["components"]][[1L]]
    # a multi-membership row belongs to several levels at once, so there
    # is no per-observation grouping factor to report; lme4's flist has
    # no representation for one, and returning the first member would
    # be a wrong answer rather than a missing one
    if (!is.null(comp$mm)) next
    lp <- object$frame[["linpreds"]][[comp$lp_key]]
    env <- object$spec$responses[[lp[["resp"]]]]$formula_env
    gv <- tryCatch(
      as.character(eval(comp$bar[[3L]], object$frame[["data_frame"]], env)),
      error = function(e) NULL
    )
    if (length(gv) != object$frame[["n_obs"]]) {
      frm_stop("getME(\"flist\"): cannot rebuild the grouping factor for `",
               bk[["term_label"]], "`", call. = FALSE)
    }
    nm <- bk[["group_name"]]
    if (is.null(out[[nm]])) out[[nm]] <- factor(gv, levels = bk[["levels"]])
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
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
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
  frm_check_dots(...)
  if (missing(name) || !is.character(name) || !length(name)) {
    frm_stop("getME() needs one or more names: ",
             paste(frmtmb_getME_vocab, collapse = ", "), call. = FALSE)
  }
  bad <- setdiff(name, frmtmb_getME_vocab)
  if (length(bad)) {
    frm_stop("getME(): unknown name(s) ", paste(bad, collapse = ", "),
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
                                dims = c(object$frame[["n_obs"]],
                                         object$frame[["n_c"]] %||% 0L))
    }
    dimnames(Z) <- list(rownames(object$frame[["data_frame"]]),
                        re_coef_labels(object$frame))
    Z
  }
  theta <- object$estimates[["theta"]] %||% numeric(0)
  names(theta) <- if (length(theta)) {
    paste0("theta_", seq_along(theta))
  }
  switch(
    name,
    X = mu_lp()$X,
    Z = mu_Z(),
    # forced before t(): an error while S4 dispatch evaluates the
    # argument is raised again as a plain simpleError, losing its class
    Zt = { z <- mu_Z(); Matrix::t(z) },
    beta = ,
    fixef = object$estimates[["beta"]],
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
