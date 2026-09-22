# Method surface for frm_sample() output: every draw is a complete
# parameter vector (tmbstan samples the random effects too), so the
# fitted-model machinery runs per draw exactly.

#' Column positions of each template component inside the draws matrix
#' (which is in template order, mapped betad entries absent, `lp__` last).
#'
#' @noRd
draws_par_index <- function(fit) {
  tpl <- fit$frame[["par_template"]]
  idx <- list()
  pos <- 0L
  for (cp in names(tpl)) {
    len <- length(tpl[[cp]])
    if (cp == "betad" && length(fx <- fit$frame[["betad_fixed_idx"]])) {
      len <- len - length(fx)
    }
    idx[[cp]] <- pos + seq_len(len)
    pos <- pos + len
  }
  idx
}

#' The originating fit stripped of the "no maximum-likelihood estimate"
#' marker, for the draws methods that use it only as a structural
#' template (its frame, its spec, its block layout). Nothing here reads
#' its `estimates` as an estimate.
#'
#' @noRd
draws_base_fit <- function(x) {
  fit <- x$fit
  class(fit) <- setdiff(class(fit), "frmtmb_unfitted")
  fit
}

#' brms's `point_estimate`: the draws collapsed to their mean or median,
#' one parameter vector, repeated `ndraws_point_estimate` times.
#'
#' It is a PARAMETER-space operation, as in brms, not a summary of the
#' output: every predictive method then runs once at that one vector.
#' It was accepted and ignored before, so `posterior_epred(ds,
#' point_estimate = "median", ndraws_point_estimate = 2)` returned all
#' 25 draws where brms returns 2 (defect S1 of
#' dev/brmsport-findings.md).
#'
#' @noRd
draws_at_point_estimate <- function(x, point_estimate,
                                    ndraws_point_estimate = 1) {
  if (is.null(point_estimate)) return(x)
  pe <- frm_match_arg(point_estimate, c("mean", "median"))
  n <- ndraws_point_estimate
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1 ||
      n != round(n)) {
    frm_stop("`ndraws_point_estimate` must be one whole number of at ",
             "least 1", call. = FALSE)
  }
  f <- if (identical(pe, "mean")) mean else stats::median
  v <- apply(x$draws, 2L, f)
  x$draws <- matrix(rep(v, each = as.integer(n)), as.integer(n),
                    ncol(x$draws),
                    dimnames = list(NULL, colnames(x$draws)))
  x
}

#' The originating fit with its estimates replaced by one draw. The draw
#' IS a parameter vector, so the object is a legitimate fit from here on
#' even when the draws came from a formula with no ML mode behind it.
#'
#' @noRd
draws_fit_at <- function(x, i, idx = draws_par_index(x$fit)) {
  fit <- draws_base_fit(x)
  est <- fit$frame[["par_template"]]   # mapped betad entries keep link(const)
  row <- draws_internal_matrix(x, i)[1L, ]
  for (cp in names(idx)) {
    if (cp == "betad" && length(fx <- fit$frame[["betad_fixed_idx"]])) {
      pos <- setdiff(seq_along(est[[cp]]), fx)
      est[[cp]][pos] <- row[idx[[cp]]]
    } else {
      est[[cp]][] <- row[idx[[cp]]]
    }
  }
  for (cp in names(fit$estimates)) {
    names(est[[cp]]) <- names(fit$estimates[[cp]])
  }
  fit$estimates <- est
  fit$cache <- new.env(parent = emptyenv())   # no stale sdreport
  fit
}

#' Row indices of an evenly spaced subsample of the draws, or of the
#' draws named by `draw_ids`. Predictive methods run the whole model per
#' draw, so a thinned set keeps them affordable without favoring one
#' part of the chain.
#'
#' @noRd
draws_subsample <- function(x, ndraws, draw_ids = NULL) {
  n <- nrow(x$draws)
  if (!is.null(draw_ids)) {
    if (!is.null(ndraws)) {
      frm_stop("`ndraws` and `draw_ids` both choose which draws to use, so ",
               "only one of them can be given: `ndraws` takes an evenly ",
               "spaced subsample of that size, `draw_ids` takes the rows ",
               "you name", call. = FALSE)
    }
    ok <- is.numeric(draw_ids) && length(draw_ids) &&
      !anyNA(draw_ids) && all(draw_ids == round(draw_ids)) &&
      all(draw_ids >= 1L) && all(draw_ids <= n)
    if (!ok) {
      frm_stop("`draw_ids` must be whole numbers between 1 and ", n,
               ", the number of draws", call. = FALSE)
    }
    return(as.integer(draw_ids))
  }
  if (is.null(ndraws) || ndraws >= n) return(seq_len(n))
  round(seq(1, n, length.out = ndraws))
}

#' @export
summary.frmtmb_draws <- function(object, priors = FALSE, prob = 0.95,
                                 robust = FALSE, mc_se = FALSE, ...) {
  frm_check_dots(...)
  check_flag(priors, "priors")
  check_flag(robust, "robust")
  check_flag(mc_se, "mc_se")
  check_probability(prob, "prob")
  if (priors) {
    frm_stop("summary(priors = TRUE) has no table to add the priors to: ",
             "this summary is a matrix of the parameters, not brms's ",
             "summary object. prior_summary(object) reports the priors ",
             "the draws were taken under", call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    frm_stop("summary() on draws needs the 'posterior' package: its ",
             "columns are brms's, which posterior computes", call. = FALSE)
  }
  keep <- draws_outer_cols(object)
  probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
  # brms:::summary.brmsfit's own `.summary()`, measure for measure and
  # in its column order, so a column carries brms's definition under
  # brms's header: Estimate, MCSE, Est.Error, the two interval ends,
  # Rhat, Bulk_ESS, Tail_ESS
  qfun <- function(x, ...) {
    qs <- posterior::quantile2(x, probs = probs, ...)
    names(qs) <- paste0(c("l-", "u-"), (probs[2] - probs[1]) * 100,
                        "% CI")
    qs
  }
  measures <- list()
  if (robust) {
    measures$Estimate <- stats::median
    if (mc_se) measures$MCSE <- posterior::mcse_median
    measures$Est.Error <- stats::mad
  } else {
    measures$Estimate <- mean
    if (mc_se) measures$MCSE <- posterior::mcse_mean
    measures$Est.Error <- stats::sd
  }
  measures$quantiles <- qfun
  measures$Rhat <- posterior::rhat
  measures$Bulk_ESS <- posterior::ess_bulk
  measures$Tail_ESS <- posterior::ess_tail
  a <- posterior::subset_draws(draws_as_array(object), variable = keep)
  out <- do.call(posterior::summarize_draws, c(list(a), measures))
  out <- as.data.frame(out)
  tab <- as.matrix(out[, -1L, drop = FALSE])
  rownames(tab) <- out$variable
  tab
}

#' @rdname draws-structure
#' @exportS3Method nlme::fixef
#' @export
fixef.frmtmb_draws <- function(object, summary = TRUE, robust = FALSE,
                               probs = c(0.025, 0.975), pars = NULL, ...) {
  # brms:::fixef.brmsfit, on this package's brms-named draws
  frm_check_dots(...)
  check_flag(summary, "summary")
  all_pars <- variables(object)
  fpars <- all_pars[grepl(draws_fixef_regex, all_pars)]
  if (!is.null(pars)) {
    pars <- as.character(pars)
    fpars <- fpars[sub("^[^_]+_", "", fpars) %in% pars]
  }
  thr <- draws_fixef_ordinal(object, all_pars)
  if (!is.null(thr)) {
    out <- thr
    if (!is.null(pars)) out <- out[, colnames(out) %in% pars, drop = FALSE]
    if (!ncol(out)) return(NULL)
  } else {
    if (!length(fpars)) return(NULL)
    out <- as.matrix(object, variable = fpars)
    colnames(out) <- gsub(draws_fixef_regex, "", fpars)
  }
  if (summary) out <- posterior_summary(out, probs, robust)
  out
}

#' An ordinal fit's population-level draws in brms's rows and order:
#' `Intercept[1]`, `Intercept[2]`, the coefficients, then `cs()` terms.
#'
#' The sampler stores the thresholds on the internal scale (`tau_raw`,
#' the first threshold and log increments for cumulative and sratio), so
#' the `b_` regex found only the slopes and `fixef(ds)` reported `x`
#' where the fit's `fixef()` and brms report three rows. Each draw is
#' mapped through the map the family declares, the same map the fit's
#' rows use. `NULL` when the fit has nothing outside the coefficients,
#' or when a column the rows need is not in the draws.
#'
#' @noRd
draws_fixef_ordinal <- function(object, all_pars) {
  fit <- draws_base_fit(object)
  rows <- brms_fixef_rows(fit)
  if (!length(rows$extra)) return(NULL)
  tab <- brms_coef_table(fit)
  need <- unique(c(tab$brms[rows$idx[!is.na(rows$idx)]],
                   unlist(lapply(rows$extra, function(e) {
                     paste0(e$comp, "_", seq_along(e$raw))
                   }))))
  if (!all(need %in% all_pars)) return(NULL)
  m <- as.matrix(object, variable = need)
  out <- matrix(NA_real_, nrow(m), length(rows$names),
                dimnames = list(NULL, rows$names))
  co <- !is.na(rows$idx)
  out[, co] <- m[, tab$brms[rows$idx[co]], drop = FALSE]
  for (b in seq_along(rows$extra)) {
    e <- rows$extra[[b]]
    raw <- m[, paste0(e$comp, "_", seq_along(e$raw)), drop = FALSE]
    vals <- matrix(apply(raw, 1L, e$map), ncol = length(e$raw),
                   byrow = TRUE)
    at <- which(rows$blk == b)
    out[, at] <- vals[, rows$pos[at], drop = FALSE]
  }
  out
}

#' @rdname draws-structure
#' @exportS3Method nlme::VarCorr
#' @export
VarCorr.frmtmb_draws <- function(x, sigma = 1, summary = TRUE,
                                 robust = FALSE,
                                 probs = c(0.025, 0.975), ...) {
  # brms:::VarCorr.brmsfit's layout and summary path. brms reads sd_ and
  # cor_ draws the sampler stored; these draws store the covariance
  # parameters on their unconstrained scale, so the standard deviations
  # and correlations are computed per draw from them first
  frm_check_dots(...)
  check_flag(summary, "summary")
  fit <- draws_base_fit(x)
  lay <- varcorr_layout(fit)
  if (!length(lay$groups) && is.null(lay$residual)) {
    frm_stop("The model does not contain covariance matrices.", call. = FALSE)
  }
  if (length(lay$groups) && draws_is_laplace(x) &&
        is.null(draws_par_index(fit)$theta)) {
    frm_stop("VarCorr() found no covariance parameters in these draws",
             call. = FALSE)
  }
  per <- draws_varcorr_values(x, lay)
  keys <- names(per[[1L]])
  out <- list()
  for (key in keys) {
    rn <- if (identical(key, "residual__")) lay$residual$rnames else
      lay$groups[[key]]$rnames
    K <- length(rn)
    sdm <- draws_derived_matrix(x, matrix(vapply(per, function(v) {
      v[[key]]$sd
    }, numeric(K)), ncol = K, byrow = TRUE), rn)
    e <- list(sd = sdm)
    if (!is.null(per[[1L]][[key]]$cor)) {
      # brms's get_cor_matrix() and get_cov_matrix(), from the sd draws
      # and the lower-triangle correlation draws, so the arithmetic is
      # brms's and not a covariance matrix read back
      lt <- draws_cor_index(K)
      cor_draws <- matrix(vapply(per, function(v) v[[key]]$cor[lt],
                                 numeric(length(lt))),
                          ncol = length(lt), byrow = TRUE)
      e$cor <- draws_cor_array(cor_draws, K)
      e$cov <- draws_cov_array(unname(e$sd), cor_draws)
      dimnames(e$cor)[2:3] <- list(rn, rn)
      dimnames(e$cov)[2:3] <- list(rn, rn)
    }
    if (summary) {
      e$sd <- posterior_summary(e$sd, probs, robust)
      if (!is.null(e$cor)) {
        e$cor <- posterior_summary(e$cor, probs, robust)
        e$cov <- posterior_summary(e$cov, probs, robust)
      }
    }
    # brms's element order: sd, then cor and cov
    out[[key]] <- e[intersect(c("sd", "cor", "cov"), names(e))]
  }
  out
}

#' @exportS3Method rstantools::prior_summary
#' @export
prior_summary.frmtmb_draws <- function(object, ...) {
  frm_check_dots(...)
  pl <- object$fit$prior
  if (is.null(pl) || (!length(unclass(pl)) &&
                        !length(attr(pl, "overrides")))) {
    cat("No priors were used (flat improper priors on the outer ",
        "parameters).\n", sep = "")
    return(invisible(NULL))
  }
  pl
}

#' @rdname draws-structure
#' @exportS3Method nlme::ranef
#' @export
ranef.frmtmb_draws <- function(object, summary = TRUE, robust = FALSE,
                               probs = c(0.025, 0.975), pars = NULL,
                               groups = NULL, ...) {
  # brms:::ranef.brmsfit's shape: per grouping factor, draws x levels x
  # coefficients, summarized to levels x statistics x coefficients
  frm_check_dots(...)
  check_flag(summary, "summary")
  fit <- draws_base_fit(object)
  lay <- draws_ranef_layout(object)
  if (!length(lay)) {
    frm_stop("The model does not contain group-level effects.", call. = FALSE)
  }
  if (!is.null(pars)) pars <- as.character(pars)
  keep <- names(lay)
  if (!is.null(groups)) keep <- intersect(keep, as.character(groups))
  m <- object$draws
  out <- list()
  for (g in keep) {
    L <- lay[[g]]
    sel <- seq_along(L$coefs)
    if (!is.null(pars)) sel <- sel[L$coef[sel] %in% pars]
    if (!length(sel)) next
    A <- array(NA_real_, c(nrow(m), length(L$levels), length(sel)))
    for (j in seq_along(sel)) {
      cols <- L$cols[, sel[j]]
      ok <- !is.na(cols)
      if (any(ok)) A[, ok, j] <- m[, cols[ok], drop = FALSE]
    }
    A <- draws_ranef_fill(object, L, sel, A, g)
    dimnames(A) <- list(NULL, L$levels, L$coefs[sel])
    # brms's draws array keeps the chain count its as.matrix() carried
    attr(A, "nchains") <- nchains(object)
    if (summary) A <- posterior_summary(A, probs, robust)
    out[[g]] <- A
  }
  out
}

#' @exportS3Method brms::hypothesis
#' @export
hypothesis.frmtmb_draws <- function(x, hypothesis, class = "b", group = "",
                                    scope = c("standard", "ranef", "coef"),
                                    alpha = 0.05, robust = FALSE,
                                    seed = NULL, ...) {
  frm_check_dots(...)
  scope <- frm_match_arg(scope)
  check_flag(robust, "robust")
  check_probability(alpha, "alpha")
  if (!is.null(seed)) set.seed(seed)
  if (!is.character(hypothesis) || !length(hypothesis)) {
    frm_stop("Argument 'hypothesis' must be a character vector.", call. = FALSE)
  }
  if (scope != "standard") {
    return(draws_hypothesis_coef(x, hypothesis, group, scope, alpha,
                                 robust))
  }
  fit <- x$fit
  vo <- hyp_vals_only(fit)
  prefix <- hyp_class_prefix(class, group)
  env_names <- names(hyp_env_vals(fit, vo$vals, vo$comp))
  # brms evaluates over variables(x), which on draws includes the stored
  # columns (r_g[1,Intercept], lp__) as well as the derived sd_ and cor_
  hp <- hyp_parse_all(hypothesis, union(env_names, colnames(x$draws)),
                      class, group)
  exs <- hp$exprs
  labels <- hyp_labels(hypothesis)
  used <- unique(unlist(lapply(exs, hyp_expr_vars)))
  from_cols <- setdiff(used, env_names)
  # `x$draws[i, name]` reads the FIRST column of a name given twice, so
  # a name that matches two columns is refused rather than half-read.
  # frm_sample() suffixes repeats as brms does, so this guards draws
  # built any other way
  cn <- colnames(x$draws)
  twice <- intersect(from_cols, cn[duplicated(cn)])
  if (length(twice)) {
    frm_stop("hypothesis() cannot read ", paste0("'", twice, "'",
                                                 collapse = ", "),
             ": the draws carry that name on more than one column",
             call. = FALSE)
  }
  need_env <- length(intersect(used, env_names)) > 0L
  idx <- draws_par_index(fit)
  n <- nrow(x$draws)
  draws <- matrix(NA_real_, n, length(exs),
                  dimnames = list(NULL, labels))
  for (i in seq_len(n)) {
    vals <- as.list(x$draws[i, from_cols])
    names(vals) <- from_cols
    if (need_env) {
      sh <- draws_fit_at(x, i, idx)
      w <- hyp_vals_only(sh)
      vals <- c(hyp_env_vals(sh, w$vals, w$comp), vals)
    }
    draws[i, ] <- vapply(exs, function(ex) hyp_eval_in(ex, vals),
                         numeric(1))
  }
  # A POINT null is the spelling that carries "=". The bare-quantity
  # spelling this package also accepts parses `two.sided` as well, and
  # it is a summary request rather than a test, so it is not weighed.
  is_point <- hp$dir == "two.sided" &
    grepl("=", hypothesis, fixed = TRUE)
  ev <- er_evidence(x, exs, hp$dir, draws, hypothesis, is_point)
  k_n <- length(exs)
  ctr <- if (robust) stats::median else mean
  spr <- if (robust) stats::mad else stats::sd
  lo_p <- ifelse(hp$dir == "two.sided", alpha / 2, alpha)
  est <- apply(draws, 2L, ctr)
  err <- apply(draws, 2L, spr)
  lwr <- vapply(seq_len(k_n), function(k) {
    unname(stats::quantile(draws[, k], lo_p[k]))
  }, 1)
  upr <- vapply(seq_len(k_n), function(k) {
    unname(stats::quantile(draws[, k], 1 - lo_p[k]))
  }, 1)
  er <- ev$evid_ratio
  pp <- ifelse(is.infinite(er), 1, er / (1 + er))
  hyp_brms_result(
    labels, est, err, lwr, upr, er, pp, hp$dir, pp > 1 - alpha,
    hyp_samples_frame(draws, k_n),
    hyp_samples_frame(matrix(NA, n, k_n), k_n), prefix, alpha,
    list(method = "posterior", evid_ratio_mcse = ev$mcse))
}

#' Expected-value and predictive draws from sampled parameters
#'
#' `posterior_epred()` evaluates the response-scale expectation per
#' draw; `posterior_predict()` additionally simulates responses from
#' the family, giving the posterior predictive distribution. Both
#' condition on each draw's own random effects (`re_formula = NA` drops
#' them; `re.form` is an accepted alias here, see *Argument
#' spellings*).
#'
#' @section Categorical outcomes:
#' An ordinal family predicts a DISTRIBUTION per observation, not one
#' number: each draw's `frm_linpred(type = "response")` is an `n x K`
#' matrix of category probabilities. Those stack into a 3-D
#' `draws x observations x categories` array. `dimnames` are
#' `list(NULL, <observation names or NULL>, <category levels>)`, so
#' `ep[, , "high"]` is the draws-by-observations matrix for one
#' category and `ep[k, , ]` is draw `k`'s own `n x K` prediction, the
#' matrix `frm_linpred(type = "response")` returns. Every `ep[k, i, ]`
#' sums to 1 for an ordinal family.
#'
#' This is brms's convention: `?brms::posterior_epred.brmsfit`
#' documents "an S x N x C array" for categorical and ordinal models
#' and an S x N matrix otherwise, and frmtmb follows brms spelling for
#' brms-origin functions. Any family whose per-draw response-scale
#' prediction is a matrix takes the array shape; every family that
#' predicts one number per observation keeps the plain
#' `draws x observations` matrix.
#'
#' `posterior_predict()` is unaffected for an ordinal or categorical
#' family (it draws one category per observation), and so is
#' `posterior_linpred()`, which is a statement about one distributional
#' parameter and stays an `n`-column matrix of the latent predictor.
#' What does take the array shape in `posterior_predict()` is a
#' matrix-valued RESPONSE: [frmtmb::multinomial()] counts, [frmtmb::mixture_mvn()]
#' draws and `frmtmb.latent::lca()` item codes give one row per observation, so the
#' draws stack into `draws x observations x columns`.
#'
#' @section Structured draws:
#' `posterior_predict()` uses the same simulator [simulate()] does,
#' including the structured families (`frmtmb.latent::hmm()`, `mixture(groups = )`,
#' [frmtmb::mixture_mvn()]) and residual correlation terms; see the Structured
#' draws section of [frmtmb::simulate.frmtmb_fit()]. Those draws index the rows
#' the model was fitted on, so `newdata` is refused for them.
#'
#' @section Argument spellings:
#' One rule decides every name in this package: where lme4 or glmmTMB
#' and brms disagree, brms wins. The random-effect switch is therefore
#' `re_formula` everywhere, on the draws methods here and on
#' [frmtmb::predict.frmtmb_fit()] and
#' [frmtmb::simulate.frmtmb_fit()] alike; lme4's `re.form` was dropped
#' from the fit surface and is refused there by name.
#'
#' Five methods take BOTH spellings, and they are exactly the five
#' where brms ITSELF accepts both. Four declare them:
#' `posterior_epred()`, `posterior_linpred()`, `posterior_predict()`
#' and `predictive_error()` carry `re_formula` and `re.form` side by
#' side on `brmsfit`. The fifth does not declare them and accepts them
#' anyway: `predictive_interval.brmsfit()`'s whole body is
#' `posterior_predict(object, ...)`, so the alias reaches a formal one
#' frame down. What brms ACCEPTS is the test, not what it declares.
#'
#' `pp_check()` is the one method that lost the alias, and the same
#' test is why it stays lost. brms DOES honor `re.form` there, through
#' the same dots forwarding, but it warns "unrecognized and ignored"
#' while doing it. Matching brms means matching what brms decided, and
#' a warn-then-honor path is a leak rather than a decision: the
#' argument changes the answer and the message says it did not.
#' `predictive_interval()` is the contrast, where brms honors the alias
#' silently and this package follows.
#'
#' Giving both at once is refused rather than resolved. Two names for
#' one setting supplied together is a question about what was meant, and
#' guessing at it would silently ignore one of them.
#'
#' The argument ORDER is brms's too, so a positional brms call means
#' the same thing here: `newdata` then `re_formula` then `re.form` then
#' `resp`, after `transform` in `posterior_linpred()` and before it in
#' `posterior_predict()`.
#'
#' The literal default of both formals is an internal "not supplied"
#' marker rather than a value, because `NULL` (keep the random effects)
#' and `NA` (drop them) are both real settings here and neither can
#' double as "unset". The behavior when neither is given is unchanged:
#' `NULL` on every draws method, `NA` on `pp_check()` for a fit.
#'
#' @param object A `frmtmb_draws` from [frm_sample()].
#' @param newdata,resp As in [frmtmb::predict.frmtmb_fit()].
#' @param re_formula The random-effect switch, in brms's spelling:
#'   `NULL` (the default) conditions on each draw's own random effects,
#'   `NA` or `~0` gives the population-level quantity. Its meaning is
#'   [frmtmb::predict.frmtmb_fit()]'s `re_formula`; see *Argument
#'   spellings*.
#' @param re.form lme4's spelling of `re_formula`, accepted on the four
#'   methods where brms accepts it too. Pass one or the other, not
#'   both.
#' @param ndraws Number of draws to use (default: all).
#' @param draw_ids The draws to use, by row index, instead of the
#'   evenly spaced subsample `ndraws` takes. Give one or the other.
#' @param sort brms's argument. Rows come back in the order of the data
#'   here, always, so `sort = TRUE` is refused by name.
#' @param ntrys,cores brms's arguments, carried so that a positional
#'   brms call lands where brms lands it. Both are refused by name: the
#'   rejection limit of a `trunc()`ed draw is the family simulator's
#'   own, and the draws are replayed in one process.
#' @param point_estimate,ndraws_point_estimate brms's arguments, which
#'   collapse the draws to their `"mean"` or `"median"` FIRST and run
#'   the method once at that one parameter vector, repeated
#'   `ndraws_point_estimate` times. It is a parameter-space operation
#'   and not a summary of the output, as in brms.
#' @param nlpar The parameter an `nlf()` body names. brms keeps it
#'   apart from `dpar`; frmtmb asks for either by the `dpar` name, so
#'   this is the same setting and the slot is here for brms's position.
#' @param incl_thres For `posterior_linpred()`: refused. brms subtracts
#'   a cumulative family's thresholds from the predictor; frmtmb
#'   returns the latent predictor itself.
#' @param negative_rt For `posterior_predict()`: refused. It is brms's
#'   sign convention for its own wiener family.
#' @param transform For `posterior_predict()`: a function applied to
#'   the finished draws, in brms's own fifth position.
#' @param ... Refused: an argument the method does not have is an error
#'   naming it, rather than a silently ignored name. The exceptions are
#'   brms's `allow_new_levels` (and `allow.new.levels`) and
#'   `sample_new_levels`. `allow_new_levels = FALSE`, and `TRUE` with
#'   levels the fit saw, answer as the call without it does. `TRUE`
#'   with a level the fit did not see, including a `newdata` that leaves
#'   the grouping column out, is refused: brms draws that level's
#'   effect from each posterior draw, which is not built here, and
#'   predicting it at the population level would drop the group
#'   variance from every draw. [frmtmb::predict.frmtmb_fit()] predicts
#'   unseen levels from the maximum-likelihood fit.
#' @return A draws-by-observations matrix; for a categorical outcome
#'   `posterior_epred()` returns a draws-by-observations-by-categories
#'   array (see the section below).
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#' set.seed(9)
#' dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
#' dd$y <- rpois(80, exp(0.3 + 0.4 * dd$x + rnorm(8, 0, 0.5)[dd$g]))
#' fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
#' ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)
#'
#' nd <- data.frame(x = c(-1, 0, 1),
#'                  g = factor(1, levels = levels(dd$g)))
#'
#' # the expected response per draw: uncertainty in the mean
#' ep <- posterior_epred(ds, newdata = nd)
#' apply(ep, 2, quantile, c(0.025, 0.5, 0.975))
#'
#' # the predictive distribution adds the family's own noise, so its
#' # intervals are wider
#' pp <- posterior_predict(ds, newdata = nd)
#' apply(pp, 2, quantile, c(0.025, 0.5, 0.975))
#'
#' # the linear predictor itself, on the link scale by default
#' head(posterior_linpred(ds, newdata = nd, ndraws = 5))
#' }
#' }
#' @export
posterior_epred <- function(object, ...) UseMethod("posterior_epred")

#' @rdname posterior_epred
#' @exportS3Method rstantools::posterior_epred
#' @export
posterior_epred.frmtmb_draws <- function(object, newdata = NULL,
                                         re_formula = arg_unset(),
                                         re.form = arg_unset(),
                                         resp = NULL, dpar = NULL,
                                         nlpar = NULL, ndraws = NULL,
                                         draw_ids = NULL, sort = FALSE,
                                         point_estimate = NULL,
                                         ndraws_point_estimate = 1, ...) {
  # unknown names are refused rather than swallowed (item 2.5e); the
  # new-level spellings pass to draws_refuse_new_levels(), which answers
  # them. Forwarding allow_new_levels to the per-draw predictor, as this
  # method did briefly, put an unseen level at a group effect of exactly
  # 0: bit-identical to re_formula = NA, the group sd (1.02 there)
  # missing from every draw. See dev/shapes-findings.md.
  frm_check_dots(..., .allow = draws_new_level_args)
  re_form <- re_form_arg(re_formula, re.form, "posterior_epred()")
  draws_refuse_new_levels(object, newdata, list(...), "posterior_epred()")
  draws_refuse_sort(sort, "posterior_epred()")
  object <- draws_at_point_estimate(object, point_estimate,
                                    ndraws_point_estimate)
  dpar <- draws_dpar_arg(dpar, nlpar, "posterior_epred()")
  idx <- draws_par_index(object$fit)
  rows <- draws_subsample(object, ndraws, draw_ids)
  out <- NULL
  cat_out <- FALSE
  for (k in seq_along(rows)) {
    sh <- draws_fit_at(object, rows[k], idx)
    # frm_linpred(), not predict(): predict() is brms's predictive
    # summary in frmtmb's development version, and what one draw
    # contributes here is the expected response at its parameters
    p <- frm_linpred(sh, newdata = newdata, resp = resp, dpar = dpar,
                     re_formula = re_form, type = "response")
    if (is.null(out)) {
      # A categorical outcome predicts a matrix per draw (an ordinal
      # family's n x K category probabilities), so the draws stack into
      # a draws x observations x categories array rather than
      # flattening the category margin into column names: that is
      # brms's posterior_epred() return for polytomous families, and it
      # keeps ep[, , "cat"] and ep[k, , ] addressable directly.
      cat_out <- is.matrix(p)
      out <- if (cat_out) {
        array(NA_real_, c(length(rows), nrow(p), ncol(p)),
              dimnames = list(NULL, rownames(p), colnames(p)))
      } else {
        matrix(NA_real_, length(rows), length(p))
      }
    }
    if (cat_out) out[k, , ] <- p else out[k, ] <- p
  }
  out
}

#' @rdname posterior_epred
#' @export
posterior_linpred <- function(object, transform = FALSE, ...) {
  UseMethod("posterior_linpred")
}

#' @rdname posterior_epred
#' @param transform For `posterior_linpred()`: if `TRUE`, apply the
#'   inverse link (the value of the `mu` dpar on its natural scale,
#'   brms's convention; unlike `posterior_epred()` this is not the
#'   response mean for zero-inflated and similar families).
#' @param dpar Which distributional parameter to evaluate: its linear
#'   predictor for `posterior_linpred()`, its response-scale value for
#'   `posterior_epred()`. The default is the family's `mu`.
#' @exportS3Method rstantools::posterior_linpred
#' @export
posterior_linpred.frmtmb_draws <- function(object, transform = FALSE,
                                           newdata = NULL,
                                           re_formula = arg_unset(),
                                           re.form = arg_unset(),
                                           resp = NULL, dpar = NULL,
                                           nlpar = NULL,
                                           incl_thres = NULL,
                                           ndraws = NULL,
                                           draw_ids = NULL, sort = FALSE,
                                           point_estimate = NULL,
                                           ndraws_point_estimate = 1, ...) {
  re_form <- re_form_arg(re_formula, re.form, "posterior_linpred()")
  draws_refuse_new_levels(object, newdata, list(...), "posterior_linpred()")
  draws_refuse_sort(sort, "posterior_linpred()")
  object <- draws_at_point_estimate(object, point_estimate,
                                    ndraws_point_estimate)
  dpar <- draws_dpar_arg(dpar, nlpar, "posterior_linpred()")
  if (!is.null(incl_thres) && !identical(incl_thres, FALSE)) {
    frm_stop("posterior_linpred(incl_thres = TRUE) subtracts an ordinal ",
             "family's thresholds from the linear predictor, which brms ",
             "supports for cumulative families alone. frmtmb keeps the ",
             "thresholds out of the predictor: frm_linpred(type = \"link\") ",
             "and this function return the latent predictor itself, and ",
             "the thresholds are coefficients you can read off ",
             "posterior_summary()", call. = FALSE)
  }
  idx <- draws_par_index(object$fit)
  rows <- draws_subsample(object, ndraws, draw_ids)
  # This function is about ONE distributional parameter, so the dpar is
  # resolved here rather than left to predict()'s type dispatch: on an
  # ordinal fit `type = "response"` with no dpar is the whole category
  # distribution (posterior_epred()'s quantity), not the mu predictor
  # this promises.
  dpar <- dpar %||% draws_default_dpar(object$fit, resp)
  out <- NULL
  for (k in seq_along(rows)) {
    sh <- draws_fit_at(object, rows[k], idx)
    p <- frm_linpred(sh, newdata = newdata, resp = resp, dpar = dpar,
                     re_formula = re_form,
                     type = if (transform) "response" else "link")
    if (is.null(out)) out <- matrix(NA_real_, length(rows), length(p))
    out[k, ] <- p
  }
  out
}

#' brms's `dpar` and `nlpar` are ONE argument here.
#'
#' brms separates the distributional parameters of a family from the
#' parameters an `nlf()` body names; frmtmb asks for either through
#' `dpar`, because `predict.frmtmb_fit()` resolves a non-linear
#' parameter by that name too. The `nlpar` slot is carried so that a
#' positional brms call lands where brms puts it, and it feeds the same
#' setting.
#'
#' @noRd
draws_dpar_arg <- function(dpar, nlpar, what) {
  if (is.null(nlpar)) return(dpar)
  if (!is.null(dpar)) {
    frm_stop(what, " was given both `dpar` and `nlpar`. frmtmb asks for a ",
             "non-linear parameter by the same `dpar` name a family's own ",
             "parameters use, so the two are one setting here. Pass one of ",
             "them", call. = FALSE)
  }
  nlpar
}

#' The dpar `predict()` defaults to for one response: `mu` when the
#' family has it, the first primary dpar otherwise (the resolution in
#' `predict.frmtmb_fit()`, kept in one place).
#'
#' @noRd
draws_default_dpar <- function(fit, resp) {
  rs <- fit$spec$responses[[resp %||% names(fit$spec$responses)[1L]]]
  if (is.null(rs)) return(NULL)
  if ("mu" %in% names(rs$dpars)) "mu" else rs$primary_dpars[1L]
}

#' @rdname posterior_epred
#' @export
posterior_predict <- function(object, ...) UseMethod("posterior_predict")

#' @rdname posterior_epred
#' @exportS3Method rstantools::posterior_predict
#' @export
posterior_predict.frmtmb_draws <- function(object, newdata = NULL,
                                           re_formula = arg_unset(),
                                           re.form = arg_unset(),
                                           transform = NULL,
                                           resp = NULL,
                                           negative_rt = FALSE,
                                           ndraws = NULL,
                                           draw_ids = NULL, sort = FALSE,
                                           ntrys = NULL, cores = NULL,
                                           point_estimate = NULL,
                                           ndraws_point_estimate = 1, ...) {
  frm_check_dots(..., .allow = draws_new_level_args)
  re_form <- re_form_arg(re_formula, re.form, "posterior_predict()")
  draws_refuse_new_levels(object, newdata, list(...), "posterior_predict()")
  draws_refuse_sort(sort, "posterior_predict()")
  draws_refuse_ntrys_cores(ntrys, cores, "posterior_predict()")
  object <- draws_at_point_estimate(object, point_estimate,
                                    ndraws_point_estimate)
  check_flag(negative_rt, "negative_rt")
  if (negative_rt) {
    frm_stop("posterior_predict(negative_rt = TRUE) is brms's sign ",
             "convention for its wiener family, which codes the lower ",
             "boundary as a negative reaction time. frmtmb's ",
             "evidence-accumulation families return the response and the ",
             "time as they declare them; see the family's own ",
             "documentation in frmtmb.eam", call. = FALSE)
  }
  fit <- object$fit
  resp <- resp %||% names(fit$spec$responses)[1L]
  rspec <- fit$spec$responses[[resp]]
  if (!sim_can(rspec$family)) {
    frm_stop("posterior_predict(): family '", rspec$family[["family"]],
             "' has no simulator yet", sim_note(rspec$family), call. = FALSE)
  }
  idx <- draws_par_index(object$fit)
  rows <- draws_subsample(object, ndraws, draw_ids)
  av <- if (is.null(newdata)) {
    fit$frame[["aterm_values"]][[resp]]
  } else if (has_trunc(rspec)) {
    # truncation bounds must follow the newdata rows, or the draws land
    # outside the support the likelihood was normalized on
    aterms_for_newdata(rspec, newdata)
  } else {
    list()
  }
  if (!is.null(newdata) &&
      sim_is_structured(sim_context(fit, rspec, list(), aterms = av))) {
    # the sequence, group and residual-correlation structures a
    # structured draw walks were built from the TRAINING rows and index
    # them; newdata rows appear in none of them
    frm_stop("posterior_predict(newdata =) is not supported for this ",
             "model: its draws are structured (a hidden state sequence, a ",
             "group-level latent class, or a correlated residual) and that ",
             "structure indexes the rows the model was fitted on. Drop ",
             "newdata to predict those rows", call. = FALSE)
  }
  if (!is.null(re_form) &&
      sim_is_structured(sim_context(fit, rspec, list(), aterms = av))) {
    # same reason from the other side: the structured draw IS a walk
    # over the fitted structure, so there is no "with the group effects
    # removed" version of it to hand back
    frm_stop("posterior_predict(re_formula =) is not supported for this ",
             "model: its draws are structured, and the structure IS the ",
             "group-level content a re_formula would remove. Drop the ",
             "argument to draw from the fitted structure", call. = FALSE)
  }
  out <- NULL
  arr <- FALSE
  for (k in seq_along(rows)) {
    sh <- draws_fit_at(object, rows[k], idx)
    dp <- if (is.null(newdata) && is.null(re_form)) {
      # the fast path IS the default: NULL keeps every random effect,
      # which is what a per-draw eval of the full model gives. A set
      # re_formula routes through predict() on the training rows, the
      # same as the newdata branch, so NA and one-sided formulas mean
      # here exactly what they mean there
      eval_dpars(sh)[[resp]]
    } else {
      dpv <- list()
      for (dnm in names(rspec$dpars)) {
        dpv[[dnm]] <- as.vector(frm_linpred(sh, newdata = newdata,
                                            dpar = dnm, resp = resp,
                                            re_formula = re_form,
                                            type = "response"))
      }
      dpv
    }
    ys <- sim_draw(sim_context(sh, rspec, dp, aterms = av,
                               n = length(dp[[1L]]),
                               extra = fit_extras(sh)))
    if (is.null(out)) {
      # a matrix-valued response (multinomial counts, mixture_mvn draws,
      # lca item codes) gives a ROW per observation, so the draws stack
      # into a draws x observations x columns array, the shape
      # posterior_epred() already uses for a category distribution
      arr <- is.matrix(ys)
      out <- if (arr) {
        array(NA_real_, c(length(rows), nrow(ys), ncol(ys)),
              dimnames = list(NULL, NULL, colnames(ys)))
      } else {
        matrix(NA_real_, length(rows), length(ys))
      }
    }
    if (arr) out[k, , ] <- ys else out[k, ] <- ys
  }
  # brms's fifth positional slot: a function applied to the finished
  # draws. It is carried so that a positional brms call means the same
  # thing here, and it does the same thing brms does with it
  if (!is.null(transform)) out <- match.fun(transform)(out)
  out
}

#' The one retired spelling `pp_check()` has to refuse itself.
#'
#' Its dots go to bayesplot, which accepts any name, so without this the
#' spelling that used to select the population check would go on
#' selecting nothing.
#'
#' @noRd
pp_check_retired_draws <- c(
  re.form = paste("lme4's spelling, no longer accepted. brms LEAKS it:",
                  "pp_check.brmsfit() forwards its dots to",
                  "posterior_predict(), so brms honors `re.form` there",
                  "while WARNING that it ignored it. Warn-then-honor is",
                  "the failure this refusal exists to stop, so the",
                  "spelling is refused here rather than copied. Pass",
                  "re_formula =")
)

#' @exportS3Method bayesplot::pp_check
#' @export
pp_check.frmtmb_draws <- function(object, type, ndraws = NULL,
                                  prefix = c("ppc", "ppd"), group = NULL,
                                  x = NULL, newdata = NULL, resp = NULL,
                                  draw_ids = NULL, nsamples = NULL,
                                  subset = NULL, ..., re_formula = NULL) {
  # brms:::pp_check.brmsfit's slots in its order, so pp_check(ds,
  # "dens_overlay", 20, "ppd") asks for a predictive-distribution plot
  # as in brms; `re_formula` used to sit where brms has `prefix`.
  #
  # bayesplot's function takes dots of its own, so an unknown name there
  # is its business. The retired lme4 spelling is not: it named a real
  # setting until the rename, so it is refused here rather than passed on
  # to be ignored.
  frm_check_dots(..., .allow = TRUE, .unsupported = pp_check_retired_draws)
  if (missing(type)) type <- "dens_overlay"
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    frm_stop("pp_check(): `type` must be a single string", call. = FALSE)
  }
  prefix <- frm_match_arg(prefix)
  ndraws_given <- !missing(ndraws) || !missing(nsamples)
  if (!is.null(nsamples)) {
    frm_warning("Argument 'nsamples' is deprecated. Please use argument ",
                "'ndraws' instead.", call. = FALSE)
    ndraws <- nsamples
  }
  if (!is.null(subset)) {
    frm_warning("Argument 'subset' is deprecated. Please use argument ",
                "'draw_ids' instead.", call. = FALSE)
    draw_ids <- subset
  }
  fun <- draws_bayesplot_fun(paste0(prefix, "_", type), "pp_check(type =)")
  fit <- draws_base_fit(object)
  rspec <- if (is.null(resp)) single_response(fit, "pp_check()") else
    fit$spec$responses[[resp]]
  if (is.null(rspec)) {
    frm_stop("pp_check(resp = \"", resp, "\") names no response of this ",
             "model; it has ", paste(names(fit$spec$responses),
                                      collapse = ", "), call. = FALSE)
  }
  resp <- rspec$resp_name
  data <- newdata %||% fit$frame[["data_frame"]]
  fargs <- names(formals(fun))
  if ("group" %in% fargs) {
    if (is.null(group)) {
      frm_stop("Argument 'group' is required for ppc type '", type, "'.",
               call. = FALSE)
    }
  }
  for (v in c(group, x)) {
    if (!is.character(v) || length(v) != 1L || !v %in% names(data)) {
      frm_stop("Variable '", v, "' could not be found in the data.",
               call. = FALSE)
    }
  }
  if (!ndraws_given) {
    # brms's own defaults, with its own messages: every draw for the
    # types that average over draws, ten otherwise
    aps_types <- c("error_scatter_avg", "error_scatter_avg_vs_x",
                   "intervals", "intervals_grouped", "loo_intervals",
                   "loo_pit", "loo_pit_overlay", "loo_pit_qq",
                   "loo_ribbon", "loo_pit_ecdf", "pit_ecdf",
                   "pit_ecdf_grouped", "ribbon", "ribbon_grouped",
                   "rootogram", "scatter_avg", "scatter_avg_grouped",
                   "stat", "stat_2d", "stat_freqpoly_grouped",
                   "stat_grouped", "violin_grouped")
    if (!is.null(draw_ids)) {
      ndraws <- NULL
    } else if (type %in% aps_types) {
      ndraws <- NULL
      frm_message("Using all posterior draws for ppc type '", type,
                  "' by default.")
    } else {
      ndraws <- 10
      frm_message("Using 10 posterior draws for ppc type '", type,
                  "' by default.")
    }
  }
  pred <- if (identical(type, "error_binned")) posterior_epred else
    posterior_predict
  yrep <- pred(object, newdata = newdata, resp = resp, ndraws = ndraws,
               draw_ids = draw_ids, re_formula = re_formula)
  if (length(dim(yrep)) > 2L) {
    frm_stop("pp_check() on draws supports vector responses", call. = FALSE)
  }
  args <- list()
  take <- rep(TRUE, ncol(yrep))
  if (prefix == "ppc") {
    y <- as.numeric(draws_response_values(fit, resp, newdata,
                                          "pp_check()"))
    if (anyNA(y)) {
      frm_warning("NA responses are not shown in 'pp_check'.", call. = FALSE)
      take <- !is.na(y)
    }
    args$y <- y[take]
    args$yrep <- yrep[, take, drop = FALSE]
  } else {
    args$ypred <- yrep
  }
  if (!is.null(group)) args$group <- data[[group]][take]
  if (!is.null(x)) {
    xv <- data[[x]][take]
    args$x <- if (is.factor(xv) || is.character(xv) || is.logical(xv)) xv else
      as.numeric(xv)
  }
  do.call(fun, c(args, list(...)))
}

#' Convert draws to a posterior draws object
#'
#' brms's converters, with brms's arguments and brms's output. The
#' `as_draws_*()` family takes `variable` (exact names unless `regex`)
#' and returns a posterior draws object; `as_draws()` is brms's
#' `as_draws_list()`. `as.matrix()`, `as.array()` and `as.data.frame()`
#' return brms's unclassed objects, with brms's deprecated `pars` and
#' `subset` accepted under brms's own warning.
#'
#' @param x A `frmtmb_draws` object.
#' @param variable Variables to keep, by exact name unless `regex`.
#' @param regex If `TRUE`, `variable` is a regular expression.
#' @param inc_warmup Only `FALSE`: the draws matrix keeps the post-warmup
#'   draws alone.
#' @param pars brms's deprecated alias of `variable`, a regular
#'   expression; it warns, as brms's does.
#' @param draw Draws to keep, by index.
#' @param subset brms's deprecated alias of `draw`; it warns.
#' @param row.names,optional Accepted for the generic and unused, as in
#'   brms.
#' @param ... For `as.matrix()`, `as.array()` and `as.data.frame()`, the
#'   `regex`, `fixed` and `inc_warmup` brms passes on; anything else is
#'   refused by name, rather than silently changing nothing.
#' @return A `posterior::draws_matrix`: one column per sampled variable
#'   and one row per draw.
#' @examples
#' \donttest{
#' if (requireNamespace("posterior", quietly = TRUE) &&
#'     requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
#'   dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
#'   fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'   ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)
#'
#'   # hands the draws to the posterior package, keeping the frmtmb
#'   # parameter names
#'   dm <- as_draws(ds)
#'   posterior::summarise_draws(dm)
#'   # which is what variables() lists
#'   head(variables(ds))
#' }
#' }
#' @name sample-as_draws
NULL

#' @rdname sample-as_draws
#' @exportS3Method posterior::as_draws
#' @export
as_draws.frmtmb_draws <- function(x, variable = NULL, regex = FALSE,
                                  inc_warmup = FALSE, ...) {
  # brms's as_draws() is its as_draws_list()
  frm_check_dots(...)
  posterior::as_draws_list(draws_as_array(x, variable, regex, inc_warmup,
                                          "as_draws()"))
}

#' @exportS3Method posterior::variables
#' @export
variables.frmtmb_draws <- function(x, ...) {
  frm_check_dots(...)
  colnames(x$draws)
}

#' @rdname sample-as_draws
#' @export
as.data.frame.frmtmb_draws <- function(x, row.names = NULL, optional = TRUE,
                                       pars = NA, variable = NULL,
                                       draw = NULL, subset = NULL, ...) {
  a <- draws_accessor_args(x, pars, variable, draw, subset,
                           "as.data.frame()",
                           format = posterior::as_draws_df, ...)
  out <- posterior::as_draws_df(a)
  out <- as.data.frame(out)
  out$.chain <- out$.iteration <- out$.draw <- NULL
  out
}

# ---- the draws matrix in other shapes --------------------------------

#' @rdname sample-as_draws
#' @export
as.array.frmtmb_draws <- function(x, pars = NA, variable = NULL,
                                  draw = NULL, subset = NULL, ...) {
  unclass(draws_accessor_args(x, pars, variable, draw, subset,
                              "as.array()", ...))
}

#' @rdname sample-as_draws
#' @exportS3Method posterior::as_draws_matrix
#' @export
as_draws_matrix.frmtmb_draws <- function(x, variable = NULL, regex = FALSE,
                                         inc_warmup = FALSE, ...) {
  frm_check_dots(...)
  posterior::as_draws_matrix(draws_as_array(x, variable, regex, inc_warmup,
                                            "as_draws_matrix()"))
}

#' @rdname sample-as_draws
#' @exportS3Method posterior::as_draws_array
#' @export
as_draws_array.frmtmb_draws <- function(x, variable = NULL, regex = FALSE,
                                        inc_warmup = FALSE, ...) {
  frm_check_dots(...)
  draws_as_array(x, variable, regex, inc_warmup, "as_draws_array()")
}

#' @rdname sample-as_draws
#' @exportS3Method posterior::as_draws_df
#' @export
as_draws_df.frmtmb_draws <- function(x, variable = NULL, regex = FALSE,
                                     inc_warmup = FALSE, ...) {
  frm_check_dots(...)
  posterior::as_draws_df(draws_as_array(x, variable, regex, inc_warmup,
                                        "as_draws_df()"))
}

#' @rdname sample-as_draws
#' @exportS3Method posterior::as_draws_list
#' @export
as_draws_list.frmtmb_draws <- function(x, variable = NULL, regex = FALSE,
                                       inc_warmup = FALSE, ...) {
  frm_check_dots(...)
  posterior::as_draws_list(draws_as_array(x, variable, regex, inc_warmup,
                                          "as_draws_list()"))
}

#' @rdname sample-as_draws
#' @exportS3Method posterior::as_draws_rvars
#' @export
as_draws_rvars.frmtmb_draws <- function(x, variable = NULL, regex = FALSE,
                                        inc_warmup = FALSE, ...) {
  frm_check_dots(...)
  posterior::as_draws_rvars(draws_as_array(x, variable, regex, inc_warmup,
                                           "as_draws_rvars()"))
}

#' @rdname sample-as_draws
#' @export
as.mcmc <- function(x, ...) UseMethod("as.mcmc")

#' @rdname sample-as_draws
#' @param pars Variables to keep, in brms's spelling: `NA` (the
#'   default) for all of them, otherwise a character vector matched as
#'   a regular expression unless `fixed = TRUE`. The argument sits in
#'   brms's own second position, so `as.mcmc(x, TRUE)` is refused here
#'   exactly as brms refuses it.
#' @param fixed If `TRUE`, `pars` is matched by exact name.
#' @param combine_chains If `TRUE`, one `mcmc` object over the pooled
#'   draws; otherwise an `mcmc.list` with one component per chain, which
#'   is what coda's diagnostics (`gelman.diag()`) need.
#' @param inc_warmup Accepted for brms's signature and only `FALSE` is
#'   supported: a `frmtmb_draws` keeps the post-warmup draws alone.
#' @exportS3Method coda::as.mcmc
#' @export
as.mcmc.frmtmb_draws <- function(x, pars = NA, fixed = FALSE,
                                 combine_chains = FALSE,
                                 inc_warmup = FALSE, ...) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    frm_stop("as.mcmc() needs the 'coda' package; as_draws() and ",
             "as.array() give the same draws without it", call. = FALSE)
  }
  check_flag(combine_chains, "combine_chains")
  check_flag(inc_warmup, "inc_warmup")
  if (inc_warmup) {
    frm_stop("as.mcmc(inc_warmup = TRUE) has nothing to include: the draws ",
             "matrix holds the post-warmup draws only, which is what ",
             "frm_sample() keeps. The warmup, if the sampler saved it, is ",
             "in `x$stanfit`", call. = FALSE)
  }
  sel <- draws_select_variables(x, pars, NULL, FALSE, fixed, "as.mcmc()")
  if (combine_chains) {
    m <- x$draws
    if (!is.null(sel)) m <- m[, sel, drop = FALSE]
    return(coda::as.mcmc(m))
  }
  a <- draws_raw_array(x)
  if (!is.null(sel)) a <- a[, , sel, drop = FALSE]
  dn <- list(NULL, dimnames(a)[[3L]])
  coda::as.mcmc.list(lapply(seq_len(dim(a)[2L]), function(ch) {
    coda::as.mcmc(array(a[, ch, ], dim(a)[c(1L, 3L)], dimnames = dn))
  }))
}

# ---- draws-matrix dimensions -----------------------------------------

#' Size of a draws object
#'
#' `ndraws()` counts the post-warmup draws (all chains pooled),
#' `niterations()` the draws per chain, `nchains()` the chains and
#' `nvariables()` the sampled parameters. The names and meanings are
#' posterior's; frmtmb registers methods with posterior so that the
#' generics work whether or not that package is attached.
#'
#' @param x A `frmtmb_draws` from [frm_sample()].
#' @param ... Unused. `nvariables()` carries it because
#'   posterior's generic does.
#' @return A single integer.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#'   ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
#'                    data = dd, chains = 1, iter = 500, refresh = 0)
#'   c(ndraws(ds), niterations(ds), nchains(ds), nvariables(ds))
#' }
#' }
#' @name draws-dimensions
NULL

# posterior's nchains(), ndraws() and niterations() generics take x
# alone and nvariables() takes (x, ...); dev/generics-audit2.R in
# core measured it. The signatures have to agree, because core now
# hands these generics back to posterior.
#' @rdname draws-dimensions
#' @exportS3Method posterior::ndraws
#' @export
ndraws.frmtmb_draws <- function(x) nrow(x$draws)

#' @rdname draws-dimensions
#' @exportS3Method posterior::nchains
#' @export
nchains.frmtmb_draws <- function(x) {
  as.integer(x$stanfit@sim$chains %||% 1L)
}

#' @rdname draws-dimensions
#' @exportS3Method posterior::niterations
#' @export
niterations.frmtmb_draws <- function(x) {
  as.integer(nrow(x$draws) %/% nchains(x))
}

#' @rdname draws-dimensions
#' @exportS3Method posterior::nvariables
#' @export
nvariables.frmtmb_draws <- function(x, ...) {
  frm_check_dots(...)
  ncol(x$draws)
}

# ---- posterior summaries ---------------------------------------------

#' Summaries and intervals of draws
#'
#' `posterior_summary()` reduces draws to estimate, error and quantiles
#' in brms's column layout (`Estimate`, `Est.Error`, `Q2.5`, `Q97.5`);
#' `posterior_interval()` gives the central interval alone, in
#' rstantools' layout. Both work on a `frmtmb_draws` object and on any
#' matrix of draws, which is what makes
#' `posterior_summary(bayes_R2(ds, summary = FALSE))` work.
#'
#' `predictive_interval()` is the same central interval of
#' [posterior_predict()] draws, and `predictive_error()` is the matrix
#' of predictive residuals `y - yrep`, one row per draw.
#'
#' @param object A `frmtmb_draws`, or a matrix of draws
#'   (variables in columns).
#' @param x The same, for `posterior_summary()`, whose generic
#'   is brms's and names its first argument `x`.
#' @param probs Quantiles for `posterior_summary()`.
#' @param prob Central interval width for `posterior_interval()` and
#'   `predictive_interval()`.
#' @param robust If `TRUE`, median and MAD instead of mean and SD.
#' @param variable Optional subset of variables, by name.
#' @param pars brms's alias of `variable`, in brms's own second
#'   position on `posterior_interval()`: `NA` (the default) for every
#'   variable, otherwise a character vector matched as a regular
#'   expression unless `fixed = TRUE`. brms refuses a `pars` that is
#'   neither `NA` nor character, and so does this, which is why
#'   `posterior_interval(x, 0.9)` is a refusal and not an interval.
#' @param regex If `TRUE`, `variable` is a regular expression.
#' @param fixed If `TRUE`, `pars` is matched by exact name.
#' @param method For `predictive_error()`, which predictive draws the
#'   error is taken against: `"posterior_predict"` (the default) or
#'   `"posterior_epred"`.
#' @param ndraws,draw_ids,newdata,resp Passed to
#'   [posterior_predict()]. `predictive_error(newdata =)` re-evaluates
#'   the response term on `newdata`, so `newdata` must carry the
#'   response.
#' @param re_formula,re.form Passed to [posterior_predict()], which
#'   takes brms's `re_formula` and accepts lme4's `re.form` as an alias
#'   of it. Pass one or the other; see the *Argument spellings* section
#'   of [posterior_epred()].
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return A matrix with one row per variable (or per observation, for
#'   the predictive functions), except `predictive_error()`, which
#'   returns a draws-by-observations matrix.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#'   ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
#'                    data = dd, chains = 1, iter = 500, refresh = 0)
#'   posterior_summary(ds, variable = c("b_Intercept", "b_x"))
#'   posterior_interval(ds, prob = 0.9, variable = "b_x")
#'   head(predictive_interval(ds))
#' }
#' }
#' @name sample-posterior_summary
NULL


#' @rdname sample-posterior_summary
#' @exportS3Method brms::posterior_summary
#' @export
posterior_summary.frmtmb_draws <- function(x, pars = NA, variable = NULL,
                                           probs = c(0.025, 0.975),
                                           robust = FALSE, ...) {
  # brms:::posterior_summary.brmsfit: every variable unless told
  # otherwise, including the group-level coefficients and lp__
  out <- as.matrix(x, pars = pars, variable = variable, ...)
  posterior_summary(out, probs = probs, robust = robust)
}


#' brms's `pars` argument, resolved to variable names.
#'
#' brms's `extract_pars()` takes `NA` for "every variable" and
#' otherwise a character vector matched as a regular expression against
#' the variable names, or by exact name when `fixed = TRUE`. Anything
#' else brms refuses with the message repeated here, which is why
#' `posterior_interval(ds, 0.9)` and `as.mcmc(ds, TRUE)` are refusals
#' and not answers: in brms that second position is `pars`.
#'
#' @noRd
draws_extract_pars <- function(pars, all_pars, fixed = FALSE) {
  if (!(anyNA(pars) || is.character(pars))) {
    frm_stop("Argument 'pars' must be NA or a character vector.",
             call. = FALSE)
  }
  if (anyNA(pars)) return(all_pars)
  check_flag(fixed, "fixed")
  if (fixed) return(intersect(pars, all_pars))
  unique(unlist(lapply(pars, function(p) all_pars[grepl(p, all_pars)])))
}

#' The variables a brms-style `pars`/`variable` pair selects, or `NULL`
#' for "the calling method's own default".
#'
#' brms's rule, which this follows: `variable` is exact unless `regex`,
#' `pars` is a regular expression unless `fixed`, and `pars` wins when
#' both are given because it is the older spelling of the same thing.
#'
#' @noRd
draws_select_variables <- function(x, pars = NA, variable = NULL,
                                   regex = FALSE, fixed = FALSE,
                                   what = "this function") {
  all_pars <- colnames(x$draws)
  if (!(is.logical(pars) && length(pars) == 1L && is.na(pars))) {
    return(draws_extract_pars(pars, all_pars, fixed = fixed))
  }
  if (is.null(variable)) return(NULL)
  if (!is.character(variable)) {
    frm_stop(what, ": `variable` names variables and must be a character ",
             "vector; variables(x) lists what is there", call. = FALSE)
  }
  check_flag(regex, "regex")
  if (regex) return(draws_extract_pars(variable, all_pars))
  miss <- setdiff(variable, all_pars)
  if (length(miss)) {
    frm_stop(what, ": variable = names ", paste(miss, collapse = ", "),
             ", which the draws do not contain. variables() lists what is ",
             "there; the draws use brms's names (b_Intercept, not ",
             "(Intercept))", call. = FALSE)
  }
  variable
}

#' @rdname sample-posterior_summary
#' @export
posterior_interval <- function(object, ...) UseMethod("posterior_interval")

#' @rdname sample-posterior_summary
#' @exportS3Method rstantools::posterior_interval
#' @export
posterior_interval.frmtmb_draws <- function(object, pars = NA,
                                            variable = NULL,
                                            prob = 0.95, regex = FALSE,
                                            fixed = FALSE, ...) {
  frm_check_dots(...)
  check_probability(prob, "prob")
  ps <- as.matrix(object, pars = pars, variable = variable, regex = regex,
                  fixed = fixed)
  # rstantools' .central_intervals(), which brms reaches: every
  # variable by default, and plain dimnames
  a <- (1 - prob) / 2
  probs <- c(a, 1 - a)
  out <- t(apply(ps, 2L, stats::quantile, probs = probs))
  structure(out, dimnames = list(colnames(ps), paste0(100 * probs, "%")))
}

#' @rdname sample-posterior_summary
#' @export
predictive_interval <- function(object, ...) UseMethod("predictive_interval")

#' @rdname sample-posterior_summary
#' @exportS3Method rstantools::predictive_interval
#' @export
predictive_interval.frmtmb_draws <- function(object, prob = 0.9,
                                             newdata = NULL, resp = NULL,
                                             re_formula = arg_unset(),
                                             re.form = arg_unset(),
                                             ndraws = NULL, ...) {
  # Both spellings, on the strength of what brms DOES rather than what
  # it declares. `predictive_interval.brmsfit` declares neither, but its
  # whole body is `posterior_predict(object, ...)` and
  # `posterior_predict.brmsfit` has a `re.form` formal, so brms accepts
  # and honors `predictive_interval(x, re.form = NA)` one frame down.
  frm_check_dots(...)
  re_form <- re_form_arg(re_formula, re.form, "predictive_interval()")
  yrep <- posterior_predict(object, newdata = newdata, resp = resp,
                            re_formula = re_form, ndraws = ndraws)
  if (length(dim(yrep)) > 2L) {
    frm_stop("predictive_interval() needs one predicted number per ",
             "observation, and this model's draws are a matrix per ",
             "observation (multinomial counts, mixture_mvn draws or lca ",
             "item codes). Take the interval of the column you want from ",
             "posterior_predict() yourself", call. = FALSE)
  }
  a <- (1 - prob) / 2
  t(apply(yrep, 2L, stats::quantile, probs = c(a, 1 - a)))
}

#' @rdname sample-posterior_summary
#' @export
predictive_error <- function(object, ...) UseMethod("predictive_error")

#' @rdname sample-posterior_summary
#' @exportS3Method rstantools::predictive_error
#' @export
predictive_error.frmtmb_draws <- function(object, newdata = NULL,
                                          re_formula = arg_unset(),
                                          re.form = arg_unset(),
                                          method = "posterior_predict",
                                          resp = NULL, ndraws = NULL,
                                          draw_ids = NULL, ...) {
  re_form <- re_form_arg(re_formula, re.form, "predictive_error()")
  method <- frm_match_arg(method,
                          c("posterior_predict", "posterior_epred"))
  fit <- draws_base_fit(object)
  resp <- resp %||% names(fit$spec$responses)[1L]
  y <- draws_response_values(fit, resp, newdata, "predictive_error()")
  if (is.matrix(y)) {
    frm_stop("predictive_error() needs a vector response; this one is a ",
             "matrix (multinomial counts, mixture_mvn columns or lca ",
             "items), and 'the' error of a row of counts is not defined. ",
             "Subtract the column you want from posterior_predict() ",
             "yourself", call. = FALSE)
  }
  # brms's `method`: the predictive draws the error is taken against,
  # either the predictive distribution or the expectation
  yrep <- if (identical(method, "posterior_epred")) {
    posterior_epred(object, newdata = newdata, re_formula = re_form,
                    resp = resp, ndraws = ndraws, draw_ids = draw_ids)
  } else {
    posterior_predict(object, newdata = newdata, re_formula = re_form,
                      resp = resp, ndraws = ndraws,
                      draw_ids = draw_ids)
  }
  # brms's convention: the error is y - yrep, one row per draw
  sweep(-yrep, 2L, as.numeric(y), "+")
}

#' The observed response the predictive error is taken against: the
#' fitted rows, or `newdata`'s own column when one is given.
#'
#' The response TERM is re-evaluated on `newdata`, not looked up by
#' column name, so a transformed response (`log(y) ~ x`) is handled the
#' same way the model frame handled it.
#'
#' @noRd
draws_response_values <- function(fit, resp, newdata, what) {
  if (is.null(newdata)) return(fit$frame[["y"]][[resp]])
  rspec <- fit$spec$responses[[resp]]
  y <- tryCatch(eval(rspec$resp_expr, newdata, rspec$formula_env),
                error = function(e) NULL)
  if (is.null(y) || length(y) != nrow(newdata)) {
    frm_stop(what, " needs the observed response to subtract from, and ",
             "newdata does not supply '", deparse1(rspec$resp_expr),
             "' for its ", nrow(newdata), " rows. Add the response column ",
             "to newdata, or call posterior_predict(newdata =) and ",
             "subtract your own", call. = FALSE)
  }
  y
}

# ---- structural delegations to the originating fit -------------------
#
# These read only the model's STRUCTURE, never its estimates, so they go
# straight to the fit even when it was assembled by the formula route
# and has no maximum-likelihood estimate behind it (draws_base_fit()).

#' Model structure behind a set of draws
#'
#' `nobs()`, `formula()`, `family()`, `getCall()` and `ngrps()` on a
#' `frmtmb_draws` report the model the sampler ran, by delegating to the
#' fit stored inside it. They read structure only, so they work on draws
#' from the formula route, which has no maximum-likelihood estimate.
#'
#' `coef()`, `fixef()`, `ranef()` and `VarCorr()` are posterior
#' quantities, not structural ones, and they are brms's methods on these
#' draws, compared `identical()` against brms's own installed methods in
#' `dev/brmsnames-findings.md`. `fixef()` is a coefficients x statistics
#' matrix; `ranef()` and `coef()` are a list keyed by grouping factor of
#' `levels x statistics x coefficients` arrays, and `coef()` broadcasts
#' every population-level coefficient over the levels, as brms's does;
#' `VarCorr()` is a list keyed by grouping factor with `sd`, and `cor` and
#' `cov` when the group has correlations, then `residual__`. With
#' `summary = FALSE` each returns brms's raw draws instead: a draws x
#' coefficients matrix, a `draws x levels x coefficients` array, or a
#' `draws x coefficients x coefficients` array.
#'
#' `VarCorr()`'s standard deviations and correlations are computed per
#' draw from the sampled covariance parameters, because the sampler
#' stores `theta` and not brms's `sd_` and `cor_` draws.
#'
#' @param object,x A `frmtmb_draws` from [frm_sample()].
#' @param summary If `TRUE` (brms's default), summaries; otherwise the
#'   draws, as above.
#' @param robust If `TRUE`, median and MAD instead of mean and SD.
#' @param probs The quantiles to report.
#' @param pars For `fixef()` and `ranef()`, the coefficients to keep, by
#'   name without the `b_` prefix, as in brms.
#' @param groups For `ranef()`, the grouping factors to keep.
#' @param sigma Ignored, as in brms.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return As for the corresponding `frmtmb_fit` method.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#'   ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
#'                    data = dd, chains = 1, iter = 500, refresh = 0)
#'   nobs(ds)
#'   ngrps(ds)
#'   coef(ds)$g[1:3, , "Intercept"]
#'   dim(ranef(ds, summary = FALSE)$g)
#'   VarCorr(ds)$g$sd
#' }
#' }
#' @name draws-structure
NULL

#' @rdname draws-structure
#' @export
nobs.frmtmb_draws <- function(object, ...) {
  frm_check_dots(...)
  stats::nobs(draws_base_fit(object))
}

#' @rdname draws-structure
#' @export
formula.frmtmb_draws <- function(x, ...) {
  frm_check_dots(...)
  stats::formula(draws_base_fit(x))
}

#' @rdname draws-structure
#' @export
family.frmtmb_draws <- function(object, ...) {
  frm_check_dots(...)
  stats::family(draws_base_fit(object))
}

#' @rdname draws-structure
#' @export
getCall.frmtmb_draws <- function(x, ...) {
  frm_check_dots(...)
  x$fit$call
}

#' @exportS3Method brms::ngrps
#' @rawNamespace S3method(lme4::ngrps,frmtmb_draws)
#' @export
ngrps.frmtmb_draws <- function(object, ...) {
  frm_check_dots(...)
  ngrps(draws_base_fit(object))
}

#' @rdname draws-structure
#' @export
coef.frmtmb_draws <- function(object, summary = TRUE, robust = FALSE,
                              probs = c(0.025, 0.975), ...) {
  # brms:::coef.brmsfit: fixef draws broadcast over each group's levels
  # with the group's own coefficient draws added, a coefficient the group
  # does not vary added as 0
  frm_check_dots(...)
  check_flag(summary, "summary")
  if (!length(draws_ranef_layout(object))) {
    frm_stop("No group-level effects detected. Call method 'fixef' to ",
             "access population-level effects.", call. = FALSE)
  }
  fe <- fixef(object, summary = FALSE)
  co <- ranef(object, summary = FALSE)
  all_ranef_names <- unique(unlist(lapply(co, function(a) dimnames(a)[[3L]])))
  fixef_names <- colnames(fe)
  no_digits <- function(v) regmatches(v, regexpr("^[^\\[]+", v))
  fixef_no_digits <- no_digits(fixef_names)
  miss_fixef <- setdiff(all_ranef_names, fixef_names)
  miss_fixef_no_digits <- no_digits(miss_fixef)
  new_fixef <- stats::setNames(vector("list", length(miss_fixef)),
                               miss_fixef)
  for (k in seq_along(miss_fixef)) {
    mk <- match(miss_fixef_no_digits[k], fixef_names)
    if (!is.na(mk)) {
      new_fixef[[k]] <- fe[, mk]
    } else if (!miss_fixef[k] %in% fixef_no_digits) {
      new_fixef[[k]] <- 0
    }
  }
  rm_fixef <- fixef_names %in% miss_fixef_no_digits
  fe <- fe[, !rm_fixef, drop = FALSE]
  fe <- do.call(cbind, c(list(fe), Filter(Negate(is.null), new_fixef)))
  for (g in names(co)) {
    ranef_names <- dimnames(co[[g]])[[3L]]
    ranef_no_digits <- no_digits(ranef_names)
    miss_ranef <- setdiff(fixef_names, ranef_names)
    miss_ranef_no_digits <- no_digits(miss_ranef)
    new_ranef <- stats::setNames(vector("list", length(miss_ranef)),
                                 miss_ranef)
    for (k in seq_along(miss_ranef)) {
      mr <- match(miss_ranef_no_digits[k], ranef_names)
      if (!is.na(mr)) {
        new_ranef[[k]] <- co[[g]][, , mr]
      } else if (!miss_ranef[k] %in% ranef_no_digits) {
        new_ranef[[k]] <- array(0, dim = dim(co[[g]])[1:2])
      }
    }
    rm_ranef <- ranef_names %in% miss_ranef_no_digits
    A <- co[[g]][, , !rm_ranef, drop = FALSE]
    add <- Filter(Negate(is.null), new_ranef)
    if (length(add)) {
      dn <- dimnames(A)
      A <- array(c(A, unlist(add)),
                 dim = c(dim(A)[1:2], dim(A)[3L] + length(add)),
                 dimnames = list(dn[[1L]], dn[[2L]],
                                 c(dn[[3L]], names(add))))
    }
    for (nm in dimnames(A)[[3L]]) {
      A[, , nm] <- fe[, nm] + A[, , nm]
    }
    if (summary) A <- posterior_summary(A, probs, robust)
    co[[g]] <- A
  }
  co
}


# ---- sampler diagnostics and plots -----------------------------------

#' Sampler diagnostics and MCMC plots
#'
#' `rhat()` and `neff_ratio()` are the convergence diagnostics brms
#' reports, computed by the posterior package on these draws:
#' `rhat()` is the rank-normalized split-R-hat and `neff_ratio()` is
#' `min(ess_bulk, ess_tail) / ndraws`. `nuts_params()` and
#' `log_posterior()` delegate to bayesplot's `stanfit` methods on the
#' `stanfit` inside the draws object, which is what brms does too, so
#' every `bayesplot::mcmc_nuts_*()` display works. `mcmc_plot()` is
#' brms's spelling for "call a bayesplot `mcmc_*` function on these
#' draws"; `pairs()` is `bayesplot::mcmc_pairs()`.
#'
#' All of these report brms's parameter names (`b_x`,
#' `r_g[1,Intercept]`), not Stan's `par[1]`, except `nuts_params()`,
#' whose rows are the sampler's own quantities and not model
#' parameters.
#'
#' @section Which R-hat this is:
#' `rhat()` follows brms, whose `rhat.brmsfit()` is
#' `posterior::summarise_draws(rhat = posterior::rhat)`. That is the
#' rank-normalized split-R-hat of Vehtari et al. (2021), the maximum of
#' the bulk and tail quantities, and it is NOT the classic split-R-hat
#' that `rstan::summary()` reports. The two disagree by about the size
#' of the excess over 1 that either of them reports, so `rhat(ds)` and
#' `ds$stanfit` do not agree and are not meant to.
#'
#' `summary(ds)` agrees with `rhat(ds)`, because it reports the same
#' three posterior quantities `brms:::summary.brmsfit()` reports, under
#' the same column names: `Rhat`, `Bulk_ESS` and `Tail_ESS`. The
#' sampler's own classic split-R-hat and `n_eff` are in
#' `rstan::summary(ds$stanfit)$summary` for anyone who wants them.
#'
#' @section Two different `pars` rules, both brms's:
#' brms spells two different selectors `pars`, and this page carries
#' both because it documents methods on either side of the line.
#'
#' `mcmc_plot()` takes brms's deprecated alias of `variable`: `NA` is
#' every variable, a string is a regular expression unless
#' `fixed = TRUE`, and anything that is neither `NA` nor character is
#' refused. `as.mcmc()` and `posterior_interval()` take the same one.
#'
#' `rhat()` and `neff_ratio()` do not. brms's `rhat.brmsfit()` passes
#' `variable = pars` straight to `as_draws_array()`, so `NULL` is every
#' variable, a string is an EXACT variable name, `regex = TRUE` makes
#' it a regular expression, and a name that is not there is an error.
#' These two follow that rule, which is why their default is `NULL` and
#' not `NA`.
#'
#' @param object,x A `frmtmb_draws` from [frm_sample()].
#' @param type The bayesplot function to call, without the `mcmc_`
#'   prefix (default `"intervals"`).
#' @param variable For `mcmc_plot()` and `pairs()`, the variables to
#'   use, by name; it defaults to everything except the group-level
#'   coefficients and `lp__` (and to four of those for `pairs()`).
#'   `rhat()` and `neff_ratio()` do not take it,
#'   because brms's do not: their selector is `pars`. Naming it on
#'   either of those two is silently ignored today; see `...`.
#' @param pars Which variables to report, in brms's spelling. The rule
#'   differs by method; see *Two different `pars` rules, both brms's*.
#' @param regex For `rhat()` and `neff_ratio()`, `TRUE` makes `pars` a
#'   regular expression; for `mcmc_plot()` and `pairs()`, it makes
#'   `variable` one.
#' @param fixed For `mcmc_plot()` and `pairs()`, `TRUE` matches `pars` by
#'   exact name rather than as a regular expression.
#' @param ... For `mcmc_plot()` and `pairs()`, passed to the bayesplot
#'   function; for `nuts_params()` and `log_posterior()`, passed to
#'   bayesplot's own `stanfit` method, so `nuts_params(x, "stepsize__")`
#'   reaches its `pars`. **`rhat()` and `neff_ratio()` read nothing
#'   from it**: they take `pars` and `regex` and no more, so an
#'   argument they do not have, such as `variable = "x"`, is accepted
#'   and IGNORED rather than refused, and the whole set of variables
#'   comes back. brms errors on that call. A refusal is coming from
#'   `frm_check_dots()` (plan item 2.5e); until it lands, this is the
#'   accurate statement of what happens.
#' @return A ggplot object, or the diagnostic data frame / vector
#'   bayesplot returns.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     requireNamespace("bayesplot", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#'   ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
#'                    data = dd, chains = 1, iter = 500, refresh = 0)
#'   mcmc_plot(ds)
#'   mcmc_plot(ds, type = "trace", variable = "b_x")
#'   head(rhat(ds))
#' }
#' }
#' @name draws-diagnostics
NULL

#' @rdname draws-diagnostics
#' @export
mcmc_plot <- function(object, ...) UseMethod("mcmc_plot")

#' @rdname draws-diagnostics
#' @exportS3Method brms::mcmc_plot
#' @export
mcmc_plot.frmtmb_draws <- function(object, pars = NA,
                                   type = "intervals", variable = NULL,
                                   regex = FALSE, fixed = FALSE, ...) {
  sel <- draws_select_variables(object, pars, variable, regex, fixed,
                                "mcmc_plot()")
  fun <- draws_bayesplot_fun(paste0("mcmc_", type), "mcmc_plot(type =)")
  a <- draws_raw_array(object)
  keep <- sel %||% draws_outer_cols(object)
  if (!length(keep)) {
    frm_stop("mcmc_plot(): no variable was selected. variables(x) lists ",
             "what is there", call. = FALSE)
  }
  fun(a[, , keep, drop = FALSE], ...)
}

#' @rdname draws-diagnostics
#' @export
pairs.frmtmb_draws <- function(x, pars = NA, variable = NULL,
                               regex = FALSE, fixed = FALSE, ...) {
  # brms's slots: pars, variable, regex, fixed. With none given, brms
  # plots its default plot variables; here, the model's outer
  # parameters, at most four of them, because a pairs grid over every
  # parameter of a model is unreadable
  sel <- draws_select_variables(x, pars, variable, regex, fixed, "pairs()")
  a <- draws_raw_array(x)
  keep <- sel %||% utils::head(draws_outer_cols(x), 4L)
  if (!length(keep)) {
    frm_stop("pairs(): no variable was selected. variables(x) lists what ",
             "is there", call. = FALSE)
  }
  draws_bayesplot_fun("mcmc_pairs", "pairs()")(a[, , keep, drop = FALSE],
                                               ...)
}

#' One bayesplot function by name, or an error that names the argument
#' that produced the name (bayesplot's own "object not found" would name
#' neither the type nor the function that asked for it).
#'
#' @noRd
draws_bayesplot_fun <- function(nm, what) {
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    frm_stop(what, " needs the 'bayesplot' package; as.array() gives the ",
             "draws in the layout bayesplot expects if you would rather ",
             "call it yourself", call. = FALSE)
  }
  ns <- asNamespace("bayesplot")
  if (!exists(nm, envir = ns, inherits = FALSE)) {
    frm_stop(what, " asks for bayesplot::", nm, "(), which does not exist. ",
             "The argument is the function name without its 'mcmc_' ",
             "prefix, for example \"intervals\", \"trace\", \"areas\" or ",
             "\"hist\"", call. = FALSE)
  }
  get(nm, envir = ns)
}

#' @rdname draws-diagnostics
#' @export
nuts_params <- function(object, ...) UseMethod("nuts_params")

#' @rdname draws-diagnostics
#' @exportS3Method bayesplot::nuts_params
#' @export
nuts_params.frmtmb_draws <- function(object, ...) {
  draws_bayesplot_ns("nuts_params()")$nuts_params(object$stanfit, ...)
}

#' @rdname draws-diagnostics
#' @export
log_posterior <- function(object, ...) UseMethod("log_posterior")

#' @rdname draws-diagnostics
#' @exportS3Method bayesplot::log_posterior
#' @export
log_posterior.frmtmb_draws <- function(object, ...) {
  draws_bayesplot_ns("log_posterior()")$log_posterior(object$stanfit, ...)
}

#' @rdname draws-diagnostics
#' @export
rhat <- function(x, ...) UseMethod("rhat")

#' @rdname draws-diagnostics
#' @exportS3Method bayesplot::rhat
#' @rawNamespace S3method(posterior::rhat,frmtmb_draws)
#' @export
rhat.frmtmb_draws <- function(x, pars = NULL, regex = FALSE, ...) {
  a <- draws_diag_array(x, pars, regex, "rhat()")
  tab <- posterior::summarise_draws(a, rhat = posterior::rhat)
  stats::setNames(tab$rhat, tab$variable)
}

#' @rdname draws-diagnostics
#' @export
neff_ratio <- function(object, ...) UseMethod("neff_ratio")

#' @rdname draws-diagnostics
#' @exportS3Method bayesplot::neff_ratio
#' @export
neff_ratio.frmtmb_draws <- function(object, pars = NULL, regex = FALSE,
                                    ...) {
  a <- draws_diag_array(object, pars, regex, "neff_ratio()")
  tab <- posterior::summarise_draws(a, ess_bulk = posterior::ess_bulk,
                                    ess_tail = posterior::ess_tail)
  # brms's neff_ratio.brmsfit: the SMALLER of the bulk and tail
  # effective sizes over the draw count, not rstan's n_eff
  stats::setNames(pmin(tab$ess_bulk, tab$ess_tail) /
                    posterior::ndraws(a), tab$variable)
}

#' The chain-separated draws array the two convergence diagnostics are
#' computed on, restricted to the variables `pars` names.
#'
#' It goes through `as_draws_array()` rather than through the
#' `stanfit`, which is what puts frmtmb's own parameter names on the
#' result: the `stanfit` carries Stan's `beta[1]`, `betad` and `theta`,
#' which `variables(x)` does not list and `rhat(x)["x"]` cannot reach.
#'
#' `pars` is handled the way brms's `rhat.brmsfit()` handles it, which
#' is NOT the way brms's `as.mcmc()`, `mcmc_plot()` and
#' `posterior_interval()` handle an argument of the same name. Those
#' three go through brms's `extract_pars()`, where `NA` means "every
#' variable" and a string is a regular expression. These two instead
#' pass `variable = pars` to `as_draws_array()`, so `NULL` means every
#' variable and a string is an EXACT name unless `regex = TRUE`. The
#' difference is brms's, not this package's, and it is measured in
#' `dev/brmsmatch-findings.md`.
#'
#' @noRd
draws_diag_array <- function(x, pars, regex, what) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    frm_stop(what, " needs the 'posterior' package: it computes the ",
             "diagnostic brms reports, which is posterior's. The sampler's ",
             "own classic split-R-hat and n_eff are in `x$stanfit` without ",
             "it", call. = FALSE)
  }
  a <- draws_as_array(x)
  if (is.null(pars)) return(a)
  check_flag(regex, "regex")
  draws_subset_variable(a, pars, regex, what)
}

#' bayesplot's namespace, or an error naming the accessor that wanted it.
#'
#' @noRd
draws_bayesplot_ns <- function(what) {
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    frm_stop(what, " needs the 'bayesplot' package: it reads the sampler ",
             "diagnostics off the stanfit, which is `ds$stanfit` if you ",
             "would rather use rstan directly", call. = FALSE)
  }
  asNamespace("bayesplot")
}

# ---- mixture membership ----------------------------------------------

#' Posterior mixture-component probabilities
#'
#' For a [frmtmb::mixture()], [frmtmb::mixture_mvn()] or `frmtmb.latent::lca()` fit, the posterior
#' probability that each observation came from each component,
#' propagating the uncertainty in the parameters: the fit-side
#' [frmtmb::mixture_probs()] computation is run at every draw. brms calls this
#' `pp_mixture()`.
#'
#' The argument order is brms's, so `summary` sits in brms's own
#' eighth position and not in the second: the second is `newdata`.
#'
#' @param x A `frmtmb_draws` from [frm_sample()].
#' @param newdata,re_formula Accepted for brms's signature and refused:
#'   [frmtmb::mixture_probs()] is a statement about the rows the model
#'   was fitted on.
#' @param resp The response whose mixture to report, for a
#'   multivariate model.
#' @param log If `TRUE`, log probabilities.
#' @param summary If `TRUE` (the default, as in brms), an
#'   `observations x statistics x components` array of summaries;
#'   otherwise the raw `draws x observations x components` array.
#' @param robust If `TRUE`, median and MAD instead of mean and SD.
#' @param probs The two quantiles the summary reports.
#' @param ndraws Number of draws to use (default: all).
#' @param draw_ids The draws to use, by row index, instead of the
#'   evenly spaced subsample `ndraws` takes.
#' @param ... Unused.
#' @return An array; see `summary`. For a group-level mixture
#'   (`mixture(groups = )`, `frmtmb.latent::lca()`) the rows are groups, as in
#'   [frmtmb::mixture_probs()].
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(4)
#'   dd <- data.frame(y = c(rnorm(60, -2), rnorm(60, 3)))
#'   fit <- frm(bf(y ~ 1), family = frmtmb::mixture(gaussian(), gaussian()),
#'              data = dd)
#'   ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0)
#'   head(pp_mixture(ds)[, "Estimate", ])
#' }
#' }
#' @export
pp_mixture <- function(x, ...) UseMethod("pp_mixture")

#' @rdname pp_mixture
#' @exportS3Method brms::pp_mixture
#' @export
pp_mixture.frmtmb_draws <- function(x, newdata = NULL,
                                    re_formula = arg_unset(),
                                    resp = NULL, ndraws = NULL,
                                    draw_ids = NULL,
                                    log = FALSE, summary = TRUE,
                                    robust = FALSE,
                                    probs = c(0.025, 0.975), ...) {
  draws_refuse_newdata(newdata, re_formula, arg_unset(), "pp_mixture()",
                       "mixture_probs() classifies the rows the model ",
                       "was fitted on: the component probability of a ",
                       "row is a statement about that row's own ",
                       "response, and newdata carries no response to ",
                       "classify")
  if (!is.null(resp) && !resp %in% names(x$fit$spec$responses)) {
    frm_stop("pp_mixture(resp = \"", resp, "\") names no response of this ",
             "model; it has ",
             paste(names(x$fit$spec$responses), collapse = ", "),
             call. = FALSE)
  }
  check_flag(log, "log")
  check_flag(summary, "summary")
  check_flag(robust, "robust")
  idx <- draws_par_index(x$fit)
  rows <- draws_subsample(x, ndraws, draw_ids)
  out <- NULL
  for (k in seq_along(rows)) {
    P <- mixture_probs(draws_fit_at(x, rows[k], idx))
    if (is.null(out)) {
      out <- array(NA_real_, c(length(rows), nrow(P), ncol(P)),
                   dimnames = list(NULL, rownames(P), colnames(P)))
    }
    out[k, , ] <- P
  }
  if (log) out <- log(out)
  if (!summary) return(out)
  draws_summarize_margin(out, probs, robust)
}

#' brms's `posterior_summary()` over the first margin of a
#' draws x rows x components array, with its statistic and quantile
#' switches.
#'
#' @noRd
draws_summarize_margin <- function(out, probs = c(0.025, 0.975),
                                   robust = FALSE) {
  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs) ||
        any(probs <= 0) || any(probs >= 1)) {
    frm_stop("probs must be two numbers strictly between 0 and 1",
             call. = FALSE)
  }
  qn <- paste0("Q", probs * 100)
  st <- array(NA_real_, c(dim(out)[2L], 4L, dim(out)[3L]),
              dimnames = list(dimnames(out)[[2L]],
                              c("Estimate", "Est.Error", qn),
                              dimnames(out)[[3L]]))
  st[, "Estimate", ] <- apply(out, c(2, 3),
                              if (robust) stats::median else mean)
  st[, "Est.Error", ] <- apply(out, c(2, 3),
                               if (robust) stats::mad else stats::sd)
  st[, qn[1L], ] <- apply(out, c(2, 3), stats::quantile, probs[1L])
  st[, qn[2L], ] <- apply(out, c(2, 3), stats::quantile, probs[2L])
  st
}

#' The refusal the methods that carry brms's `newdata` and
#' `re_formula` slots without supporting them share.
#'
#' They carry the slots so that a positional brms call asks the same
#' question it asks in brms; they refuse the values so that it does not
#' get a different question's answer instead.
#'
#' @noRd
draws_refuse_newdata <- function(newdata, re_formula, re.form, what,
                                 ...) {
  if (!is.null(newdata)) {
    frm_stop(what, " does not take newdata. ", ..., call. = FALSE)
  }
  # core's marker class, read directly: `is_arg_unset()` is internal to
  # frmtmb while `arg_unset()` is the exported half of the pair.
  #
  # `re_formula = NULL` is NOT refused. It is brms's own default and it
  # means "condition on the group-level values", which is what these
  # methods already do, so refusing it would refuse the no-op a
  # positional brms call spells out. Only a value that asks for
  # something else is refused.
  asks <- function(v) {
    !inherits(v, "frmtmb_arg_unset") && !is.null(v)
  }
  if (asks(re_formula) || asks(re.form)) {
    frm_stop(what, " does not take re_formula. ", ..., call. = FALSE)
  }
  invisible(NULL)
}

# ---- refusals and renamed spellings ----------------------------------

#' Methods a ported brms script may call that frmtmb does not have
#'
#' These `brmsfit` methods either describe machinery frmtmb does not use
#' (Stan code and Stan data) or are brms spellings that have been
#' renamed. They are defined so that a ported script gets the reason and
#' the replacement rather than "could not find function", which is what
#' the vignette-port audit measured most of its post-processing failures
#' as.
#'
#' @param object,x,... Ignored; these functions always stop.
#' @return These functions never return; they signal an error.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(1)
#'   dd <- data.frame(x = rnorm(40))
#'   dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
#'   fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
#'   ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0)
#'   # each refusal names its reason and the replacement
#'   try(stancode(ds))
#'   try(standata(ds))
#' }
#' }
#' @name frmtmb-draws-refusals
NULL

#' @rdname frmtmb-draws-refusals
#' @export
stancode <- function(object, ...) UseMethod("stancode")

#' @rdname frmtmb-draws-refusals
#' @exportS3Method brms::stancode
#' @export
stancode.frmtmb_draws <- function(object, ...) {
  frm_stop("stancode() has no meaning for frmtmb: there is no Stan ",
           "program. The model is an R closure built by build_objective() ",
           "from the assembled frame and differentiated by RTMB, and the ",
           "closure IS the source: print `ds$fit$obj$fn` for the ",
           "evaluator and `ds$fit$frame` for everything baked into it",
           call. = FALSE)
}

#' @rdname frmtmb-draws-refusals
#' @export
standata <- function(object, ...) UseMethod("standata")

#' @rdname frmtmb-draws-refusals
#' @exportS3Method brms::standata
#' @export
standata.frmtmb_draws <- function(object, ...) {
  frm_stop("standata() has no meaning for frmtmb: nothing is exported to a ",
           "Stan data list. The assembled frame `ds$fit$frame` holds the ",
           "same content (the response, the design matrices, the sparse Z, ",
           "the addition terms), and model.matrix(), getME() and ",
           "model.frame() read the pieces of it individually",
           call. = FALSE)
}

#' @rdname frmtmb-draws-refusals
#' @exportS3Method brms::expose_functions
#' @export
expose_functions.frmtmb_draws <- function(x, ...) {
  frm_stop("expose_functions() has nothing to expose: brms compiles Stan ",
           "functions and this makes them callable from R, while a frmtmb ",
           "custom family is already plain R: the lpdf you passed to ",
           "custom_family() is an R function you can call directly",
           call. = FALSE)
}


#' @rdname frmtmb-draws-refusals
#' @export
plot.frmtmb_draws <- function(x, ...) {
  frm_stop("plot() has no display for frmtmb draws: brms's default panel ",
           "is the trace-and-density view, which mcmc_plot(x) renders ",
           "here (mcmc_plot(x, type = \"trace\") for the traces alone)",
           call. = FALSE)
}

#' @rdname frmtmb-draws-refusals
#' @export
update.frmtmb_draws <- function(object, ...) {
  frm_stop("update() has no method for draws: the sampled object carries ",
           "no formula to revise. update() the underlying frmtmb fit and ",
           "frm_sample() the result", call. = FALSE)
}

#' @rdname frmtmb-draws-refusals
#' @export
restructure <- function(x, ...) UseMethod("restructure")

#' @rdname frmtmb-draws-refusals
#' @exportS3Method brms::restructure
#' @export
restructure.frmtmb_draws <- function(x, ...) {
  frm_stop("restructure() is brms's upgrade path for objects saved by an ",
           "older brms; frmtmb has no such conversion. A draws object from ",
           "an older frmtmb is re-created by re-running frm_sample()",
           call. = FALSE)
}

# brms's deprecated draws accessors (item 2.6f).
#
# These three were refused as "the deprecated brms spelling". brms keeps
# all three LIVE, with a deprecation warning, and the user's rule of
# 2026-09-15 is that frmtmb.sample matches brms. A refusal is not a
# deprecation warning: a ported script stops at it. So each one now
# does what brms's does and warns the way brms's warns.

#' Deprecated brms draws accessors
#'
#' @description
#' `posterior_samples()`, `nsamples()` and `parnames()` are brms
#' spellings that brms itself deprecated and still answers. Each warns
#' and then does what brms's does, so a ported script runs.
#'
#' * `posterior_samples(x)` is `as.data.frame(x)`; `pars` is a regular
#'   expression unless `fixed = TRUE`, as in brms.
#' * `nsamples(x)` is [ndraws()].
#' * `parnames(x)` is [variables()].
#'
#' Use `as_draws(x)` for a posterior draws object, `as.data.frame(x)`,
#' [ndraws()] and [variables()] in new code.
#'
#' @param x,object A `frmtmb_draws`.
#' @param pars Variable names. A regular expression unless
#'   `fixed = TRUE`; `NA` (the default) takes all of them.
#' @param fixed If `TRUE`, `pars` names variables exactly.
#' @param add_chain If `TRUE`, add the `chain` and `iter` columns brms
#'   adds.
#' @param subset Draw indices to keep.
#' @param as.matrix,as.array Return a matrix or a
#'   draws-by-chains-by-variables array instead of a data frame.
#' @param incl_warmup Refused: `frm_sample()` discards the warmup, so
#'   there is no warmup draw to count.
#' @param ... Refused: an argument the method does not have is an error
#'   naming it.
#' @return A data frame (or matrix, or array) of draws for
#'   `posterior_samples()`, one integer for `nsamples()`, and a
#'   character vector for `parnames()`.
#' @name frmtmb-draws-deprecated
NULL

#' @rdname frmtmb-draws-deprecated
#' @export
posterior_samples <- function(x, pars = NA, ...) {
  UseMethod("posterior_samples")
}

#' @rdname frmtmb-draws-deprecated
#' @exportS3Method brms::posterior_samples
#' @rawNamespace S3method(gratia::posterior_samples,frmtmb_draws)
#' @export
posterior_samples.frmtmb_draws <- function(x, pars = NA, fixed = FALSE,
                                           add_chain = FALSE, subset = NULL,
                                           as.matrix = FALSE,
                                           as.array = FALSE, ...) {
  frm_check_dots(...)
  check_flag(fixed, "fixed")
  check_flag(add_chain, "add_chain")
  check_flag(as.matrix, "as.matrix")
  check_flag(as.array, "as.array")
  frm_warning("Method 'posterior_samples' is deprecated. Please see ",
              "?as_draws for recommended alternatives.", call. = FALSE)
  if (as.matrix && as.array) {
    frm_stop("posterior_samples(): 'as.matrix' and 'as.array' cannot both ",
             "be TRUE", call. = FALSE)
  }
  # `pars` here is a regular expression by default, which is brms's rule
  # and NOT the `variable` argument's: as.data.frame(pars = ) would warn
  # about a deprecated alias of its own and match exactly
  variable <- if (anyNA(pars)) NULL else pars
  out <- if (as.matrix) {
    as.matrix(x, variable = variable, regex = !fixed, draw = subset)
  } else if (as.array) {
    as.array(x, variable = variable, regex = !fixed, draw = subset)
  } else {
    as.data.frame(x, variable = variable, regex = !fixed, draw = subset)
  }
  if (add_chain && !as.array) {
    nc <- nchains(x)
    ni <- nrow(as.matrix(out)) %/% max(nc, 1L)
    out <- as.data.frame(out)
    out$chain <- factor(rep(seq_len(nc), each = ni))
    out$iter <- rep(seq_len(ni), nc)
  }
  out
}

#' @rdname frmtmb-draws-deprecated
#' @export
nsamples <- function(object, ...) UseMethod("nsamples")

#' @rdname frmtmb-draws-deprecated
#' @exportS3Method rstantools::nsamples
#' @export
nsamples.frmtmb_draws <- function(object, subset = NULL,
                                  incl_warmup = FALSE, ...) {
  frm_check_dots(...)
  check_flag(incl_warmup, "incl_warmup")
  frm_warning("'nsamples.frmtmb_draws' is deprecated. Please use 'ndraws' ",
              "instead.", call. = FALSE)
  if (incl_warmup) {
    frm_stop("nsamples(incl_warmup = TRUE) has nothing to count: ",
             "frm_sample() discards the warmup rather than storing it, so ",
             "the object carries post-warmup draws only. ndraws(x) is that ",
             "count", call. = FALSE)
  }
  if (!is.null(subset)) return(length(subset))
  ndraws(object)
}

#' @rdname frmtmb-draws-deprecated
#' @export
parnames <- function(x, ...) UseMethod("parnames")

#' @rdname frmtmb-draws-deprecated
#' @exportS3Method brms::parnames
#' @export
parnames.frmtmb_draws <- function(x, ...) {
  frm_check_dots(...)
  frm_warning("'parnames' is deprecated. Please use 'variables' instead.",
              call. = FALSE)
  variables(x)
}


# ---- brms's summarizing post-processing methods ----------------------
#
# brms's fitted(), predict() and residuals() on a fit are summaries of
# posterior_epred(), posterior_predict() and predictive_error(). Without
# them `fitted(ds)` reached stats::fitted.default(), which reads
# `object$fitted.values` and returns NULL: a wrong answer with nothing
# said (defect S4 of dev/brmsport-findings.md).

#' Summaries of the posterior predictive quantities
#'
#' @description
#' brms's three summarizing methods, each a summary of the draws method
#' beside it:
#'
#' * `fitted()` summarizes [posterior_epred()] (or
#'   [posterior_linpred()] under `scale = "linear"`);
#' * `predict()` summarizes [posterior_predict()];
#' * `residuals()` summarizes [predictive_error()].
#'
#' `summary = FALSE` returns the draws themselves, which is what the
#' method it wraps returns.
#'
#' @param object A `frmtmb_draws` from `frm_sample()`.
#' @param newdata,re_formula,resp,dpar,nlpar,ndraws,draw_ids Passed to
#'   the draws method, where they are documented.
#' @param scale `"response"` for the expected response, `"linear"` for
#'   the linear predictor.
#' @param transform Applied to the draws before they are summarized.
#' @param negative_rt,ntrys,cores,sort brms's remaining arguments,
#'   carried so that a positional brms call asks the same question; each
#'   is passed on or refused by the method it belongs to.
#' @param type For `residuals()`, `"ordinary"` (brms's name for the raw
#'   error) or `"pearson"`.
#' @param method For `residuals()`, which predictive distribution the
#'   error is taken against.
#' @param summary If `FALSE`, the draws instead of their summary.
#' @param robust If `TRUE`, the median and MAD instead of the mean and
#'   standard deviation.
#' @param probs Probabilities of the quantile columns.
#' @param ... Passed to the draws method, which is where brms's
#'   `point_estimate` and `ndraws_point_estimate` are answered.
#' @return With `summary = TRUE` an observations-by-four matrix, or an
#'   observations-by-four-by-K array for a category distribution. With
#'   `summary = FALSE` the draws.
#' @name frmtmb-draws-summaries
NULL

#' @rdname frmtmb-draws-summaries
#' @export
fitted.frmtmb_draws <- function(object, newdata = NULL,
                                re_formula = arg_unset(),
                                scale = c("response", "linear"),
                                resp = NULL, dpar = NULL, nlpar = NULL,
                                ndraws = NULL, draw_ids = NULL,
                                sort = FALSE, summary = TRUE,
                                robust = FALSE, probs = c(0.025, 0.975),
                                ...) {
  scale <- frm_match_arg(scale)
  draws_refuse_sort(sort, "fitted()")
  # checked here as well as in posterior_epred(), so a refusal names the
  # function the caller called
  draws_refuse_new_levels(object, newdata, list(...), "fitted()")
  out <- if (identical(scale, "response")) {
    posterior_epred(object, newdata = newdata, re_formula = re_formula,
                    resp = resp, dpar = dpar, nlpar = nlpar,
                    ndraws = ndraws, draw_ids = draw_ids, ...)
  } else {
    posterior_linpred(object, newdata = newdata, re_formula = re_formula,
                      resp = resp, dpar = dpar, nlpar = nlpar,
                      ndraws = ndraws, draw_ids = draw_ids, ...)
  }
  draws_summarize_or_not(out, summary, probs, robust)
}

#' @rdname frmtmb-draws-summaries
#' @export
predict.frmtmb_draws <- function(object, newdata = NULL,
                                 re_formula = arg_unset(),
                                 transform = NULL, resp = NULL,
                                 negative_rt = FALSE, ndraws = NULL,
                                 draw_ids = NULL, sort = FALSE,
                                 ntrys = NULL, cores = NULL,
                                 summary = TRUE, robust = FALSE,
                                 probs = c(0.025, 0.975), ...) {
  draws_refuse_sort(sort, "predict()")
  draws_refuse_ntrys_cores(ntrys, cores, "predict()")
  draws_refuse_new_levels(object, newdata, list(...), "predict()")
  out <- posterior_predict(object, newdata = newdata,
                           re_formula = re_formula, transform = transform,
                           resp = resp, negative_rt = negative_rt,
                           ndraws = ndraws, draw_ids = draw_ids, ...)
  if (!summary) return(out)
  # a category-valued response has no mean to summarize, so brms reports
  # the simulated proportion of each category, the same shape the fit
  # method returns
  fit <- draws_base_fit(object)
  rspec <- fit$spec$responses[[resp %||% names(fit$spec$responses)[1L]]]
  if (fam_is_category_valued(rspec$family)) {
    return(predict_category_props(fit, rspec, out))
  }
  draws_summarize_or_not(out, summary, probs, robust)
}

#' @rdname frmtmb-draws-summaries
#' @export
residuals.frmtmb_draws <- function(object, newdata = NULL,
                                   re_formula = arg_unset(),
                                   method = "posterior_predict",
                                   type = c("ordinary", "pearson"),
                                   resp = NULL, ndraws = NULL,
                                   draw_ids = NULL, sort = FALSE,
                                   summary = TRUE, robust = FALSE,
                                   probs = c(0.025, 0.975), ...) {
  type <- frm_match_arg(type)
  draws_refuse_sort(sort, "residuals()")
  draws_refuse_new_levels(object, newdata, list(...), "residuals()")
  out <- predictive_error(object, newdata = newdata,
                          re_formula = re_formula, method = method,
                          resp = resp, ndraws = ndraws,
                          draw_ids = draw_ids, ...)
  if (identical(type, "pearson")) {
    # brms divides the error draws by the predictive standard deviation
    # of the same draws, which is what makes it a Pearson residual
    pp <- posterior_predict(object, newdata = newdata,
                            re_formula = re_formula, resp = resp,
                            ndraws = ndraws, draw_ids = draw_ids, ...)
    sdv <- apply(pp, 2L, stats::sd)
    out <- sweep(out, 2L, sdv, "/")
  }
  draws_summarize_or_not(out, summary, probs, robust)
}

#' brms's `ntrys` and `cores`, carried so that a positional brms call
#' lands where brms lands it and refused with the reason.
#'
#' @noRd
draws_refuse_ntrys_cores <- function(ntrys, cores, what) {
  if (!is.null(ntrys)) {
    frm_stop(what, " cannot honor `ntrys` yet: the rejection limit of a ",
             "trunc()ed draw is the family simulator's own", call. = FALSE)
  }
  if (!is.null(cores)) {
    frm_stop(what, " cannot honor `cores`: the draws are replayed in one ",
             "process. Lower `ndraws` if a call is too slow", call. = FALSE)
  }
  invisible(NULL)
}

#' @noRd
draws_refuse_sort <- function(sort, what) {
  check_flag(sort, "sort")
  if (sort) {
    frm_stop(what, " cannot honor sort = TRUE: rows come back in the order ",
             "of the data, always", call. = FALSE)
  }
  invisible(NULL)
}

#' @noRd
draws_summarize_or_not <- function(out, summary, probs, robust) {
  check_flag(summary, "summary")
  check_flag(robust, "robust")
  if (!summary) return(out)
  posterior_summary(out, probs = probs, robust = robust)
}

# The spellings draws_refuse_new_levels() reads from the dots, so a
# method that refuses unknown dots names lets these through to it.
draws_new_level_args <- c("allow_new_levels", "allow.new.levels",
                          "sample_new_levels")

#' Refuse an unseen grouping level on draws, with a message that is true,
#' and ONLY when the call asks for one.
#'
#' brms answers `allow_new_levels = TRUE` on draws by drawing a new
#' level's effect per posterior draw (`sample_new_levels = "gaussian"`)
#' or by resampling the draws of the levels it saw (the default,
#' `"uncertainty"`); neither is built for `frmtmb_draws`. Before this the
#' argument vanished into `...` and the design builder refused an unseen
#' level with "Use allow_new_levels = TRUE", which the caller had just
#' done.
#'
#' A refusal on the argument's mere PRESENCE is too wide:
#' `allow_new_levels = FALSE` is brms's default, and `TRUE` with no
#' `newdata`, or with only levels the fit saw, changes nothing in brms
#' either. Those answer, as they did before. What is refused is `TRUE`
#' together with a `newdata` that holds a level the fit did not see,
#' which is decided by asking the fit: the rows build without the flag,
#' or build only with it.
#'
#' @noRd
draws_refuse_new_levels <- function(object, newdata, dots, fn) {
  anl <- dots[["allow_new_levels"]] %||% dots[["allow.new.levels"]]
  if (!isTRUE(anl) || is.null(newdata)) return(invisible(NULL))
  if (!draws_has_unseen_level(draws_base_fit(object), newdata)) {
    return(invisible(NULL))
  }
  frm_stop(fn, " on draws cannot predict a grouping level the fit did not ",
           "see: brms draws that level's effect from each posterior draw ",
           "or resamples the levels it saw, and neither is implemented ",
           "for frmtmb_draws yet. predict(fit, allow_new_levels = TRUE) on ",
           "the maximum-likelihood fit draws it, and re_formula = NA ",
           "predicts at the population level", call. = FALSE)
}

#' Whether `newdata` holds a grouping level the fit did not see: its
#' rows fail to build as they are and build once unseen levels are
#' allowed. A failure both ways is some other fault, left to the call
#' itself to report.
#'
#' @noRd
draws_has_unseen_level <- function(fit, newdata) {
  known <- tryCatch({
    frm_linpred(fit, newdata = newdata)
    TRUE
  }, error = function(e) FALSE)
  if (known) return(FALSE)
  tryCatch({
    frm_linpred(fit, newdata = newdata, allow_new_levels = TRUE)
    TRUE
  }, error = function(e) FALSE)
}
