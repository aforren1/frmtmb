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

#' The originating fit with its estimates replaced by one draw. The draw
#' IS a parameter vector, so the object is a legitimate fit from here on
#' even when the draws came from a formula with no ML mode behind it.
#'
#' @noRd
draws_fit_at <- function(x, i, idx = draws_par_index(x$fit)) {
  fit <- draws_base_fit(x)
  est <- fit$frame[["par_template"]]   # mapped betad entries keep link(const)
  row <- x$draws[i, ]
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
      stop("`ndraws` and `draw_ids` both choose which draws to use, so ",
           "only one of them can be given: `ndraws` takes an evenly ",
           "spaced subsample of that size, `draw_ids` takes the rows ",
           "you name", call. = FALSE)
    }
    ok <- is.numeric(draw_ids) && length(draw_ids) &&
      !anyNA(draw_ids) && all(draw_ids == round(draw_ids)) &&
      all(draw_ids >= 1L) && all(draw_ids <= n)
    if (!ok) {
      stop("`draw_ids` must be whole numbers between 1 and ", n,
           ", the number of draws", call. = FALSE)
    }
    return(as.integer(draw_ids))
  }
  if (is.null(ndraws) || ndraws >= n) return(seq_len(n))
  round(seq(1, n, length.out = ndraws))
}

#' @export
summary.frmtmb_draws <- function(object, ...) {
  m <- object$draws
  keep <- setdiff(colnames(m),
                  c("lp__", grep("^b\\[", colnames(m), value = TRUE)))
  tab <- t(vapply(keep, function(nm) {
    c(mean = mean(m[, nm]), sd = stats::sd(m[, nm]),
      `2.5%` = unname(stats::quantile(m[, nm], 0.025)),
      `97.5%` = unname(stats::quantile(m[, nm], 0.975)))
  }, numeric(4)))
  # The SAME three quantities `brms:::summary.brmsfit` reports, from the
  # same posterior functions: `Rhat = posterior::rhat`, `Bulk_ESS =
  # posterior::ess_bulk`, `Tail_ESS = posterior::ess_tail`. They used to
  # be rstan's `n_eff` and classic split-Rhat, which put a column headed
  # `Rhat` beside `rhat(object)` carrying a different definition, below
  # 1 in one and above it in the other on a measured fit.
  if (requireNamespace("posterior", quietly = TRUE)) {
    a <- posterior::subset_draws(as_draws_array(object), variable = keep)
    d <- posterior::summarise_draws(a, Rhat = posterior::rhat,
                                    Bulk_ESS = posterior::ess_bulk,
                                    Tail_ESS = posterior::ess_tail)
    i <- match(keep, d$variable)
    tab <- cbind(tab, Rhat = d$Rhat[i], Bulk_ESS = d$Bulk_ESS[i],
                 Tail_ESS = d$Tail_ESS[i])
  }
  tab
}

#' @exportS3Method nlme::fixef
#' @export
fixef.frmtmb_draws <- function(object, ...) {
  # the draws-side spelling: parenthesis-free, matching the draws
  # matrix, summary(), variables() and hypothesis()
  nm <- par_name_bare(estimated_coef_names(object$fit))
  idx <- draws_par_index(object$fit)
  cols <- c(idx$beta, idx$betad)
  m <- object$draws[, cols, drop = FALSE]
  out <- cbind(
    Estimate = colMeans(m),
    Est.Error = apply(m, 2, stats::sd),
    Q2.5 = apply(m, 2, stats::quantile, 0.025),
    Q97.5 = apply(m, 2, stats::quantile, 0.975)
  )
  rownames(out) <- nm
  out
}

#' @exportS3Method nlme::VarCorr
#' @export
VarCorr.frmtmb_draws <- function(x, sigma = 1, ...) {
  # a structural template only: every number in `base` is replaced by a
  # posterior summary below, so the starting values of a formula-sampled
  # object never reach the result
  fit <- draws_base_fit(x)
  idx <- draws_par_index(fit)
  if (is.null(idx$theta)) return(NULL)
  th_draws <- x$draws[, idx$theta, drop = FALSE]
  per_draw <- lapply(seq_len(nrow(th_draws)), function(i) {
    fit$estimates[["theta"]] <- th_draws[i, ]
    # the exported generic, not core's method by name: the method is
    # core's to move or rename and this only ever wanted the dispatch
    as.data.frame(VarCorr(fit))$sdcor
  })
  base <- as.data.frame(VarCorr(fit))
  M <- do.call(rbind, per_draw)
  base$estimate <- colMeans(M)
  base$lwr <- apply(M, 2, stats::quantile, 0.025)
  base$upr <- apply(M, 2, stats::quantile, 0.975)
  base$vcov <- NULL
  base$sdcor <- NULL
  base
}

#' @exportS3Method rstantools::prior_summary
#' @export
prior_summary.frmtmb_draws <- function(object, ...) {
  pl <- object$fit$prior
  if (is.null(pl) || (!length(unclass(pl)) &&
                        !length(attr(pl, "overrides")))) {
    cat("No priors were used (flat improper priors on the outer ",
        "parameters).\n", sep = "")
    return(invisible(NULL))
  }
  pl
}

#' @exportS3Method nlme::ranef
#' @export
ranef.frmtmb_draws <- function(object, ...) {
  fit <- object$fit
  if (!length(fit$frame[["re_blocks"]])) return(list())
  idx <- draws_par_index(fit)
  n <- nrow(object$draws)
  per <- vector("list", n)
  for (i in seq_len(n)) {
    per[[i]] <- ranef(draws_fit_at(object, i, idx))
  }
  # brms shape: per term, a levels x statistics x coefficients array.
  # BY POSITION, not by name: core's ranef() keys its list by the
  # GROUPING FACTOR now (brms's and lme4's key), so two blocks on one
  # factor share a name and assigning by name wrote one entry twice,
  # silently dropping the second block. core's own print() and
  # as.data.frame() index by position for the same reason.
  out <- vector("list", length(per[[1]]))
  names(out) <- names(per[[1]])
  for (bi in seq_along(per[[1]])) {
    M0 <- per[[1]][[bi]]
    A <- vapply(per, function(r) r[[bi]], M0)
    st <- array(NA_real_, c(nrow(M0), 4L, ncol(M0)),
                dimnames = list(rownames(M0),
                                c("Estimate", "Est.Error",
                                  "Q2.5", "Q97.5"),
                                colnames(M0)))
    st[, "Estimate", ] <- apply(A, c(1, 2), mean)
    st[, "Est.Error", ] <- apply(A, c(1, 2), stats::sd)
    st[, "Q2.5", ] <- apply(A, c(1, 2), stats::quantile, 0.025)
    st[, "Q97.5", ] <- apply(A, c(1, 2), stats::quantile, 0.975)
    # the block label rides along exactly as core's ranef() carries it
    attr(st, "term") <- attr(M0, "term")
    out[[bi]] <- st
  }
  out
}

#' @exportS3Method brms::hypothesis
#' @export
hypothesis.frmtmb_draws <- function(x, hypothesis, alpha = 0.05,
                                    class = NULL, group = NULL, ...) {
  # Arm core's reserved-name shadowing note for this call. It used to
  # be armed by core's GENERIC, and this method relied on that. The
  # exported `hypothesis` now resolves to brms's generic whenever brms
  # is loaded, so core moved the arming into its methods, and this one
  # lost the note in every session until it armed it too. Measured on
  # a covariate named `sigma`: 1 note on the fit, 0 on its draws.
  old <- hyp_shadow_arm()
  on.exit(hyp_shadow_disarm(old), add = TRUE)
  fit <- x$fit
  vo <- hyp_vals_only(fit)
  hp <- hyp_parse_all(hypothesis,
                      names(hyp_env_vals(fit, vo$vals, vo$comp)),
                      class, group)
  exs <- hp$exprs
  idx <- draws_par_index(fit)
  n <- nrow(x$draws)
  draws <- matrix(NA_real_, n, length(exs),
                  dimnames = list(NULL, hypothesis))
  for (i in seq_len(n)) {
    sh <- draws_fit_at(x, i, idx)
    w <- hyp_vals_only(sh)
    draws[i, ] <- vapply(exs, function(ex) {
      hyp_eval(sh, ex, w$vals, w$comp)
    }, numeric(1))
  }
  # brms's Evid.Ratio and Post.Prob, in this frame's own lower-case
  # spelling. Only the draws surface carries them: they are posterior
  # quantities, and a Wald or profile hypothesis has no posterior
  # density to put over a prior one.
  #
  # A POINT null is the spelling that carries "=". The bare-quantity
  # spelling this package also accepts parses `two.sided` as well, and
  # it is a summary request rather than a test, so it must not be
  # weighed or complained about.
  is_point <- hp$dir == "two.sided" &
    grepl("=", hypothesis, fixed = TRUE)
  ev <- er_evidence(x, exs, hp$dir, draws, hypothesis, is_point)
  rows <- lapply(seq_along(exs), function(k) {
    t_k <- draws[, k]
    dir <- hp$dir[k]
    lo <- if (dir == "two.sided") alpha / 2 else alpha
    data.frame(hypothesis = hypothesis[k], estimate = mean(t_k),
               se = stats::sd(t_k),
               lwr = if (dir == "less") -Inf else
                 unname(stats::quantile(t_k, lo)),
               upr = if (dir == "greater") Inf else
                 unname(stats::quantile(t_k, 1 - lo)),
               z = mean(t_k) / stats::sd(t_k),
               p = hyp_tail_p(t_k, dir),
               evid_ratio = ev$evid_ratio[k],
               post_prob = ev$post_prob[k])
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "method") <- "posterior"
  attr(out, "alpha") <- alpha
  attr(out, "direction") <- hp$dir
  attr(out, "draws") <- draws
  attr(out, "nsim") <- n
  attr(out, "converged") <- rep(TRUE, n)
  # an attribute rather than a column: the frame stays the shape a
  # ported brms script indexes, and the Monte Carlo error of a kernel
  # density ratio is exactly what a reader needs before believing one
  attr(out, "evid_ratio_mcse") <- ev$mcse
  class(out) <- c("frmtmb_hypothesis", "data.frame")
  out
}

#' Expected-value and predictive draws from sampled parameters
#'
#' `posterior_epred()` evaluates the response-scale expectation per
#' draw; `posterior_predict()` additionally simulates responses from
#' the family, giving the posterior predictive distribution. Both
#' condition on each draw's own random effects (`re_formula = NA` drops
#' them; `re.form` is the accepted alias, see *Argument spellings*).
#'
#' @section Categorical outcomes:
#' An ordinal family predicts a DISTRIBUTION per observation, not one
#' number: each draw's `predict(type = "response")` is an `n x K`
#' matrix of category probabilities. Those stack into a 3-D
#' `draws x observations x categories` array. `dimnames` are
#' `list(NULL, <observation names or NULL>, <category levels>)`, so
#' `ep[, , "high"]` is the draws-by-observations matrix for one
#' category and `ep[k, , ]` is draw `k`'s own `n x K` prediction, the
#' matrix `predict(type = "response")` returns. Every `ep[k, i, ]`
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
#' frmtmb answers to two dialects, and this family sits on the seam.
#' The rule is that a brms-NAMED function speaks brms's argument names,
#' while frmtmb's own fit surface ([frmtmb::predict.frmtmb_fit()],
#' [frmtmb::simulate.frmtmb_fit()], [frmtmb::frm_bootstrap()]) keeps lme4's, because
#' that is the heritage each name comes from and a reader should be able
#' to tell which library a call was written against.
#'
#' `posterior_epred()` and its relatives are brms functions, so the
#' random-effect switch is `re_formula`. They also SHIPPED taking
#' lme4's `re.form`, so that spelling keeps working and means exactly
#' the same thing: both names feed one internal setting, and whichever
#' one is given wins. brms does the same on `posterior_epred.brmsfit()`,
#' which carries `re_formula` and `re.form` side by side.
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
#'   [frmtmb::predict.frmtmb_fit()]'s `re.form`; see *Argument spellings*.
#' @param re.form lme4's spelling of `re_formula`, accepted as an alias.
#'   Pass one or the other, not both.
#' @param ndraws Number of draws to use (default: all).
#' @param draw_ids The draws to use, by row index, instead of the
#'   evenly spaced subsample `ndraws` takes. Give one or the other.
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
#' @param ... Unused.
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
                                         draw_ids = NULL, ...) {
  re_form <- re_form_arg(re_formula, re.form, "posterior_epred()")
  dpar <- draws_dpar_arg(dpar, nlpar, "posterior_epred()")
  idx <- draws_par_index(object$fit)
  rows <- draws_subsample(object, ndraws, draw_ids)
  out <- NULL
  cat_out <- FALSE
  for (k in seq_along(rows)) {
    sh <- draws_fit_at(object, rows[k], idx)
    p <- predict(sh, newdata = newdata, resp = resp, dpar = dpar,
                 re.form = re_form, type = "response")
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
                                           draw_ids = NULL, ...) {
  re_form <- re_form_arg(re_formula, re.form, "posterior_linpred()")
  dpar <- draws_dpar_arg(dpar, nlpar, "posterior_linpred()")
  if (!is.null(incl_thres) && !identical(incl_thres, FALSE)) {
    stop("posterior_linpred(incl_thres = TRUE) subtracts an ordinal ",
         "family's thresholds from the linear predictor, which brms ",
         "supports for cumulative families alone. frmtmb keeps the ",
         "thresholds out of the predictor: predict(type = \"link\") ",
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
    p <- predict(sh, newdata = newdata, resp = resp, dpar = dpar,
                 re.form = re_form,
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
    stop(what, " was given both `dpar` and `nlpar`. frmtmb asks for a ",
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
                                           draw_ids = NULL, ...) {
  re_form <- re_form_arg(re_formula, re.form, "posterior_predict()")
  check_flag(negative_rt, "negative_rt")
  if (negative_rt) {
    stop("posterior_predict(negative_rt = TRUE) is brms's sign ",
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
    stop("posterior_predict(): family '", rspec$family[["family"]],
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
    stop("posterior_predict(newdata =) is not supported for this ",
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
    stop("posterior_predict(re_formula =) is not supported for this ",
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
        dpv[[dnm]] <- as.vector(predict(sh, newdata = newdata,
                                        dpar = dnm, resp = resp,
                                        re.form = re_form,
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

#' @exportS3Method bayesplot::pp_check
#' @export
pp_check.frmtmb_draws <- function(object, type = "dens_overlay",
                                  ndraws = 50,
                                  re_formula = arg_unset(),
                                  re.form = arg_unset(), ...) {
  # the draws method's default is NULL, not the fit method's NA: a draw
  # already CARRIES its random effects, so conditioning on them is the
  # posterior predictive check, while the fit method has one point
  # estimate and has to simulate new levels to get a spread at all
  re_form <- re_form_arg(re_formula, re.form, "pp_check()")
  fit <- object$fit
  rspec <- single_response(fit, "pp_check()")
  y <- fit$frame[["y"]][[1L]]
  if (is.matrix(y)) {
    stop("pp_check() on draws supports vector responses", call. = FALSE)
  }
  yrep <- posterior_predict(object, ndraws = ndraws, re_formula = re_form)
  fun <- get(paste0("ppc_", type), envir = asNamespace("bayesplot"))
  fun(as.numeric(y), yrep, ...)
}

#' Convert draws to a posterior draws object
#'
#' @param x A `frmtmb_draws` object.
#' @param ... Unused.
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
as_draws.frmtmb_draws <- function(x, ...) {
  posterior::as_draws_matrix(x$draws)
}

#' @exportS3Method posterior::variables
#' @export
variables.frmtmb_draws <- function(x, ...) {
  colnames(x$draws)
}

#' @export
as.data.frame.frmtmb_draws <- function(x, ...) {
  as.data.frame(x$draws)
}

# ---- the draws matrix in other shapes --------------------------------

#' @rdname sample-as_draws
#' @export
as.array.frmtmb_draws <- function(x, ...) {
  # iterations x chains x parameters, the layout bayesplot's mcmc_*
  # functions and posterior's draws_array both read. frm_sample()
  # rbinds the chains in order, so the draws matrix is already
  # chain-major and reshapes without a permutation.
  m <- x$draws
  nc <- x$stanfit@sim$chains %||% 1L
  if (nc <= 1L || nrow(m) %% nc != 0L) nc <- 1L
  array(as.vector(m), c(nrow(m) %/% nc, nc, ncol(m)),
        dimnames = list(NULL, paste0("chain:", seq_len(nc)),
                        colnames(m)))
}

#' @rdname sample-as_draws
#' @exportS3Method posterior::as_draws_matrix
#' @export
as_draws_matrix.frmtmb_draws <- function(x, ...) {
  posterior::as_draws_matrix(x$draws)
}

#' @rdname sample-as_draws
#' @exportS3Method posterior::as_draws_array
#' @export
as_draws_array.frmtmb_draws <- function(x, ...) {
  # through as.array(), so the chains stay separate: a draws_array
  # built from the flattened matrix would claim one chain and every
  # convergence diagnostic computed on it would be wrong
  posterior::as_draws_array(as.array(x))
}

#' @rdname sample-as_draws
#' @exportS3Method posterior::as_draws_df
#' @export
as_draws_df.frmtmb_draws <- function(x, ...) {
  posterior::as_draws_df(as_draws_array(x))
}

#' @rdname sample-as_draws
#' @exportS3Method posterior::as_draws_list
#' @export
as_draws_list.frmtmb_draws <- function(x, ...) {
  posterior::as_draws_list(as_draws_array(x))
}

#' @rdname sample-as_draws
#' @exportS3Method posterior::as_draws_rvars
#' @export
as_draws_rvars.frmtmb_draws <- function(x, ...) {
  posterior::as_draws_rvars(as_draws_array(x))
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
    stop("as.mcmc() needs the 'coda' package; as_draws() and ",
         "as.array() give the same draws without it", call. = FALSE)
  }
  check_flag(combine_chains, "combine_chains")
  check_flag(inc_warmup, "inc_warmup")
  if (inc_warmup) {
    stop("as.mcmc(inc_warmup = TRUE) has nothing to include: the draws ",
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
  a <- as.array(x)
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
nvariables.frmtmb_draws <- function(x, ...) ncol(x$draws)

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
#' @param ... Unused.
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
#'   posterior_summary(ds, variable = c("Intercept", "x"))
#'   posterior_interval(ds, prob = 0.9, variable = "x")
#'   head(predictive_interval(ds))
#' }
#' }
#' @name sample-posterior_summary
NULL


#' @rdname sample-posterior_summary
#' @exportS3Method brms::posterior_summary
#' @export
posterior_summary.frmtmb_draws <- function(x, probs = c(0.025, 0.975),
                                           robust = FALSE,
                                           variable = NULL, ...) {
  # the generic dispatches to core's default method on a matrix
  posterior_summary(draws_columns(x, variable),
                    probs = probs, robust = robust)
}

#' The requested columns of the draws matrix, defaulting to the ones
#' `summary()` and `print()` show: the group-level modes are thousands
#' of columns on a large fit and nobody asked for them by writing
#' `posterior_summary(ds)`.
#'
#' @noRd
draws_columns <- function(x, variable = NULL) {
  m <- x$draws
  if (is.null(variable)) {
    keep <- setdiff(colnames(m),
                    c("lp__", grep("^b\\[", colnames(m), value = TRUE)))
    return(m[, keep, drop = FALSE])
  }
  miss <- setdiff(variable, colnames(m))
  if (length(miss)) {
    stop("variable = names ", paste(miss, collapse = ", "),
         ", which the draws do not contain. variables() lists what is ",
         "there; note the draws-side spelling drops parentheses ",
         "(Intercept, not (Intercept))", call. = FALSE)
  }
  m[, variable, drop = FALSE]
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
    stop("Argument 'pars' must be NA or a character vector.",
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
    stop(what, ": `variable` names variables and must be a character ",
         "vector; variables(x) lists what is there", call. = FALSE)
  }
  check_flag(regex, "regex")
  if (regex) return(draws_extract_pars(variable, all_pars))
  miss <- setdiff(variable, all_pars)
  if (length(miss)) {
    stop(what, ": variable = names ", paste(miss, collapse = ", "),
         ", which the draws do not contain. variables() lists what is ",
         "there; note the draws-side spelling drops parentheses ",
         "(Intercept, not (Intercept))", call. = FALSE)
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
  sel <- draws_select_variables(object, pars, variable, regex, fixed,
                                "posterior_interval()")
  check_probability(prob, "prob")
  m <- draws_columns(object, sel)
  a <- (1 - prob) / 2
  t(apply(m, 2L, stats::quantile, probs = c(a, 1 - a)))
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
  re_form <- re_form_arg(re_formula, re.form, "predictive_interval()")
  yrep <- posterior_predict(object, newdata = newdata, resp = resp,
                            re_formula = re_form, ndraws = ndraws)
  if (length(dim(yrep)) > 2L) {
    stop("predictive_interval() needs one predicted number per ",
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
  method <- match.arg(method, c("posterior_predict", "posterior_epred"))
  fit <- draws_base_fit(object)
  resp <- resp %||% names(fit$spec$responses)[1L]
  y <- draws_response_values(fit, resp, newdata, "predictive_error()")
  if (is.matrix(y)) {
    stop("predictive_error() needs a vector response; this one is a ",
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
    stop(what, " needs the observed response to subtract from, and ",
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
#' `coef()` is a posterior quantity, not a structural one: it summarizes
#' the per-group coefficients (fixed effects plus that group's own
#' random effects) over the draws, in the same nested shape
#' [frmtmb::coef.frmtmb_fit()] returns, with a `levels x statistics x
#' coefficients` array in place of each data frame. That is brms's
#' `coef.brmsfit` layout.
#'
#' @param object,x A `frmtmb_draws` from [frm_sample()].
#' @param ... Unused.
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
#'   coef(ds)$g[1:3, , "(Intercept)"]
#' }
#' }
#' @name draws-structure
NULL

#' @rdname draws-structure
#' @export
nobs.frmtmb_draws <- function(object, ...) {
  stats::nobs(draws_base_fit(object))
}

#' @rdname draws-structure
#' @export
formula.frmtmb_draws <- function(x, ...) {
  stats::formula(draws_base_fit(x))
}

#' @rdname draws-structure
#' @export
family.frmtmb_draws <- function(object, ...) {
  stats::family(draws_base_fit(object))
}

#' @rdname draws-structure
#' @export
getCall.frmtmb_draws <- function(x, ...) x$fit$call

#' @exportS3Method brms::ngrps
#' @rawNamespace S3method(lme4::ngrps,frmtmb_draws)
#' @export
ngrps.frmtmb_draws <- function(object, ...) {
  ngrps(draws_base_fit(object))
}

#' @rdname draws-structure
#' @export
coef.frmtmb_draws <- function(object, ...) {
  idx <- draws_par_index(object$fit)
  n <- nrow(object$draws)
  per <- lapply(seq_len(n), function(i) coef(draws_fit_at(object, i, idx)))
  draws_summarize_coef(per)
}

#' Turn a list of per-draw `coef()` results into brms's array shape.
#' `coef.frmtmb_fit()` returns either a nested list of data frames (per
#' dpar, per grouping factor) or, for a fit with no group-level terms, a
#' plain named vector; both are summarized here so the draws method
#' mirrors the fit method exactly.
#'
#' @noRd
draws_summarize_coef <- function(per) {
  first <- per[[1L]]
  if (!is.list(first)) {
    M <- do.call(rbind, per)
    return(cbind(Estimate = colMeans(M),
                 Est.Error = apply(M, 2L, stats::sd),
                 Q2.5 = apply(M, 2L, stats::quantile, 0.025),
                 Q97.5 = apply(M, 2L, stats::quantile, 0.975)))
  }
  if (is.data.frame(first)) {
    A <- vapply(per, function(d) as.matrix(d), as.matrix(first))
    st <- array(NA_real_, c(nrow(first), 4L, ncol(first)),
                dimnames = list(rownames(first),
                                c("Estimate", "Est.Error",
                                  "Q2.5", "Q97.5"),
                                colnames(first)))
    st[, "Estimate", ] <- apply(A, c(1, 2), mean)
    st[, "Est.Error", ] <- apply(A, c(1, 2), stats::sd)
    st[, "Q2.5", ] <- apply(A, c(1, 2), stats::quantile, 0.025)
    st[, "Q97.5", ] <- apply(A, c(1, 2), stats::quantile, 0.975)
    return(st)
  }
  out <- list()
  for (nm in names(first)) {
    out[[nm]] <- draws_summarize_coef(lapply(per, `[[`, nm))
  }
  out
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
#' All of these report the frmtmb draws-side parameter names (no
#' parentheses), not Stan's `par[1]`, except `nuts_params()`, whose
#' rows are the sampler's own quantities and not model parameters.
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
#'   modes and `lp__`. `rhat()` and `neff_ratio()` do not take it,
#'   because brms's do not: their selector is `pars`. Naming it on
#'   either of those two is silently ignored today; see `...`.
#' @param pars Which variables to report, in brms's spelling. The rule
#'   differs by method; see *Two different `pars` rules, both brms's*.
#' @param regex For `rhat()` and `neff_ratio()`, `TRUE` makes `pars` a
#'   regular expression; for `mcmc_plot()`, it makes `variable` one.
#' @param fixed For `mcmc_plot()`, `TRUE` matches `pars` by exact name
#'   rather than as a regular expression.
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
#'   mcmc_plot(ds, type = "trace", variable = "x")
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
  a <- as.array(object)
  keep <- sel %||%
    setdiff(dimnames(a)[[3L]],
            c("lp__", grep("^b\\[", dimnames(a)[[3L]], value = TRUE)))
  if (!length(keep)) {
    stop("mcmc_plot(): no variable was selected. variables(x) lists ",
         "what is there", call. = FALSE)
  }
  fun(a[, , keep, drop = FALSE], ...)
}

#' @rdname draws-diagnostics
#' @export
pairs.frmtmb_draws <- function(x, variable = NULL, ...) {
  # the same validation mcmc_plot() does, so a name that is not there
  # gets the message that lists variables() rather than base R's
  # "subscript out of bounds" out of the array index below
  sel <- draws_select_variables(x, NA, variable, FALSE, FALSE, "pairs()")
  a <- as.array(x)
  keep <- sel %||%
    utils::head(setdiff(dimnames(a)[[3L]],
                        c("lp__",
                          grep("^b\\[", dimnames(a)[[3L]],
                               value = TRUE))), 4L)
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
    stop(what, " needs the 'bayesplot' package; as.array() gives the ",
         "draws in the layout bayesplot expects if you would rather ",
         "call it yourself", call. = FALSE)
  }
  ns <- asNamespace("bayesplot")
  if (!exists(nm, envir = ns, inherits = FALSE)) {
    stop(what, " asks for bayesplot::", nm, "(), which does not exist. ",
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
    stop(what, " needs the 'posterior' package: it computes the ",
         "diagnostic brms reports, which is posterior's. The sampler's ",
         "own classic split-R-hat and n_eff are in `x$stanfit` without ",
         "it", call. = FALSE)
  }
  a <- as_draws_array(x)
  if (is.null(pars)) return(a)
  check_flag(regex, "regex")
  posterior::subset_draws(a, variable = pars, regex = regex)
}

#' bayesplot's namespace, or an error naming the accessor that wanted it.
#'
#' @noRd
draws_bayesplot_ns <- function(what) {
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop(what, " needs the 'bayesplot' package: it reads the sampler ",
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
    stop("pp_mixture(resp = \"", resp, "\") names no response of this ",
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
    stop("probs must be two numbers strictly between 0 and 1",
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
    stop(what, " does not take newdata. ", ..., call. = FALSE)
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
    stop(what, " does not take re_formula. ", ..., call. = FALSE)
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
#' @param object,x,pars,... Ignored; these functions always stop.
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
#'   try(nsamples(ds))
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
  stop("stancode() has no meaning for frmtmb: there is no Stan ",
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
  stop("standata() has no meaning for frmtmb: nothing is exported to a ",
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
  stop("expose_functions() has nothing to expose: brms compiles Stan ",
       "functions and this makes them callable from R, while a frmtmb ",
       "custom family is already plain R: the lpdf you passed to ",
       "custom_family() is an R function you can call directly",
       call. = FALSE)
}


#' @rdname frmtmb-draws-refusals
#' @export
plot.frmtmb_draws <- function(x, ...) {
  stop("plot() has no display for frmtmb draws: brms's default panel ",
       "is the trace-and-density view, which mcmc_plot(x) renders ",
       "here (mcmc_plot(x, type = \"trace\") for the traces alone)",
       call. = FALSE)
}

#' @rdname frmtmb-draws-refusals
#' @export
update.frmtmb_draws <- function(object, ...) {
  stop("update() has no method for draws: the sampled object carries ",
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
  stop("restructure() is brms's upgrade path for objects saved by an ",
       "older brms; frmtmb has no such conversion. A draws object from ",
       "an older frmtmb is re-created by re-running frm_sample()",
       call. = FALSE)
}

#' @rdname frmtmb-draws-refusals
#' @export
posterior_samples <- function(x, pars = NA, ...) {
  UseMethod("posterior_samples")
}

#' @rdname frmtmb-draws-refusals
#' @exportS3Method brms::posterior_samples
#' @rawNamespace S3method(gratia::posterior_samples,frmtmb_draws)
#' @export
posterior_samples.frmtmb_draws <- function(x, pars = NA, ...) {
  stop("posterior_samples() is the deprecated brms spelling. Use ",
       "as_draws(x) for a posterior draws_matrix, as.matrix(x) for a ",
       "plain matrix, or as.data.frame(x)", call. = FALSE)
}

#' @rdname frmtmb-draws-refusals
#' @export
nsamples <- function(object, ...) UseMethod("nsamples")

#' @rdname frmtmb-draws-refusals
#' @exportS3Method rstantools::nsamples
#' @export
nsamples.frmtmb_draws <- function(object, ...) {
  stop("nsamples() is the deprecated brms spelling. Use ndraws(x) for ",
       "the pooled draw count, niterations(x) for the per-chain count",
       call. = FALSE)
}

#' @rdname frmtmb-draws-refusals
#' @export
parnames <- function(x, ...) UseMethod("parnames")

#' @rdname frmtmb-draws-refusals
#' @exportS3Method brms::parnames
#' @export
parnames.frmtmb_draws <- function(x, ...) {
  stop("parnames() is the deprecated brms spelling. Use variables(x), ",
       "which lists the same names in the same draws-side spelling ",
       "(no parentheses)", call. = FALSE)
}
