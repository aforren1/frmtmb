# S3 accessor and print methods for frmtmb_fit.

#' Template-shaped Estimate / Std. Error lists. as.list(sdreport) fills the
#' full parameter template regardless of map or REML, which keeps indexing
#' by lp$idx valid; mapped-out entries get SE = NA.
#'
#' @noRd
par_est_se <- function(fit, vcov = NULL) {
  est <- fit$estimates
  if (!is.null(vcov)) {
    # a supplied covariance (vcov_cluster(), say) replaces the SEs of
    # the outer parameters only; the inner ones have none to replace
    rv <- resolve_vcov_arg(fit, vcov, "summary")
    om <- outer_par_map(fit)
    se <- lapply(est, function(v) rep(NA_real_, length(v)))
    sd_all <- sqrt(pmax(0, diag(rv$V)))
    for (cp in unique(om$comp)) {
      i <- om$comp == cp
      pos <- seq_along(est[[cp]])
      if (cp == "betad" && length(fit$frame[["betad_fixed_idx"]])) {
        pos <- setdiff(pos, fit$frame[["betad_fixed_idx"]])
      }
      se[[cp]][pos] <- sd_all[i]
    }
    return(list(est = est, se = se))
  }
  # A degenerate fit (no free outer parameters, no random effects) gives
  # sdreport an empty summary, and as.list() has no rows to reshape
  degenerate <- !length(fit$opt$par) && !length(fit$obj$env$random)
  se <- if (degenerate) {
    lapply(est, function(v) rep(NA_real_, length(v)))
  } else {
    as.list(sdr_of(fit), what = "Std. Error")
  }
  if (length(fit$frame[["betad_fixed_idx"]])) {
    se$betad[fit$frame[["betad_fixed_idx"]]] <- NA_real_
  }
  list(est = est, se = se)
}

#' The text of the ` Links:` line for one or more responses, or `""`
#' when no family has a link to name.
#'
#' A reader looks at a coefficient table and has to know which scale
#' the numbers are on. The family name settled that only while every
#' family had one link; it now does not. A multivariate fit prefixes
#' each response, because its families need not agree.
#'
#' @noRd
family_links_str <- function(responses) {
  s <- vapply(responses, function(r) family_link_str(r$family), "")
  # filter BEFORE the prefix: prefixing first makes a family with
  # nothing to say into the non-empty string "y: "
  keep <- nzchar(s)
  if (length(s) > 1L) s <- paste0(names(responses), ": ", s)
  s <- s[keep]
  if (!length(s)) "" else paste(s, collapse = "; ")
}

#' Print that text under the `Family:` line, wrapped.
#'
#' Wrapped because a mixture names every component's every dpar: two
#' gaussian-ish components run to 106 characters on one line, and a
#' third component is worse. Unwrapped, the terminal breaks it wherever
#' the width happens to land, in the middle of a name. `strwrap()`
#' breaks it between entries instead and indents the continuations
#' under the first one.
#'
#' @noRd
cat_family_links <- function(responses) {
  txt <- if (is.character(responses)) responses else {
    family_links_str(responses)
  }
  if (!nzchar(txt)) return(invisible(NULL))
  w <- max(40L, getOption("width", 80L))
  # break BETWEEN entries and never inside one. strwrap() splits on any
  # space, which puts "nu2" at the end of a line and "= logm1" at the
  # start of the next; a reader then cannot see which dpar that is.
  parts <- strsplit(txt, "; ", fixed = TRUE)[[1L]]
  parts <- paste0(parts, c(rep(";", length(parts) - 1L), ""))
  out <- character(0)
  cur <- " Links:"
  for (p in parts) {
    cand <- paste(cur, p)
    if (nchar(cand) > w && !identical(cur, " Links:")) {
      out <- c(out, cur)
      cur <- paste("       ", p)
    } else {
      cur <- cand
    }
  }
  cat(c(out, cur), sep = "\n")
  cat("\n")
  invisible(NULL)
}

#' The brms arguments the point-estimate accessors cannot answer.
#'
#' Named rather than reported as unknown, for the reason
#' `fitted_no_draws` gives: each of these is a REAL argument of the
#' `brmsfit` method of the same generic, so "there is no such argument"
#' sends the caller looking for a typo that is not there. Every one of
#' them was accepted in silence before the dots were refused.
#'
#' @noRd
brms_draws_summary_args <- c(
  summary = paste("brms summarizes posterior draws; this is one",
                  "estimate, so the return value is already what",
                  "`summary = TRUE` would produce. Draw from the fit",
                  "with frmtmb.sample::frm_sample() for the posterior"),
  robust = paste("a median and MAD over draws need draws.",
                 "frmtmb.sample's posterior_summary() takes robust"),
  probs = paste("quantiles over draws need draws. For an interval",
                "around a coefficient here use confint()"),
  pars = paste("brms's `pars` filters parameter names with a regular",
               "expression. Subset the returned value instead;",
               "variables() lists the names this package uses")
)

#' The brms arguments the two print methods and summary() cannot honor.
#'
#' @noRd
brms_print_args <- c(
  digits = paste("not implemented here yet. brms's print() takes it",
                 "and this one does not; summary(object) returns the",
                 "coefficient table, which round() accepts"),
  short = paste("brms's short display has no counterpart here.",
                "print() already prints the short form")
)

#' @noRd
brms_summary_args <- c(
  mc_se = paste("Monte Carlo standard errors describe a sampler, and",
                "a maximum likelihood fit has none")
)

# brms's print() of a fit IS its summary, and the two used to disagree
# here: print() showed the estimates alone while summary() showed the
# table. One renderer, so a section can only be added once.
#' @export
print.frmtmb_fit <- function(x, ...) {
  frm_check_dots(..., .unsupported = brms_print_args)
  require_fitted(x, "print()")
  print(summary(x))
  invisible(x)
}

#' The name a linear predictor's coefficient block gets in output. A
#' multivariate fit needs the response in the key to stay unambiguous;
#' a univariate fit shows the dpar alone.
#'
#' @noRd
coef_block_key <- function(fit, lp) {
  if (length(fit$spec$responses) > 1) {
    paste(lp[["resp"]], lp[["dpar"]], sep = "_")
  } else {
    lp[["dpar"]]
  }
}

# `vcov` takes a covariance over the whole outer parameter vector -
# vcov_cluster(full = TRUE), or a function of the fit returning one -
# and reports its standard errors in the coefficient table, with a t
# reference when the matrix carries degrees of freedom. Documented on
# the vcov_cluster() page; the variance components keep the
# model-based standard errors either way.
#' @export
summary.frmtmb_fit <- function(object, priors = FALSE, prob = 0.95,
                               robust = FALSE, mc_se = FALSE, ...,
                               vcov = NULL) {
  frm_check_dots(..., .unsupported = c(brms_summary_args,
                                       brms_draws_summary_args["pars"]))
  check_flag(priors, "priors")
  check_flag(robust, "robust")
  check_flag(mc_se, "mc_se")
  if (robust) {
    frm_stop("summary() cannot honor robust = TRUE: brms's robust summary ",
             "is the median and MAD of the draws, and a maximum-likelihood ",
             "fit has no draws. Sample with frmtmb.sample::frm_sample() and ",
             "summarize the draws for that", call. = FALSE)
  }
  if (mc_se) {
    frm_stop("summary() cannot honor mc_se = TRUE: ",
             brms_summary_args[["mc_se"]], call. = FALSE)
  }
  if (!is.numeric(prob) || length(prob) != 1L || is.na(prob) ||
      prob <= 0 || prob >= 1) {
    frm_stop("`prob` must be one probability strictly between 0 and 1, ",
             "not ", arg_desc(prob), call. = FALSE)
  }
  rdf <- NULL
  if (!is.null(vcov)) {
    # resolve once: `vcov` may be a function of the fit
    rv <- resolve_vcov_arg(object, vcov, "summary")
    vcov <- rv$V
    rdf <- rv$df
  }
  ps <- par_est_se(object, vcov)
  coefs <- list()
  for (lp in object$frame[["linpreds"]]) {
    est <- ps$est[[lp[["par"]]]][lp[["idx"]]]
    se <- ps$se[[lp[["par"]]]][lp[["idx"]]]
    z <- est / se
    cm <- if (is.null(rdf)) {
      cbind(Estimate = est, `Std. Error` = se, `z value` = z,
            `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
    } else {
      cbind(Estimate = est, `Std. Error` = se, `t value` = z,
            `Pr(>|t|)` = 2 * stats::pt(-abs(z), rdf))
    }
    rownames(cm) <- colnames(lp[["X"]])
    coefs[[coef_block_key(object, lp)]] <- cm
  }
  structure(
    list(call = object$call, family = family(object),
         # family(object) is the SINGLE family and is empty for a
         # multivariate fit, so the links are rendered here off the
         # responses, the way print.frmtmb_fit() does. Reading
         # x$family for them printed no Links line at all on the one
         # kind of fit whose families need not agree.
         links = family_links_str(object$spec$responses),
         formula = formula(object), nobs = stats::nobs(object),
         ngrps = ngrps(object),
         data_name = summary_data_name(object),
         group = names(ngrps(object) %||% list()),
         # the (MAP) marker says the fit was penalized by a prior, which
         # is the one thing "ML" would otherwise hide
         algorithm = paste0(if (object$REML) "REML" else "ML",
                            if (!is.null(object$prior)) " (MAP)"),
         prob = prob,
         loglik = logLik(object), AIC = stats::AIC(object),
         BIC = stats::BIC(object), REML = object$REML,
         importance = object$importance,
         # brms's slots. `fixed` and `random` are what a ported script
         # indexes; the frmtmb slots below them are unchanged.
         fixed = summary_fixed_frame(object, ps, rdf, prob),
         spec_pars = summary_spec_frame(object, prob),
         cor_pars = summary_cor_pars_frame(object, prob),
         random = summary_random_list(object, prob),
         # the FLAG is stored beside the value: a plain ML fit has no
         # priors, and a summary asked for them still says so rather
         # than dropping the section
         priors = priors,
         # prior_summary() says "no priors" by printing, which belongs
         # to the summary's own Priors section and not to the middle of
         # a summary() call
         prior = if (priors) {
           utils::capture.output(p <- prior_summary(object))
           p
         } else NULL,
         coefficients = coefs, varcor = varcorr_matrices(object),
         rescor = rescor_matrix(object),
         # R-side residual correlation, on the natural scale with the
         # same delta-method interval confint_varcorr() reports
         autocor = local({
           tr <- autocor_trans_rows(object)
           if (is.null(tr)) NULL else {
             z <- stats::qnorm(0.975)
             m <- cbind(
               Estimate = varcorr_untrans(tr$type, tr$est_t),
               `2.5 %` = varcorr_untrans(tr$type, tr$est_t - z * tr$se_t),
               `97.5 %` = varcorr_untrans(tr$type, tr$est_t + z * tr$se_t))
             rownames(m) <- if (length(unique(tr$block)) > 1L) {
               paste(tr$block, tr$term)
             } else tr$term
             attr(m, "label") <- tr$block[1L]
             m
           }
         }),
         smooth_edf = smooth_edf(object),
         extras = local({
           ex <- list()
           for (nm in object$frame[["extra_names"]] %||% character(0)) {
             cm <- cbind(Estimate = ps$est[[nm]],
                         `Std. Error` = ps$se[[nm]])
             rownames(cm) <- paste0(nm, "_", seq_len(nrow(cm)))
             ex[[nm]] <- cm
           }
           ex
         }),
         fixed_dpars = local({
           fx <- Filter(function(lp) !is.null(lp[["constant"]]),
                        object$frame[["linpreds"]])
           stats::setNames(vapply(fx, `[[`, numeric(1), "constant"),
                           vapply(fx, function(lp) {
                             coef_block_key(object, lp)
                           }, character(1)))
         })),
    class = "summary.frmtmb_fit"
  )
}

#' The name brms prints on its `Data:` line: the expression the fit was
#' given, not the frame itself.
#'
#' @noRd
summary_data_name <- function(object) {
  d <- object$call[["data"]]
  if (is.null(d)) return("")
  deparse1(d)
}

#' The two interval column names brms writes at coverage `prob`.
#'
#' @noRd
brms_ci_cols <- function(prob) {
  p <- format(prob * 100, trim = TRUE, scientific = FALSE,
              drop0trailing = TRUE)
  c(paste0("l-", p, "% CI"), paste0("u-", p, "% CI"))
}

#' brms's `$fixed`: the population-level coefficients, brms's four
#' columns first, then the Wald test this package reports and brms has
#' no counterpart for (brms puts `Rhat`, `Bulk_ESS` and `Tail_ESS`
#' there, and those describe a sampler).
#'
#' @noRd
summary_fixed_frame <- function(object, ps, rdf, prob) {
  rows <- brms_fixef_rows(object)
  est <- unname(brms_fixef_values(object, rows))
  se <- rep(NA_real_, length(rows$names))
  co <- !is.na(rows$idx)
  se[co] <- par_est_se_flat(ps, object)[rows$idx[co]]
  if (any(!co)) {
    # an ordinal fit's thresholds: their standard error is the delta
    # method of brms_fixef_extra_vcov(), which is model-based. A
    # `vcov` argument replaces the COEFFICIENT errors above and has no
    # threshold block to replace, so those rows keep the model-based
    # error rather than being dropped from the table.
    Vx <- tryCatch(suppressWarnings(brms_fixef_extra_vcov(object, rows)),
                   error = function(e) NULL)
    if (!is.null(Vx)) se[!co] <- sqrt(diag(Vx))[!co]
  }
  z <- est / se
  q <- if (is.null(rdf)) stats::qnorm(1 - (1 - prob) / 2) else {
    stats::qt(1 - (1 - prob) / 2, rdf)
  }
  out <- data.frame(Estimate = est, `Est.Error` = se,
                    lo = est - q * se, hi = est + q * se,
                    stat = z,
                    p = if (is.null(rdf)) 2 * stats::pnorm(-abs(z)) else {
                      2 * stats::pt(-abs(z), rdf)
                    },
                    check.names = FALSE)
  names(out) <- c("Estimate", "Est.Error", brms_ci_cols(prob),
                  if (is.null(rdf)) c("z value", "Pr(>|z|)") else {
                    c("t value", "Pr(>|t|)")
                  })
  rownames(out) <- rows$names
  out
}

#' Standard errors of the estimated coefficients in
#' `brms_coef_table()` order.
#'
#' @noRd
par_est_se_flat <- function(ps, object) {
  bd <- ps$se[["betad"]]
  if (length(fx <- object$frame[["betad_fixed_idx"]])) bd <- bd[-fx]
  unname(c(ps$se[["beta"]], bd))
}

#' brms's `$spec_pars`: a distributional parameter nobody wrote a
#' formula for, on its own natural scale.
#'
#' The interval is the coefficient's Wald interval mapped through the
#' link inverse, so it cannot leave the parameter's range; the standard
#' error is the delta method on the same map. brms reports the posterior
#' quantiles, which are inside the range for the same reason.
#'
#' @noRd
summary_spec_frame <- function(object, prob) {
  tab <- brms_coef_table(object)
  inv <- attr(tab, "linkinv")
  smp <- attr(tab, "simplex")
  in_smp <- unlist(lapply(smp, `[[`, "pos"))
  keep <- setdiff(which(tab$natural), in_smp)
  V <- tryCatch(suppressWarnings(vcov(object, full = TRUE)),
                error = function(e) NULL)
  cf <- fixef_estimated(object)
  se_in <- if (is.null(V)) rep(NA_real_, length(cf)) else {
    sqrt(diag(V))[seq_along(cf)]
  }
  q <- stats::qnorm(1 - (1 - prob) / 2)
  est <- err <- lo <- hi <- numeric(0)
  nm <- character(0)
  for (i in keep) {
    f <- inv[[i]]
    e0 <- f(cf[i])
    h <- 1e-5 * max(1, abs(cf[i]))
    d <- (f(cf[i] + h) - f(cf[i] - h)) / (2 * h)
    est <- c(est, e0)
    err <- c(err, abs(d) * se_in[i])
    lo <- c(lo, f(cf[i] - q * se_in[i]))
    hi <- c(hi, f(cf[i] + q * se_in[i]))
    nm <- c(nm, tab$brms[i])
  }
  for (s in smp %||% list()) {
    p <- s$to_simplex(cf[s$pos])
    for (k in seq_along(s$names)) {
      est <- c(est, p[1L, k])
      err <- c(err, NA_real_)
      lo <- c(lo, NA_real_)
      hi <- c(hi, NA_real_)
      nm <- c(nm, s$names[k])
    }
  }
  if (!length(nm)) return(summary_empty_block(prob))
  out <- data.frame(Estimate = est, `Est.Error` = err, lo = lo, hi = hi,
                    check.names = FALSE)
  names(out) <- c("Estimate", "Est.Error", brms_ci_cols(prob))
  rownames(out) <- nm
  # a simplex is set as a whole, so its interval is not a one-parameter
  # Wald interval; the estimates are reported and the bounds are NA
  out
}

#' brms's `$cor_pars`: the residual correlation structure's own
#' parameters, on their natural scale under brms's names.
#'
#' @noRd
summary_cor_pars_frame <- function(object, prob) {
  tr <- autocor_trans_rows(object)
  if (is.null(tr)) return(summary_empty_block(prob))
  summary_nat_frame(tr, prob, tr$term)
}

#' One natural-scale table from transformed-scale rows: the estimate,
#' its delta-method standard error, and the transformed-scale Wald
#' interval mapped back.
#'
#' @noRd
summary_nat_frame <- function(tr, prob, rn) {
  q <- stats::qnorm(1 - (1 - prob) / 2)
  est <- varcorr_untrans(tr$type, tr$est_t)
  out <- data.frame(
    Estimate = est,
    `Est.Error` = varcorr_nat_deriv(tr$type, tr$est_t) * tr$se_t,
    lo = varcorr_untrans(tr$type, tr$est_t - q * tr$se_t),
    hi = varcorr_untrans(tr$type, tr$est_t + q * tr$se_t),
    check.names = FALSE)
  names(out) <- c("Estimate", "Est.Error", brms_ci_cols(prob))
  rownames(out) <- rn
  out
}

#' Derivative of `varcorr_untrans()` at the transformed estimate, for
#' the delta-method standard error on the natural scale.
#'
#' @noRd
varcorr_nat_deriv <- function(type, v) {
  p <- 1 / (1 + exp(-v))
  ifelse(type == "raw", 1,
         ifelse(type == "cor", 1 - tanh(v)^2,
                ifelse(type == "prop", p * (1 - p), exp(v))))
}

#' brms's `$random`: one data frame per GROUPING FACTOR, with
#' `sd(<term>)` and `cor(<t1>,<t2>)` rows under brms's term names.
#'
#' `varcorr_trans_rows()` returns the random-effect blocks followed by
#' the residual correlation rows, which belong in `$cor_pars`; the tail
#' is dropped by count rather than by label, because a block label and
#' an autocorrelation label are both user text and could collide.
#'
#' @noRd
summary_random_list <- function(object, prob) {
  bks <- object$frame[["re_blocks"]]
  # brms gives NULL, not an empty list, when there is nothing here; a
  # ported script writes `if (is.null(s$random))` and `list()` is not
  # NULL (measured, dev/shapes-rev-brmsref.rds)
  if (!length(bks)) return(NULL)
  tr <- tryCatch(suppressWarnings(varcorr_trans_rows(object)),
                 error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  acr <- autocor_trans_rows(object)
  if (!is.null(acr)) tr <- utils::head(tr, nrow(tr) - nrow(acr))
  out <- list()
  for (bk in bks) {
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp", "car", "spde")) {
      next
    }
    rows <- which(tr$block == bk[["term_label"]])
    if (!length(rows)) next
    tn <- brms_re_rnames(object, bk)
    lab <- character(length(rows))
    isd <- tr$type[rows] == "sd"
    lab[isd] <- paste0("sd(", tn[seq_len(sum(isd))], ")")
    if (any(!isd)) {
      pairs <- which(lower.tri(diag(bk[["dim"]])), arr.ind = TRUE)
      lab[!isd] <- paste0("cor(", tn[pairs[, 2]], ",", tn[pairs[, 1]], ")")
    }
    g <- brms_group_name(bk)
    df <- summary_nat_frame(tr[rows, , drop = FALSE], prob, lab)
    out[[g]] <- if (is.null(out[[g]])) df else rbind(out[[g]], df)
  }
  if (!length(out)) NULL else out
}

#' The empty version of a brms summary block: brms gives a frame with
#' its own columns and no rows, not NULL, when a model has no
#' `spec_pars` or no `cor_pars` (measured, `dev/shapes-rev-brmsref.rds`
#' shows `$cor_pars` as a 0-row frame on a plain gaussian fit). A
#' ported script reads `nrow()` on it.
#'
#' @noRd
summary_empty_block <- function(prob) {
  out <- data.frame(Estimate = numeric(0), `Est.Error` = numeric(0),
                    lo = numeric(0), hi = numeric(0), check.names = FALSE)
  names(out) <- c("Estimate", "Est.Error", brms_ci_cols(prob))
  out
}

#' Print one of brms's summary blocks, rounded as brms rounds it.
#'
#' @noRd
print_summary_block <- function(df, digits = 2) {
  d <- as.data.frame(df, check.names = FALSE)
  for (j in seq_along(d)) {
    if (is.numeric(d[[j]])) {
      d[[j]] <- if (grepl("^Pr[(]", names(d)[j])) {
        format.pval(d[[j]], digits = digits)
      } else {
        round(d[[j]], digits)
      }
    }
  }
  print(d)
  invisible(NULL)
}

#' @export
print.summary.frmtmb_fit <- function(x, ...) {
  frm_check_dots(..., .unsupported = brms_print_args)
  cat(" Family:", x$family[["family"]], "\n")
  cat_family_links(x$links %||% family_link_str(x$family))
  cat("Formula:", deparse1(x$formula), "\n")
  cat("   Data:", x$data_name,
      paste0("(Number of observations: ", x$nobs, ")"), "\n")
  cat(" Method:", x$algorithm,
      "  logLik:", format(as.numeric(x$loglik), digits = 6),
      "  AIC:", format(x$AIC, digits = 6),
      "  BIC:", format(x$BIC, digits = 6), "\n")
  if (!is.null(x$importance)) {
    cat(imp_report_line(x$importance), "\n")
  }
  if (length(x$random)) {
    cat("\nMultilevel Hyperparameters:\n")
    for (g in names(x$random)) {
      # $random is keyed by brms's group name and ngrps() by frmtmb's,
      # which differ where brms's renaming rewrites one; an unmatched
      # key reports the count as unknown rather than printing nothing
      n_g <- x$ngrps[[g]] %||% NA_integer_
      cat("~", g, " (Number of levels: ", n_g, ") \n", sep = "")
      print_summary_block(x$random[[g]])
    }
  }
  # an empty block is a 0-row frame now, brms's shape, so the test
  # is on the rows and not on NULL
  if (NROW(x$cor_pars)) {
    cat("\nCorrelation Structures:\n")
    print_summary_block(x$cor_pars)
  }
  if (!is.null(x$smooth_edf)) {
    cat("\nSmoothing Spline Hyperparameters (edf of the penalized part):\n")
    print(round(x$smooth_edf, 2))
  }
  if (length(x$varcor_special %||% list())) {
    cat("\nGaussian Process Terms:\n")
    print_summary_block(x$varcor_special)
  }
  cat("\nRegression Coefficients:\n")
  print_summary_block(x$fixed)
  if (NROW(x$spec_pars)) {
    cat("\nFurther Distributional Parameters:\n")
    print_summary_block(x$spec_pars)
  }
  for (i in seq_along(x$fixed_dpars)) {
    cat("\nFixed dpar: ", names(x$fixed_dpars)[i], " = ",
        x$fixed_dpars[i], "\n", sep = "")
  }
  if (!is.null(x$rescor)) {
    cat("\nResidual correlation:\n")
    print(signif(x$rescor, 4))
  }
  if (isTRUE(x$priors)) {
    cat("\nPriors:\n")
    if (is.null(x$prior)) {
      cat("No priors were set (plain maximum likelihood).\n")
    } else {
      print(x$prior)
    }
  }
  invisible(x)
}

#' Estimated outer parameters: profiled betas (control profile = TRUE)
#' leave opt$par but are still estimated and must count toward df.
#'
#' @noRd
n_outer_est <- function(object) {
  n <- length(object$opt$par)
  if (isTRUE(object$control$profile)) {
    n <- n + length(object$frame[["par_template"]][["beta"]])
  }
  n
}

#' @export
logLik.frmtmb_fit <- function(object, ...) {
  frm_check_dots(...)
  require_fitted(object, "logLik() (and AIC(), BIC(), anova())")
  structure(-object$opt$objective,
            df = n_outer_est(object),
            nobs = object$frame[["n_obs"]],
            REML = object$REML,
            class = "logLik")
}

#' @export
nobs.frmtmb_fit <- function(object, ...) {
  # `use.fallback` is not named here: it is in `s3_contract_args`, with
  # every other name R's own machinery passes through a generic.
  frm_check_dots(..., .unsupported = c(resp = paste(
                   "a multivariate fit here shares one set of rows, so",
                   "every response has the same nobs()")))
  object$frame[["n_obs"]]
}

#' @export
df.residual.frmtmb_fit <- function(object, ...) {
  frm_check_dots(...)
  object$frame[["n_obs"]] - n_outer_est(object)
}

#' @export
family.frmtmb_fit <- function(object, ...) {
  frm_check_dots(..., .unsupported = c(resp = paste(
    "family() returns a NAMED LIST of families for a multivariate fit,",
    "so index it by response name instead")))
  fams <- lapply(object$spec$responses, `[[`, "family")
  if (length(fams) == 1) fams[[1]] else fams
}

#' Estimated residual correlation matrix (rescor fits), else NULL
#' @param fit A `frmtmb_fit`.
#' @return A correlation matrix or `NULL`.
#' @examples
#' set.seed(2)
#' n <- 80
#' dd <- data.frame(x = rnorm(n))
#' # two responses that share a residual disturbance
#' e <- rnorm(n)
#' dd$y1 <- 1 + 0.5 * dd$x + e + rnorm(n, 0, 0.5)
#' dd$y2 <- 2 - 0.3 * dd$x + e + rnorm(n, 0, 0.5)
#' fit <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE) + gaussian(),
#'            data = dd)
#' rescor_matrix(fit)
#'
#' # a fit without rescor has no residual correlation to report
#' fit0 <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + gaussian(), data = dd)
#' rescor_matrix(fit0)
#' @export
rescor_matrix <- function(fit) {
  if (inherits(fit, "frmtmb_draws")) {
    frm_stop("rescor_matrix() reads the fitted point estimate, so it ",
             "takes the frmtmb_fit, not draws: rescor_matrix(ds$fit). For ",
             "the posterior of the correlation, subset_draws() on the ",
             "rescor columns of as_draws(ds)", call. = FALSE)
  }
  if (!isTRUE(fit$spec$rescor)) return(NULL)
  K <- length(fit$spec$responses)
  C <- us_chol_cor(fit$estimates[["thetar"]], K)
  dimnames(C) <- list(names(fit$spec$responses),
                      names(fit$spec$responses))
  C
}

#' @export
formula.frmtmb_fit <- function(x, ...) {
  frm_check_dots(...)
  if (inherits(x$bform, "frmtmb_mvformula")) {
    x$bform$forms[[1]]$formula
  } else {
    x$bform$formula
  }
}

#' Names of the estimated (non-mapped) fixed coefficients, in template order.
#'
#' @noRd
estimated_coef_names <- function(fit) {
  tpl <- fit$frame[["par_template"]]
  nm_betad <- names(tpl[["betad"]])
  if (length(fit$frame[["betad_fixed_idx"]])) {
    nm_betad <- nm_betad[-fit$frame[["betad_fixed_idx"]]]
  }
  c(names(tpl[["beta"]]), nm_betad)
}

#' Covariance matrix of the fixed-effect estimates
#'
#' Covers brms's population-level coefficients and names its rows as
#' brms does (`Intercept`, `sigma_Intercept`, `x`). A distributional
#' parameter nobody wrote a formula for is not one of them: frmtmb
#' estimates `sigma` or `nu` as an intercept-only linear predictor,
#' brms reports it as a parameter of its own, and `summary()` shows it
#' under `Further Distributional Parameters`. Dpars fixed to constants
#' are excluded too.
#'
#' `full = TRUE` is the joint covariance of the whole outer parameter
#' vector on its internal scale: the fixed-effect coefficients, the
#' covariance parameters `theta` (log standard deviations, Fisher-z
#' correlations, and whatever else a structure keeps there), and any
#' extra parameters such as the ordinal thresholds. It is the matrix a
#' delta-method calculation on a variance component needs, and it is
#' what [hypothesis()] uses for `method = "wald"` - so an ICC or a
#' heritability is usually easier to ask for through `hypothesis()`,
#' which names the components for you, than to assemble by hand from
#' this matrix.
#'
#' Under `REML = TRUE` (or `frmtmb_control(profile = TRUE)`) the fixed
#' effects are integrated out of the outer problem, so they are not
#' part of `full = TRUE`; the block comes from the joint precision and
#' carries exactly the parameters [confint.frmtmb_fit()] reports.
#' `vcov(object)`
#' is still the fixed-effect covariance there.
#'
#' Passing `cluster` forwards to [vcov_cluster()] for the
#' cluster-robust (sandwich) covariance, in the `sandwich::vcovCL()`
#' spelling: `vcov(fit, cluster = ~ g, type = "CR1")`.
#'
#' @param object A `frmtmb_fit`.
#' @param correlation If `TRUE`, the correlation matrix instead, as in
#'   brms.
#' @param pars Row and column names to keep, in the order given, as in
#'   brms.
#' @param full If `TRUE`, include covariance parameters (`theta`),
#'   named as in `confint()` (the glmmTMB `vcov(full = TRUE)`
#'   convention). `full = TRUE` keeps the INTERNAL names, because it is
#'   the matrix a delta-method calculation on `confint()`'s rows needs;
#'   the default block takes brms's.
#' @param cluster Optional clustering factor. When given, the result is
#'   [vcov_cluster()]'s cluster-robust covariance instead of the
#'   model-based one.
#' @param type Small-sample correction for `cluster`, see
#'   [vcov_cluster()].
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return A covariance matrix.
#' @seealso [confint_varcorr()] for natural-scale intervals on the same
#'   covariance parameters, and [hypothesis()] for delta-method tests of
#'   expressions in them.
#'
#' @srrstats {RE4.6} The variance-covariance matrix of the model
#'   parameters is returned by `vcov()`: the fixed-effect block by
#'   default, and the covariance parameters as well under `full = TRUE`,
#'   named as in `confint()`. It comes from the inverse observed
#'   information produced by `RTMB::sdreport()`, or from the joint
#'   precision for a REML or profiled fit. A covariance that could not be
#'   recovered from the Hessian warns rather than returning silent `NaN`.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # standard errors of the fixed effects, which is fixef()'s Est.Error
#' sqrt(diag(vcov(fit)))
#' # the covariance parameters join the block on their internal scale,
#' # under the internal names, which is confint()'s vocabulary
#' rownames(vcov(fit, full = TRUE))
#'
#' # the matrix is what a delta-method calculation needs
#' V <- vcov(fit)
#' a <- c(1, 2)                       # prediction at x = 2, no group
#' sqrt(drop(t(a) %*% V[1:2, 1:2] %*% a))
#' @export
vcov.frmtmb_fit <- function(object, correlation = FALSE, pars = NULL, ...,
                            full = FALSE, cluster = NULL, type = "CR0") {
  frm_check_dots(...)
  check_flag(correlation, "correlation")
  V <- if (is.null(cluster)) {
    vcov_estimated(object, full = full)
  } else {
    vcov_cluster(object, cluster, type = type, full = full)
  }
  if (full) return(V)
  V <- vcov_brms_block(object, V, model_based = is.null(cluster))
  if (correlation) V <- stats::cov2cor(V)
  brms_pars_filter(t(brms_pars_filter(t(V), pars, "vcov()")), pars,
                   "vcov()")
}

#' The covariance of EVERY estimated coefficient, under the internal
#' names of `estimated_coef_names()`.
#'
#' This is what `vcov()` returned before item 2.6f moved it to brms's
#' population-level block. It stayed, because the interop seams need a
#' covariance that lines up with `get_coef()` and
#' `insight::get_parameters()`, both of which report every estimated
#' coefficient: marginaleffects builds a numeric jacobian by perturbing
#' that vector, so a covariance one row shorter than it is not merely
#' differently named, it is the wrong matrix.
#'
#' @noRd
vcov_estimated <- function(object, full = FALSE) {
  nm <- estimated_coef_names(object)
  if (!object$REML && !isTRUE(object$control$profile)) {
    V <- sdr_of(object)$cov.fixed
    # This branch reads an already-inverted covariance and used to
    # return whatever sdreport put there. When sdreport could not
    # invert the outer Hessian that is a matrix of NaN, on a fit whose
    # own convergence checks all passed - so say so here, or nothing
    # does. (The REML/profile branch below warns through
    # solve_joint_precision().)
    if (any(!is.finite(V))) warn_nonfinite_cov(object$cache, object)
    if (full) {
      # cov.fixed rows repeat the component names; the per-parameter
      # names (confint rows) are the useful labels
      dimnames(V) <- list(outer_par_names(object),
                          outer_par_names(object))
      return(V)
    }
    ord <- c(which(rownames(V) == "beta"), which(rownames(V) == "betad"))
    V <- V[ord, ord, drop = FALSE]
  } else {
    Q <- sdr_of(object)$jointPrecision
    Vall <- solve_joint_precision(Q, object$cache, object)
    rn <- rownames(Q)
    if (full) {
      # The outer parameter vector under REML (or control profile =
      # TRUE) does not contain beta: it is integrated out. So
      # full = TRUE returns exactly the parameters confint() reports,
      # which is what the naming invariant asks for, and the block
      # comes out of the joint precision rather than cov.fixed. Use
      # vcov(object) for the fixed-effect covariance.
      comps <- setdiff(names(object$frame[["par_template"]]),
                       c("b", "miss", "beta"))
      keep <- unlist(lapply(comps, function(cp) which(rn == cp)))
      onm <- outer_par_names(object)
      if (length(keep) == length(onm)) {
        Vf <- as.matrix(Vall[keep, keep, drop = FALSE])
        dimnames(Vf) <- list(onm, onm)
        return(Vf)
      }
      frm_warning("full = TRUE could not align the joint-precision blocks ",
                  "with the outer parameter names; returning the ",
                  "fixed-effect block", call. = FALSE)
    }
    ord <- c(which(rn == "beta"), which(rn == "betad"))
    V <- as.matrix(Vall[ord, ord, drop = FALSE])
  }
  dimnames(V) <- list(nm, nm)
  V
}

#' The population-level sub-block of a coefficient covariance, named as
#' brms names it.
#'
#' brms's `vcov()` covers the coefficients its `fixef()` reports, so a
#' distributional parameter nobody wrote a formula for is NOT in it: it
#' is a `spec_par` reported on its own natural scale. frmtmb estimates
#' such a parameter as an intercept-only linear predictor, which is why
#' this block used to be one row wider than brms's on the same model.
#'
#' @noRd
vcov_brms_block <- function(object, V, model_based = TRUE) {
  rows <- brms_fixef_rows(object)
  if (length(rows$extra)) {
    # an ordinal fit's thresholds and cs() coefficients are brms
    # population-level parameters and belong in this matrix, but they
    # are not rows of V: they come from the joint covariance through a
    # delta method. A cluster-robust V has no such joint matrix behind
    # it, so there the block is the coefficient rows it does have and
    # the omission is said out loud rather than left to be read off a
    # dimension.
    if (model_based) {
      out <- brms_fixef_extra_vcov(object, rows)
      if (!is.null(out)) {
        if (!is.null(df <- attr(V, "df"))) attr(out, "df") <- df
        return(out)
      }
      frm_warning("the ordinal thresholds could not be aligned with the ",
                  "joint covariance; vcov() reports the coefficient rows ",
                  "only", call. = FALSE)
    } else {
      frm_warning("a cluster-robust covariance covers the estimated ",
                  "coefficients only, so the ordinal threshold rows are ",
                  "not in it; use vcov(object) for those", call. = FALSE)
    }
    rows$names <- rows$names[!is.na(rows$idx)]
    rows$idx <- rows$idx[!is.na(rows$idx)]
  }
  out <- V[rows$idx, rows$idx, drop = FALSE]
  dimnames(out) <- list(rows$names, rows$names)
  # a cluster-robust matrix carries its degrees of freedom, which the
  # coefficient table reads to choose a t reference; subsetting a matrix
  # drops every attribute but dim and dimnames
  if (!is.null(df <- attr(V, "df"))) attr(out, "df") <- df
  out
}

#' Per-group coefficients (fixed effects plus conditional modes)
#'
#' Follows the lme4/glmmTMB/brms convention: for each random-effect
#' grouping factor, the fixed effects of its linear predictor broadcast
#' over the group levels, with the conditional modes added to the
#' matching columns. Use [fixef()] for the fixed effects alone.
#'
#' The result is a list of data frames keyed by grouping factor. When
#' random effects appear in more than one dpar (or response), an outer
#' layer keyed like [fixef()] is added. Smooth terms are excluded. A fit
#' without random effects returns the coefficient vector of the location
#' predictor (when there is one linear predictor), or one vector per
#' predictor.
#'
#' The coefficients are named as brms names them, which is how
#' [fixef()] and [vcov()] name their rows: `Intercept`, not
#' `(Intercept)`. Anything that pairs `coef()` with `vcov()` by name,
#' `lmtest::coeftest()` among them, needs the two to agree. Use
#' [fixef_by_dpar()] for the design-column spelling.
#'
#' @param object A `frmtmb_fit`.
#' @param summary,robust,probs brms's arguments, in brms's
#'   positions so that a positional brms call asks the same question.
#'   brms answers `summary = FALSE` with the posterior draws and
#'   `robust = TRUE` with their median and MAD, and a maximum-likelihood
#'   fit has no draws, so both are refused by name with the reason. The
#'   default of each is accepted and changes nothing.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return A named list of data frames, one per grouping factor, each
#'   with one row per group level and one column per coefficient. When
#'   random effects appear in more than one linear predictor, the list is
#'   nested one level deeper, keyed as in [fixef()]. A fit without random
#'   effects returns the per-predictor coefficient vectors instead.
#'
#' @srrstats {RE4.2} Model coefficients are returned by `coef()`, in the
#'   lme4 and glmmTMB sense of per-group coefficients (fixed effects
#'   broadcast over the group levels with the conditional modes added),
#'   and by [fixef()] for the fixed effects alone. [ranef()] returns the
#'   conditional modes and [VarCorr()] the variance components.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # one row per group: the fixed effects with the modes added in
#' head(coef(fit)$g)
#' # which is fixef() plus ranef(), the lme4 identity, under the same
#' # names fixef() and vcov() use
#' all.equal(coef(fit)$g[["Intercept"]],
#'           fixef(fit)["Intercept", "Estimate"] + ranef(fit)$g[, 1],
#'           check.attributes = FALSE)
#'
#' # without random effects there are no groups, so coef() is the
#' # coefficient vector of the location predictor
#' coef(frm(bf(y ~ x) + gaussian(), data = dd))
#' @export
coef.frmtmb_fit <- function(object, summary = TRUE, robust = FALSE,
                            probs = c(0.025, 0.975), ...) {
  frm_check_dots(..., .unsupported = brms_draws_summary_args["pars"])
  # brms's positions, so coef(fit, FALSE) asks brms's question and gets
  # a refusal rather than an argument landing nowhere
  fit_refuse_draws_args("coef()", summary = summary, robust = robust,
                        probs = probs)
  # brms's coefficient names, the same ones fixef() and vcov() put on
  # their rows. They used to be the DESIGN COLUMN names, and the split
  # was not cosmetic: lmtest::coeftest() intersects names(coef()) with
  # rownames(vcov()), so `(Intercept)` here against `Intercept` there
  # dropped the intercept row out of a printed significance table with
  # no warning. brms names both the same way.
  tab <- brms_coef_table(object)
  fe <- list()
  cmap <- list()
  for (lp in object$frame[["linpreds"]]) {
    key <- coef_block_key(object, lp)
    v <- object$estimates[[lp[["par"]]]][lp[["idx"]]]
    bn <- brms_lp_coef_names(lp, tab)
    names(v) <- bn
    fe[[key]] <- v
    # the conditional modes are keyed by DESIGN column, which is the
    # vocabulary the random-effect blocks carry, so the columns they
    # are added to have to be looked up
    cmap[[key]] <- stats::setNames(bn, colnames(lp[["X"]]))
  }
  # an ordinal fit's thresholds and cs() coefficients are population-
  # level parameters with no design column of their own, and brms's
  # coef() carries them; without them vcov() had rows coef() did not
  rows <- brms_fixef_rows(object)
  for (e in rows$extra) {
    fe[[e$key]] <- c(fe[[e$key]], stats::setNames(e$values, e$names))
  }
  for (k in names(fe)) {
    fe[[k]] <- fe[[k]][order(match(names(fe[[k]]), rows$names))]
  }
  cvec <- coef_b(object)
  out <- list()
  for (bk in object$frame[["re_blocks"]]) {
    if (bk[["covstruct"]] == "smooth") next
    bmat <- t(matrix(cvec[bk[["c_idx"]]], nrow = bk[["dim"]]))
    for (cp in bk[["components"]]) {
      lp <- object$frame[["linpreds"]][[cp$lp_key]]
      key <- coef_block_key(object, lp)
      # a second term on the same factor adds its modes to the same
      # frame; the fixed effects are broadcast only once
      gname <- bk[["group_name"]]
      df <- out[[key]][[gname]]
      if (!is.null(df) && !identical(rownames(df), bk[["levels"]])) {
        gname <- bk[["term_label"]]
        df <- out[[key]][[gname]]
      }
      if (is.null(df)) {
        fev <- fe[[key]]
        df <- as.data.frame(
          matrix(fev, nrow = bk[["n_levels"]], ncol = length(fev),
                 byrow = TRUE, dimnames = list(bk[["levels"]], names(fev))),
          optional = TRUE
        )
      }
      thr <- Filter(function(e) {
        identical(e$comp, "tau_raw") && identical(e$key, key)
      }, rows$extra)
      for (j in seq_len(cp$dim)) {
        cn <- cp$cnms[j]
        bv <- bmat[, cp$offset + j]
        if (length(thr) && identical(cn, "(Intercept)")) {
          # an ordinal predictor has no intercept column: its intercept
          # IS the thresholds, so a random intercept moves each of them,
          # as brms:::coef.brmsfit moves them. The sign is the one the
          # family's likelihood gives eta against a threshold. That was
          # a stray `(Intercept)` column holding the mode alone beside
          # thresholds repeated unchanged across groups.
          df <- coef_shift_thresholds(df, thr[[1L]]$names, bv,
                                      brms_lp_family(object, lp))
          next
        }
        # a random-effect column the fixed part does not have (an
        # `0 + x | g` whose x is not a fixed term) takes the name brms
        # gives the group-level coefficient
        if (!is.na(mapped <- cmap[[key]][cn])) {
          cn <- unname(mapped)
        } else {
          cn <- brms_usc(brms_lp_prefix(object, lp), brms_rename(cn))
        }
        if (cn %in% colnames(df)) {
          df[[cn]] <- df[[cn]] + bv
        } else {
          df[[cn]] <- bv
        }
      }
      out[[key]][[gname]] <- df
    }
  }
  if (!length(out)) {
    # GLM-style fits: the mu vector alone, unless another dpar is
    # actually modeled with covariates
    if (length(object$spec$responses) == 1L && "mu" %in% names(fe)) {
      aux <- Filter(function(lp) lp[["dpar"]] != "mu",
                    object$frame[["linpreds"]])
      simple <- all(vapply(aux, function(lp) {
        !is.null(lp[["constant"]]) ||
          (ncol(lp[["X"]]) == 1L && identical(colnames(lp[["X"]]),
                           "(Intercept)"))
      }, TRUE))
      if (simple) return(fe[["mu"]])
    }
    return(fe)
  }
  if (length(out) == 1L) out[[1L]] else out
}

#' Move an ordinal predictor's thresholds by a group's random intercept,
#' with brms's sign.
#'
#' brms:::coef.brmsfit keys the sign on the family's specials:
#' `thres_minus_eta` families (cumulative, sratio) evaluate
#' `threshold - eta`, so the group's threshold is the threshold minus
#' the group's mode; `eta_minus_thres` families (cratio, acat) evaluate
#' `eta - threshold`, and brms reports the mode minus the threshold
#' there. frmtmb's four ordinal likelihoods use the same two forms
#' (`ord_cat_probs()`). Any other family gets brms's default, the
#' threshold plus the mode.
#'
#' @noRd
coef_shift_thresholds <- function(df, nms, bv, fam) {
  f <- fam[["family"]] %||% ""
  for (nm in nms) {
    df[[nm]] <- if (f %in% c("cumulative", "sratio")) {
      df[[nm]] - bv
    } else if (f %in% c("cratio", "acat")) {
      bv - df[[nm]]
    } else {
      df[[nm]] + bv
    }
  }
  df
}

#' Extract fixed effects
#'
#' Two shapes, and the difference is the NAMES. The default is brms's
#' summary matrix: one row per population-level coefficient under
#' brms's own name (`Intercept`, `sigma_Intercept`, `x`), and the
#' columns `Estimate`, `Est.Error`, `Q2.5` and `Q97.5`. A
#' maximum-likelihood fit has no draws, so `Est.Error` is the standard
#' error and the `Q` columns are the Wald interval at those
#' probabilities, which is the interval [confint()] reports.
#' `flatten = TRUE` is one vector of estimates in the INTERNAL parameter
#' name, which is what `confint()`, `par_template()`, `start` and
#' `newparams` use: `dpar_column`, with the location parameter's own
#' coefficients left unprefixed (`x`, not `mu_x`), and the response
#' prefixed ahead of that in a multivariate fit.
#'
#' A distributional parameter nobody wrote a formula for is not a
#' population-level coefficient in brms and is not a row here: `sigma`
#' of a plain gaussian fit is reported by `summary()` under
#' `Further Distributional Parameters`, on its own natural scale.
#'
#' `hypothesis()` and `variables()` use a THIRD vocabulary, brms's
#' parameter names: `b_Intercept`, `b_sigma_Intercept` for a `sigma`
#' formula, and `sigma` on its natural scale when no formula was written
#' for it. `hypothesis()` puts `class = "b"`'s `b_` in front of a bare
#' name, so its default spelling of those is `Intercept` and
#' `sigma_Intercept`.
#'
#' `unlist(fixef(fit, flatten = TRUE))` is none of them. Use
#' `flatten = TRUE` to line coefficients up with a prior or with
#' `confint()`, and the default matrix to read an estimate with its
#' standard error.
#'
#' @param object A `frmtmb_fit`.
#' @param summary,robust brms's arguments, in brms's positions so that
#'   a positional brms call asks the same question. brms answers
#'   `summary = FALSE` with the posterior draws and `robust = TRUE` with
#'   their median and MAD, and a maximum-likelihood fit has no draws, so
#'   both are refused by name with the reason. The default of each is
#'   accepted and changes nothing.
#' @param probs Probabilities of the two quantile columns. brms takes
#'   the posterior quantiles there; here they are the ends of the Wald
#'   interval at those probabilities.
#' @param pars Row names to keep, in the order given, as in brms. A name
#'   the fit does not have is an error listing the ones it does.
#' @param flatten If `TRUE`, one named vector of ESTIMATES in the
#'   `confint()` spelling instead of the summary matrix. Coefficients
#'   come in linear-predictor order, which need not be the matrix's row
#'   order; index by name. Distributional parameters held at a constant
#'   are included here and are absent from `vcov()` and from the matrix,
#'   which cover the ESTIMATED population-level coefficients only.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return A coefficients-by-four matrix, or with `flatten = TRUE` a
#'   single named vector of estimates.
#' @seealso [vcov.frmtmb_fit()], which names its rows as this matrix
#'   does, and [confint.frmtmb_fit()], which names them as
#'   `flatten = TRUE` names its entries.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#'
#' # brms's summary matrix, one row per population-level coefficient
#' fit <- frm(bf(y ~ x + (1 | g), sigma ~ x) + gaussian(), data = dd)
#' fixef(fit)
#' exp(fixef(fit)["sigma_Intercept", "Estimate"])  # sigma is on the log
#'
#' # flatten to the vector confint() names its rows by
#' fixef(fit, flatten = TRUE)
#' all(names(fixef(fit, flatten = TRUE)) %in% rownames(confint(fit)))
#'
#' # the standard error column is vcov()'s diagonal
#' all.equal(fixef(fit)[, "Est.Error"], sqrt(diag(vcov(fit))))
#' @rdname fixef
#' @aliases fixef
#' @export
fixef.frmtmb_fit <- function(object, summary = TRUE, robust = FALSE,
                             probs = c(0.025, 0.975), pars = NULL, ...,
                             flatten = FALSE) {
  frm_check_dots(...)
  # brms's positions ahead of `...`: fixef(fit, FALSE) used to set
  # `flatten` and return the default shape, identical() to fixef(fit)
  fit_refuse_draws_args("fixef()", summary = summary, robust = robust)
  require_fitted(object, "fixef()")
  check_flag(flatten, "flatten")
  est <- object$estimates
  # the estimates already CARRY the canonical names (the parameter
  # template builds both), so flatten reads them rather than rebuilding
  # the dpar-prefix rule a second place where it could drift
  if (flatten) {
    if (!is.null(pars)) {
      frm_stop("fixef(flatten = TRUE) names its entries in the internal ",
               "vocabulary, where brms's `pars` names rows of the summary ",
               "matrix. Subset the vector, or drop flatten", call. = FALSE)
    }
    out <- numeric(0)
    for (lp in object$frame[["linpreds"]]) {
      out <- c(out, est[[lp[["par"]]]][lp[["idx"]]])
    }
    return(out)
  }
  rows <- brms_fixef_rows(object)
  cf <- unname(brms_fixef_values(object, rows))
  V <- tryCatch(suppressWarnings(vcov(object)), error = function(e) NULL)
  se <- if (is.null(V)) NULL else unname(sqrt(diag(V))[rows$names])
  out <- brms_summary_matrix(cf, se, probs, rownames = rows$names)
  brms_pars_filter(out, pars, "fixef()")
}

#' The estimated fixed-effect coefficients in `brms_coef_table()` order,
#' which is `estimated_coef_names()` order.
#'
#' @noRd
fixef_estimated <- function(object) {
  est <- object$estimates
  bd <- est[["betad"]]
  if (length(fx <- object$frame[["betad_fixed_idx"]])) bd <- bd[-fx]
  unname(c(est[["beta"]], bd))
}

#' Fixed effects per linear predictor
#'
#' The coefficients of each linear predictor, as a named list keyed by
#' distributional parameter (and by response on a multivariate fit),
#' each entry named by DESIGN COLUMN: `fixef_by_dpar(fit)$sigma`,
#' `fixef_by_dpar(fit)$mu[["(Intercept)"]]`.
#'
#' This is the shape [fixef()] returned before frmtmb took brms's
#' summary matrix. It is kept because nothing else has it: `fixef()`
#' names its rows as brms does and flattens every predictor into one
#' table, and `fixef(flatten = TRUE)` is one vector in the [confint()]
#' spelling. A model with several linear predictors (a mixture, a
#' multivariate response, a nonlinear body) is the case this reads
#' cleanly.
#'
#' It carries ESTIMATES only. For a standard error or an interval use
#' [fixef()], whose rows are the same coefficients under brms's names.
#'
#' @param object A `frmtmb_fit`.
#' @return A named list of coefficient vectors, one per linear
#'   predictor.
#' @seealso [fixef()], [coef.frmtmb_fit()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x, exp(0.2 + 0.1 * dd$x))
#' fit <- frm(bf(y ~ x, sigma ~ x) + gaussian(), data = dd)
#'
#' # one entry per distributional parameter, named by design column
#' fixef_by_dpar(fit)
#' exp(fixef_by_dpar(fit)$sigma[["(Intercept)"]])
#'
#' # the same coefficients under brms's names, with their errors
#' fixef(fit)
#' @export
fixef_by_dpar <- function(object) {
  require_fitted(object, "fixef_by_dpar()")
  est <- object$estimates
  out <- list()
  for (lp in object$frame[["linpreds"]]) {
    v <- est[[lp[["par"]]]][lp[["idx"]]]
    names(v) <- colnames(lp[["X"]])
    out[[coef_block_key(object, lp)]] <- v
  }
  out
}

#' Extract random-effect modes
#' @param object A `frmtmb_fit`.
#' @param summary,robust,probs,pars,groups brms's arguments, in brms's
#'   positions so that a positional brms call asks the same question.
#'   brms answers `summary = FALSE` with the posterior draws and
#'   `robust = TRUE` with their median and MAD, and a maximum-likelihood
#'   fit has no draws, so both are refused by name with the reason. The
#'   default of each is accepted and changes nothing.
#' @param condVar If `TRUE`, attach the conditional SDs of the modes
#'   (from the Laplace posterior) as a `"condSD"` attribute on each
#'   matrix, in matching layout.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return A named list of levels-by-coefficients matrices, one per
#'   random-effect term, KEYED BY THE GROUPING FACTOR as brms and lme4
#'   key it (so `ranef(fit)$g` and `coef(fit)$g` name the same group).
#'   Each matrix carries its block label in a `"term"` attribute, which
#'   is also the `grp` column of `as.data.frame()`. That long form
#'   (with a `condsd` column when
#'   `condVar = TRUE` was used) is what broom.mixed-style code reads.
#'
#'   Several terms on ONE grouping factor give several entries under one
#'   name, and `$` and `[[` then take the block label instead:
#'   `ranef(fit)[["chi: 1 | id"]]`. The bare factor name is refused
#'   there, and names the labels, rather than returning the first block.
#'   Positional indexing reaches every block either way.
#'
#'   Three ordinary spellings put more than one block on one factor, so
#'   this is not a nonlinear-model corner: an uncorrelated slope
#'   `(1 + x || g)`, which desugars to `1 | g` and `0 + x | g`; two bars
#'   written out on the same factor, `(1 | g) + (0 + x | g)`, the same
#'   two blocks; and a random effect on more than one distributional or
#'   nonlinear parameter, such as `bf(y ~ x + (1 | g), sigma ~ (1 | g))`
#'   (`1 | g` and `sigma: 1 | g`) or `(1 | id)` on three nonlinear
#'   parameters. A CORRELATED slope `(1 + x | g)` is one block and is
#'   addressed by the factor name as before.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # one matrix per random-effect term, levels by coefficients
#' ranef(fit)
#' ranef(fit)$g[1:3, ]                 # keyed by the grouping factor
#' attr(ranef(fit)$g, "term")          # the block it came from
#' ranef(fit)[["1 | g"]][1:3, ]        # or by that label
#'
#' # condVar adds the conditional SDs a caterpillar plot needs
#' re <- as.data.frame(ranef(fit, condVar = TRUE))
#' head(re)
#' with(re[order(re$condval), ],
#'      plot(condval, seq_along(condval), pch = 16,
#'           xlim = range(condval - 2 * condsd, condval + 2 * condsd),
#'           xlab = "conditional mode", ylab = "group"))
#' @rdname ranef
#' @aliases ranef
#' @export
ranef.frmtmb_fit <- function(object, summary = TRUE, robust = FALSE,
                             probs = c(0.025, 0.975), pars = NULL,
                             groups = NULL, ..., condVar = FALSE) {
  frm_check_dots(...)
  # brms's positions ahead of `...`: ranef(fit, FALSE) used to set
  # `condVar` and return the default shape, identical() to ranef(fit)
  fit_refuse_draws_args("ranef()", summary = summary, robust = robust,
                        probs = probs, pars = pars, groups = groups)
  require_fitted(object, "ranef()")
  check_flag(condVar, "condVar")
  cvec <- coef_b(object)
  # variances, not SDs: an esicar block's field is the CENTERED b, so
  # what it reports is the PROJECTED variance, and the projection
  # subtracts a variance
  cvr <- NULL
  if (condVar) {
    sdr <- sdr_of(object)
    dcr <- sdr$diag.cov.random
    if (!is.null(dcr)) {
      cvr <- pmax(dcr[names(sdr$par.random) == "b"], 0)
    }
  }
  out <- list()
  for (bk in object$frame[["re_blocks"]]) {
    M <- t(matrix(cvec[bk[["c_idx"]]], nrow = bk[["dim"]]))
    # brms's coefficient names, which is what coef() puts on the same
    # columns: `Intercept`, not `(Intercept)`, and `sigma_Intercept` for
    # a block on sigma. A block brms has no r_ for (a reduced-rank
    # factor, a smooth basis) keeps the names it carries.
    rn <- bk[["cnms"]]
    if (brms_block_has_r(bk)) {
      bn <- tryCatch(brms_re_rnames(object, bk), error = function(e) NULL)
      if (length(bn) == length(rn)) rn <- bn
    }
    dimnames(M) <- list(bk[["levels"]], rn)
    if (!is.null(cvr)) {
      # rr factors live in a different space than the displayed
      # coefficients; no conditional SDs for those blocks. An esicar
      # block's coefficients are a different space too, but a knowable
      # one: ranef() shows the centered field, whose variance is var(b)
      # with the inert component mean projected out.
      S <- if (bk[["covstruct"]] == "rr") {
        matrix(NA_real_, nrow(M), ncol(M))
      } else if (block_is_esicar(bk)) {
        t(matrix(car_center_condsd(cvr[bk[["b_idx"]]], bk[["aux_car"]]),
                 nrow = bk[["dim"]]))
      } else {
        t(matrix(sqrt(cvr[bk[["b_idx"]]]), nrow = bk[["dim"]]))
      }
      dimnames(S) <- dimnames(M)
      attr(M, "condSD") <- S
    }
    # the block's own label stays reachable: it is what VarCorr() keys
    # by and the only thing that tells two blocks on one factor apart
    attr(M, "term") <- bk[["term_label"]]
    # appended, then named: `out[[label]] <- M` would DROP a block whose
    # key repeats ((1 | g) in mu and in sigma both key the factor g)
    out[[length(out) + 1L]] <- M
  }
  # keyed by the GROUPING FACTOR, as brms and lme4 key it, and as this
  # package's own coef() already did: ranef(fit)$Subject used to be NULL
  # in a model where coef(fit)$Subject was a data frame
  names(out) <- vapply(object$frame[["re_blocks"]], function(bk) {
    bk[["group_name"]] %||% bk[["term_label"]]
  }, "")
  structure(out, class = "ranef_frmtmb")
}

#' Blocks that share a grouping factor stay addressable.
#'
#' 0.52.0 re-keyed this list by the GROUPING FACTOR, which is brms's and
#' lme4's key and the one `coef()` already used. That is the right key
#' and it does not change here. What it left behind is a list whose
#' names repeat: three `(1 | id)` terms on three nonlinear parameters
#' are three entries all called `id`, and `[["id"]]` reached the first
#' of them silently - a wrong answer rather than an error, and no
#' spelling reached the other two at all.
#'
#' So the KEY grows a second accepted form rather than the list changing
#' shape: the block label, which is what `attr(blk, "term")` carries,
#' and what `as.data.frame()` puts in `grp`.
#' A grouping factor with one block is addressed as before; a factor
#' with several is either addressed by its block label or refused by
#' name. Positional indexing is untouched, so `print()`,
#' `as.data.frame()` and frmtmb.sample's `ranef.frmtmb_draws()` - all of
#' which index by position on purpose - go through the fast path.
#'
#' @noRd
ranef_pick <- function(x, i) {
  y <- unclass(x)
  hit <- which(names(y) == i)
  if (length(hit) == 1L) return(list(found = TRUE, value = y[[hit]]))
  tl <- vapply(y, function(m) attr(m, "term") %||% NA_character_, "")
  ti <- which(tl == i)
  if (length(ti) == 1L) return(list(found = TRUE, value = y[[ti]]))
  if (length(hit) > 1L) {
    frm_stop("ranef() has ", length(hit), " random-effect blocks on ",
             "grouping factor '", i, "', so '", i, "' does not name one of ",
             "them. Address a block by its term label - ",
             paste0("[[\"", tl[hit], "\"]]", collapse = ", "),
             " - or by position. The label is also each block's \"term\" ",
             "attribute and the `grp` column of as.data.frame()",
             call. = FALSE)
  }
  list(found = FALSE, value = NULL)
}

#' @export
`[[.ranef_frmtmb` <- function(x, i, ...) {
  if (!is.character(i) || length(i) != 1L) return(unclass(x)[[i, ...]])
  got <- ranef_pick(x, i)
  if (got$found) return(got$value)
  # whatever base does with a name this list does not carry, which on a
  # list is NULL rather than an error
  unclass(x)[[i, ...]]
}

#' @export
`$.ranef_frmtmb` <- function(x, name) {
  got <- ranef_pick(x, name)
  if (got$found) return(got$value)
  y <- unclass(x)
  # base `$` on a list partial-matches and returns NULL for a miss;
  # keep both, or a name that used to work would silently become NULL
  j <- pmatch(name, names(y))
  if (is.na(j)) NULL else y[[j]]
}

#' @export
print.ranef_frmtmb <- function(x, ...) {
  frm_check_dots(...)
  for (i in seq_along(x)) {           # by position: names can repeat
    tl <- attr(x[[i]], "term")
    cat("$", names(x)[i], if (!is.null(tl)) paste0("   (", tl, ")"),
        "\n", sep = "")
    print(`attr<-`(`attr<-`(x[[i]], "condSD", NULL), "term", NULL))
    cat("\n")
  }
  invisible(x)
}

#' @export
as.data.frame.ranef_frmtmb <- function(x, ...) {
  frm_check_dots(...)
  rows <- lapply(seq_along(x), function(i) {   # by position: see print()
    M <- x[[i]]
    # `grp` names the BLOCK, which is what tells two terms on one factor
    # apart; the list itself is keyed by the factor, brms's and lme4's
    # key, so the two are not the same string any more
    nm <- attr(M, "term") %||% names(x)[i]
    S <- attr(M, "condSD")
    lv <- rownames(M) %||% as.character(seq_len(nrow(M)))
    df <- data.frame(
      grp = nm,
      term = rep(colnames(M), each = nrow(M)),
      level = rep(lv, ncol(M)),
      condval = as.vector(M)
    )
    if (!is.null(S)) df$condsd <- as.vector(S)
    df
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @export
as.data.frame.VarCorr_frmtmb <- function(x, ...) {
  frm_check_dots(...)
  rows <- list()
  # by position, not by name: two blocks can share a term label (an
  # animal model's (1 | gr(id, cov = A)) and its permanent-environment
  # (1 | id) both deparse to "1 | id"), and x[[nm]] would then return
  # the first block once per duplicate name
  for (i in seq_along(x)) {
    nm <- names(x)[i]
    V <- x[[i]]
    sds <- sqrt(diag(V))
    cn <- colnames(V)
    for (i in seq_along(sds)) {
      rows[[length(rows) + 1L]] <- data.frame(
        grp = nm, var1 = cn[i], var2 = NA_character_,
        vcov = V[i, i], sdcor = sds[i]
      )
    }
    if (ncol(V) > 1L) {
      C <- stats::cov2cor(V)
      for (i in seq_len(ncol(V) - 1L)) {
        for (j in seq(i + 1L, ncol(V))) {
          rows[[length(rows) + 1L]] <- data.frame(
            grp = nm, var1 = cn[i], var2 = cn[j],
            vcov = V[i, j], sdcor = C[i, j]
          )
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' The covariance matrix of every random-effect block, one per block in
#' block order and named by term label, at a given `theta`.
#'
#' This is what `VarCorr()` returned before it took brms's shape, and it
#' is still what the printed fit, `summary()` and the singularity check
#' read: a covariance matrix per BLOCK is the natural unit for those,
#' while brms's `VarCorr()` is keyed by GROUP. A smooth, GP, CAR or SPDE
#' block contributes its one variance under its own label, as before.
#'
#' @noRd
varcorr_matrices <- function(fit, theta = fit$estimates[["theta"]]) {
  th <- theta
  out <- lapply(fit$frame[["re_blocks"]], function(bk) {
    if (bk[["covstruct"]] == "smooth") {
      # one smoothing variance; the k x k identity blowup is noise
      matrix(exp(th[bk[["theta_idx"]]])^2, 1, 1,
             dimnames = list("sd(wiggle)", "sd(wiggle)"))
    } else if (bk[["covstruct"]] %in% c("gp", "hsgp")) {
      # marginal GP sd; the range lives in confint_varcorr
      matrix(exp(th[bk[["theta_idx"]][1]])^2, 1, 1,
             dimnames = list("sd(gp)", "sd(gp)"))
    } else {
      V <- covstruct_registry[[bk[["covstruct"]]]]$vcov(th[bk[["theta_idx"]]],
        bk)
      if (is_student_block(bk)) {
        # On a gr(dist = "student") block this is the SCALE matrix, not
        # the covariance: brms names the same quantity `sd_<group>__...`
        # and frmtmb keeps that name, so the matrix is tagged instead of
        # silently converted. The variance is `nu/(nu-2)` times it.
        attr(V, "dist_nu") <- bk[["dist_nu"]]
      }
      V
    }
  })
  names(out) <- vapply(fit$frame[["re_blocks"]], `[[`, "", "term_label")
  structure(out, class = "VarCorr_frmtmb")
}

#' The layout of brms's `VarCorr()` for a fit: one entry per grouping
#' factor, in first-appearance order, then `residual__`.
#'
#' Each group entry is keyed by brms's group name (`brms_group_name()`,
#' `g:h2` for an interaction) and lists its coefficient names (brms's
#' `rnames`, `brms_re_rnames()`), which block each came from and at which
#' position, and whether any block in it carries correlations. Blocks on
#' one factor merge, as brms merges every term on a group; a coefficient
#' given twice on one group is refused when the model is assembled, as
#' brms refuses it. Smooth, GP, CAR and SPDE blocks are not in brms's
#' `VarCorr()` and are not here; [confint_varcorr()] reports them.
#'
#' `residual__` follows `brms:::VarCorr.brmsfit()`. On a univariate
#' model it is there when `sigma` is a parameter nobody wrote a formula
#' for, one row named by brms's `bterms$resp`, which is `""`. On a
#' multivariate model it is there when at least one response has such a
#' `sigma` and no response predicts it: one row per such response, named
#' by the response, with the residual correlations when `rescor = TRUE`.
#'
#' @noRd
varcorr_layout <- function(fit) {
  groups <- list()
  for (bi in seq_along(fit$frame[["re_blocks"]])) {
    bk <- fit$frame[["re_blocks"]][[bi]]
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp", "car", "spde")) {
      next
    }
    rn <- brms_re_rnames(fit, bk)
    key <- brms_group_name(bk)
    g <- groups[[key]] %||% list(rnames = character(0), block = integer(0),
                                 pos = integer(0), cor = FALSE)
    g$rnames <- c(g$rnames, rn)
    g$block <- c(g$block, rep(bi, length(rn)))
    g$pos <- c(g$pos, seq_along(rn))
    g$cor <- g$cor || (bk[["dim"]] > 1L &&
                         !bk[["covstruct"]] %in% c("diag", "homdiag"))
    groups[[key]] <- g
  }
  list(groups = groups, residual = varcorr_residual_layout(fit))
}

#' brms's `residual__` entry of `VarCorr()`, see `varcorr_layout()`:
#' the coefficient positions of the natural-scale `sigma` rows, their
#' inverse links, the row names, and whether residual correlations
#' follow. `NULL` when brms has no such entry.
#'
#' @noRd
varcorr_residual_layout <- function(fit) {
  tab <- brms_coef_table(fit)
  is_sigma <- tab$dpar %in% "sigma"
  simple <- which(is_sigma & tab$natural)
  if (!length(simple)) return(NULL)
  mv <- length(fit$spec$responses) > 1L
  if (!mv) {
    return(list(rnames = "", pos = simple[1L],
                linkinv = attr(tab, "linkinv")[simple[1L]], rescor = FALSE))
  }
  # brms: any simple sigma and no predicted one
  if (any(is_sigma & !tab$natural)) return(NULL)
  list(rnames = brms_stan_name(tab$resp[simple]), pos = simple,
       linkinv = attr(tab, "linkinv")[simple],
       rescor = isTRUE(fit$spec$rescor) &&
         length(simple) == length(fit$spec$responses))
}

#' The numbers behind brms's `VarCorr()` at one parameter vector, in
#' the `hyp_vals_only()` layout: per entry of `varcorr_layout()`, the
#' standard deviations and, when the group has correlations, the full
#' correlation and covariance matrices. A pair of coefficients from two
#' different blocks is uncorrelated by construction, which is the 0
#' brms fills in for a correlation it does not find.
#'
#' @noRd
varcorr_values <- function(fit, vals, comp, layout = varcorr_layout(fit)) {
  th <- vals[comp == "theta"]
  mats <- if (length(layout$groups)) varcorr_matrices(fit, th) else list()
  out <- list()
  for (key in names(layout$groups)) {
    g <- layout$groups[[key]]
    K <- length(g$rnames)
    S <- matrix(0, K, K)
    for (i in seq_len(K)) {
      for (j in seq_len(K)) {
        if (g$block[i] == g$block[j]) {
          S[i, j] <- mats[[g$block[i]]][g$pos[i], g$pos[j]]
        }
      }
    }
    sd <- sqrt(diag(S))
    e <- list(sd = sd)
    if (g$cor) {
      C <- S / tcrossprod(sd)
      diag(C) <- 1
      e$cor <- C
      e$cov <- S
    }
    out[[key]] <- e
  }
  if (!is.null(r <- layout$residual)) {
    cf <- c(vals[comp == "beta"], vals[comp == "betad"])
    sd <- vapply(seq_along(r$pos), function(i) {
      r$linkinv[[i]](cf[r$pos[i]])
    }, 1)
    e <- list(sd = sd)
    if (r$rescor) {
      C <- us_chol_cor(vals[comp == "thetar"], length(sd))
      # exactly 1, so the diagonal's error is exactly 0 as in a group
      diag(C) <- 1
      e$cor <- C
      e$cov <- C * tcrossprod(sd)
    }
    out[["residual__"]] <- e
  }
  out
}

#' Extract random-effect standard deviations and correlations
#'
#' brms's `VarCorr()`, with frequentist content. The return value is a
#' list with one entry per GROUPING FACTOR, named by the factor as brms
#' and lme4 name it (`"patient"`, not the term label `"1 | patient"`),
#' followed by `residual__` when `sigma` is a parameter nobody wrote a
#' formula for (one row per response on a multivariate model, with the
#' residual correlations under `rescor = TRUE`). Each entry has:
#'
#' - `sd`: a matrix with one row per coefficient and the columns
#'   `Estimate`, `Est.Error` and one quantile column per `probs`
#'   (`Q2.5`, `Q97.5`);
#' - `cor` and `cov`, when the group has correlations: arrays of
#'   `coefficient x statistic x coefficient`, the same statistics.
#'
#' The names are brms's: the entry is the group as brms spells it
#' (`g:h` for an interaction), and the coefficients are `Intercept`, `x`,
#' and for a distributional or nonlinear parameter `sigma_Intercept`, as
#' in `sd_<group>__sigma_Intercept`.
#'
#' @section What the columns mean on a maximum-likelihood fit:
#' brms summarizes posterior draws. A fit has one estimate and its
#' sampling distribution, so the same columns carry:
#'
#' - `Estimate`: the maximum-likelihood (or REML) estimate.
#' - `Est.Error`: its delta-method standard error, from the joint
#'   covariance of the covariance parameters ([vcov()] with
#'   `full = TRUE`).
#' - `Q<p>`: the Wald quantile `Estimate + qnorm(p) * Est.Error`, on the
#'   natural scale. That is the interval [hypothesis()] reports for the
#'   same quantity. It can cross zero for a standard deviation near its
#'   boundary; [confint_varcorr()] with `method = "profile"` gives an
#'   interval that respects the boundary.
#'
#' A diagonal entry of `cor` is 1 with error 0.
#'
#' `summary = FALSE` and `robust = TRUE` are refused by name: brms
#' returns the draws, or their median and MAD, and a fit has no draws.
#' `frmtmb.sample::frm_sample()` gives an object whose `VarCorr()` does
#' both.
#'
#' Several random-effect blocks on one grouping factor merge into its
#' entry, as brms merges every term on a group, and a pair across two
#' blocks has correlation 0. Two terms that give one group the same
#' coefficient, as in `(1 | gr(id, cov = A)) + (1 | id)`, are refused
#' when the model is built, as brms refuses them: give the second term a
#' copy of the grouping column under another name.
#' Smooth, Gaussian-process, CAR and SPDE blocks are not in brms's
#' `VarCorr()` and are not here: [confint_varcorr()] reports them. On a
#' `gr(dist = "student")` block the `sd` row is the SCALE, as brms's
#' `sd_` is.
#'
#' @param x A `frmtmb_fit`.
#' @param sigma Ignored, as in brms. It is carried by nlme's generic,
#'   which frmtmb shares.
#' @param summary Must be `TRUE`; see above.
#' @param robust Must be `FALSE`; see above.
#' @param probs The quantiles to report.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return A named list, as above.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
#' u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.4))
#' dd$y <- rnorm(200, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
#' fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)
#'
#' vc <- VarCorr(fit)
#' names(vc)                       # "g" and "residual__", as in brms
#' vc$g$sd                         # standard deviations
#' vc$g$cor["Intercept", "Estimate", "x"]
#' vc$residual__$sd
#' @rdname VarCorr
#' @aliases VarCorr
#' @export
VarCorr.frmtmb_fit <- function(x, sigma = 1, summary = TRUE,
                               robust = FALSE,
                               probs = c(0.025, 0.975), ...) {
  frm_check_dots(...)
  require_fitted(x, "VarCorr()")
  fit_refuse_draws_args("VarCorr()", summary = summary, robust = robust)
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) ||
        any(probs < 0) || any(probs > 1)) {
    frm_stop("VarCorr(): `probs` must be numbers between 0 and 1",
             call. = FALSE)
  }
  lay <- varcorr_layout(x)
  if (!length(lay$groups) && is.null(lay$residual)) {
    frm_stop("The model does not contain covariance matrices.", call. = FALSE)
  }
  pc <- hyp_par_cov(x)
  flat <- function(v) {
    unlist(lapply(varcorr_values(x, v, pc$comp, lay), function(e) {
      c(e$sd, e$cor, e$cov)
    }), use.names = FALSE)
  }
  q0 <- flat(pc$vals)
  rel <- which(pc$comp %in% c("theta", "betad", "thetar"))
  J <- matrix(0, length(q0), length(rel))
  for (k in seq_along(rel)) {
    i <- rel[k]
    step <- max(1e-5, 1e-5 * abs(pc$vals[i]))
    vp <- pc$vals; vp[i] <- vp[i] + step
    vm <- pc$vals; vm[i] <- vm[i] - step
    J[, k] <- (flat(vp) - flat(vm)) / (2 * step)
  }
  Vr <- pc$V[rel, rel, drop = FALSE]
  se <- sqrt(pmax(0, rowSums((J %*% Vr) * J)))
  stats_nm <- c("Estimate", "Est.Error", paste0("Q", probs * 100))
  tab <- cbind(q0, se, outer(se, stats::qnorm(probs)) + q0)
  colnames(tab) <- stats_nm
  vals <- varcorr_values(x, pc$vals, pc$comp, lay)
  pos <- 0L
  take <- function(n) {
    r <- tab[pos + seq_len(n), , drop = FALSE]
    pos <<- pos + n
    r
  }
  rn_of <- function(key) {
    if (identical(key, "residual__")) lay$residual$rnames else
      lay$groups[[key]]$rnames
  }
  out <- list()
  for (key in names(vals)) {
    rn <- rn_of(key)
    K <- length(rn)
    sdm <- take(K)
    rownames(sdm) <- rn
    e <- list(sd = sdm)
    if (!is.null(vals[[key]]$cor)) {
      for (el in c("cor", "cov")) {
        blk <- take(K * K)
        # the flattened matrix is column-major, so row i of `blk` is
        # entry (row = (i - 1) %% K + 1, column = (i - 1) %/% K + 1)
        arr <- array(NA_real_, c(K, ncol(tab), K),
                     dimnames = list(rn, stats_nm, rn))
        for (s in seq_len(ncol(tab))) {
          arr[, s, ] <- matrix(blk[, s], K, K)
        }
        e[[el]] <- arr
      }
    }
    out[[key]] <- e
  }
  out
}

#' The refusal a maximum-likelihood method owes brms's draws arguments,
#' by name and with the reason. brms answers `summary = FALSE` with the
#' posterior draws and `robust = TRUE` with their median and MAD; a fit
#' has one estimate and no chains, so returning the summary anyway
#' would answer a different question under the same call.
#'
#' @noRd
fit_refuse_draws_args <- function(what, summary = TRUE, robust = FALSE,
                                  probs = NULL, pars = NULL,
                                  groups = NULL) {
  check_flag(summary, "summary")
  check_flag(robust, "robust")
  tail <- paste0(" Sample with frmtmb.sample::frm_sample() and call ",
                 what, " on the draws for that")
  if (!summary) {
    frm_stop(what, " cannot honor summary = FALSE: brms returns the ",
             "posterior draws there, and a maximum-likelihood fit has none. ",
             "It carries one estimate and no chains, so the only answer it ",
             "has is the summary, which is a different return shape.", tail,
             call. = FALSE)
  }
  if (robust) {
    frm_stop(what, " cannot honor robust = TRUE: brms's robust summary is ",
             "the median and MAD of the draws, and a maximum-likelihood fit ",
             "has no draws.", tail, call. = FALSE)
  }
  if (!is.null(probs) && !isTRUE(all.equal(probs, c(0.025, 0.975)))) {
    frm_stop(what, " cannot honor `probs`: brms reports quantiles of the ",
             "draws, and this method returns point estimates with no ",
             "quantile columns. For an interval here use confint().", tail,
             call. = FALSE)
  }
  if (!is.null(pars)) {
    frm_stop(what, " cannot honor `pars`: ",
             brms_draws_summary_args[["pars"]], call. = FALSE)
  }
  if (!is.null(groups)) {
    frm_stop(what, " cannot honor `groups`: brms's `groups` selects ",
             "grouping factors; the return value here is a named list, so ",
             "index it", call. = FALSE)
  }
  invisible(NULL)
}

#' @export
print.VarCorr_frmtmb <- function(x, ...) {
  # by position: duplicate term labels are legal (see
  # as.data.frame.VarCorr_frmtmb), and name lookup would print the
  # first block once per duplicate and never print the others
  frm_check_dots(...)
  for (i in seq_along(x)) {
    nm <- names(x)[i]
    V <- x[[i]]
    sdv <- sqrt(diag(V))
    nu <- attr(V, "dist_nu")
    cat(" ", nm, "\n")
    # a t block's diagonal is the SCALE, so the column is not headed
    # Std.Dev.: sd = scale * sqrt(nu/(nu-2)), printed alongside
    tab <- if (is.null(nu)) {
      data.frame(Name = colnames(V), `Std.Dev.` = signif(sdv, 5),
                 check.names = FALSE)
    } else {
      data.frame(Name = colnames(V), Scale = signif(sdv, 5),
                 `Std.Dev.` = signif(sdv * sqrt(nu / (nu - 2)), 5),
                 check.names = FALSE)
    }
    if (ncol(V) > 1) {
      C <- stats::cov2cor(V)
      corr <- format(signif(C, 3))
      corr[upper.tri(corr, diag = TRUE)] <- ""
      tab <- cbind(tab, Corr = corr[, -ncol(corr), drop = FALSE])
    }
    print(tab, row.names = FALSE)
    if (!is.null(nu)) {
      cat("   Student-t latent, nu = ", format(nu),
          " (fixed); the stored matrix is the scale\n", sep = "")
    }
  }
  invisible(x)
}

# ---- brms methods with nothing to do here ----------------------------

#' Expose a model's compiled functions
#'
#' brms compiles Stan functions and `expose_functions()` makes them
#' callable from R. frmtmb compiles nothing: the model is an R closure
#' built by `build_objective()` and differentiated by RTMB, and a
#' custom family's density is the plain R function handed to
#' [custom_family()]. The method is defined so that a ported brms
#' script gets that reason rather than "could not find function".
#'
#' @param x A `frmtmb_fit`.
#' @param ... Ignored; this method always stops.
#' @return This function never returns; it signals an error.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(40))
#' dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
#' fit <- frm(bf(y ~ x) + gaussian(), data = dd)
#' try(expose_functions(fit))
#' @export
expose_functions <- function(x, ...) UseMethod("expose_functions")

#' @rdname expose_functions
#' @exportS3Method brms::expose_functions
#' @export
expose_functions.frmtmb_fit <- function(x, ...) {
  frm_stop("expose_functions() has no Stan program to read on a frmtmb ",
           "fit: a custom family's lpdf is the plain R function handed to ",
           "custom_family(), callable as it is", call. = FALSE)
}
