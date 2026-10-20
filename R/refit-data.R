# The data a refit reads, kept on the fit.
#
# drop1(), update() without data, influence() and frmtmb.latent's
# hmm_starts() refit the model by evaluating its formula again. The model
# frame alone cannot serve them: it holds `log(x2)`, not `x2`. Reading
# the data the call names from the workspace instead gives a different
# model, or none, once that object is changed or removed, and always
# after saveRDS() and readRDS() in a new session.
# dev/drop1raw-findings.md has the measurements behind the choice.

#' The data a refit of `fit` reads: the model frame's columns and the
#' raw variables the formula reads through a transform
#' (`frame$raw_vars`), on the rows after `na.action`. brms keeps the
#' same thing as `fit$data`. `rows = "fit"` gives the rows the fit
#' used; `rows = "data"` adds back the rows a univariate `subset()`
#' cut, as brms's stored data holds them, so a refit re-applies the
#' cut and `nobs()` counts them. `NULL` when the fit could not keep the
#' raw variables (data given as an environment).
#'
#' @noRd
fit_refit_data <- function(fit, rows = c("fit", "data")) {
  rows <- match.arg(rows)
  fr <- fit[["frame"]]
  mf <- fr[["data_frame"]]
  raw <- fr[["raw_vars"]]
  if (is.null(mf) || is.null(raw)) return(NULL)
  cols <- as.list(mf)
  for (v in setdiff(names(raw), names(cols))) cols[[v]] <- raw[[v]]
  rn <- attr(mf, "row.names")
  cut <- fr[["subset_cut"]]
  if (rows == "data" && !is.null(cut)) {
    keep <- cut[["keep"]]
    ord <- order(c(keep, setdiff(seq_len(cut[["n"]]), keep)))
    both <- function(a, b) {
      if (is.null(dim(a))) c(a, b)[ord] else rbind(a, b)[ord, , drop = FALSE]
    }
    cols <- stats::setNames(lapply(names(cols), function(v) {
      both(cols[[v]], cut[["data"]][[v]])
    }), names(cols))
    rn <- c(rownames(mf), rownames(cut[["data"]]))[ord]
  }
  # a plain frame: the model frame's `terms` and the `data_name` that
  # frm() records would otherwise travel into the refit
  structure(cols, class = "data.frame", row.names = rn)
}

#' The variables a formula, a [bf()] or an [mvbf()] reads, over every
#' formula it holds.
#'
#' @noRd
formula_vars <- function(x) {
  if (is.language(x)) return(all.vars(x))
  if (is.list(x) && !is.data.frame(x)) {
    return(unique(unlist(lapply(unclass(x), formula_vars),
                         use.names = FALSE)))
  }
  character(0)
}

#' `data` with the channel `frm()` reads from it, the "frmtmb_refit"
#' attribute: `call` is the call the refit records, which names the
#' data as the original call did; `predvar_map` is the fit's frozen
#' bases; `na_action` is the rows the fit's `na.action` removed. A
#' refit on the fit's own rows (`fit` given) takes all three, so it
#' reproduces the fit: computed again on those rows, a spline's knots
#' or `scale()`'s center would differ when `na.action` removed rows. A
#' refit on the call's data takes the call alone. `formula`, `family`
#' and the row count name the call the channel was made for
#' (`refit_channel_for()`).
#'
#' @noRd
refit_channel <- function(data, formula, family = NULL, call = NULL,
                          fit = NULL) {
  attr(data, "frmtmb_refit") <- list(
    formula = formula, family = family, n = NROW(data), call = call,
    predvar_map = if (!is.null(fit)) fit[["frame"]][["predvar_map"]],
    na_action = if (!is.null(fit)) fit[["frame"]][["na_action"]])
  data
}

#' The channel `data` carries, when it was made for this `frm()` call:
#' the same formula and family values and the same number of rows.
#' A call edited after `update(evaluate = FALSE)` returned it, or its
#' data given to another `frm()` call, is another model; the call it
#' would record names the old one, and the old fit's bases and NA rows
#' do not belong to it. `NULL` then, and `frm()` records its own call.
#'
#' @noRd
refit_channel_for <- function(data, formula, family) {
  ch <- attr(data, "frmtmb_refit")
  if (is.null(ch)) return(NULL)
  if (!identical(formula, ch[["formula"]]) ||
        !identical(family, ch[["family"]]) ||
        !identical(NROW(data), ch[["n"]])) {
    return(NULL)
  }
  ch
}

#' Does the objective a fit maximized hold a prior density? Read from the
#' resolved prior entries: a prior that only sets bounds adds nothing
#' to the objective, and `logLik()` is then the log-likelihood.
#'
#' @noRd
fit_prior_density <- function(fit) {
  pr <- fit[["prior"]]
  if (is.null(pr) || !length(pr)) return(FALSE)
  ri <- tryCatch(resolve_prior_input(list(frame = fit[["frame"]],
                                          spec = fit[["spec"]]), pr),
                 error = function(e) NULL)
  # a prior that cannot be resolved again cannot be shown to be flat
  is.null(ri) || length(ri$entries) > 0L
}

#' Refuse a likelihood-ratio statistic between fits whose objective
#' holds a prior density. Such a fit maximizes the log-likelihood plus
#' the log prior, so a difference of two maxima has a term of prior
#' densities at two different estimates: it is not a likelihood ratio
#' and it can be negative.
#'
#' @noRd
prior_lrt_refusal <- function(what, which) {
  frm_stop(what, " compares log-likelihoods, and ", which, " a prior ",
           "that adds its log density to the objective (frm(prior = )). ",
           "Its logLik() is the log-likelihood plus the log prior, so a ",
           "difference of two is not a likelihood ratio and has no ",
           "chi-square reference. For a likelihood-ratio test, refit ",
           "without the prior; a prior that only sets bounds can stay. ",
           "To compare the models under the prior, use loo() or ",
           "bayes_factor() of frmtmb.sample", call. = FALSE)
}

#' Does `data`, the object the stored call names now, still hold the
#' data the fit was fitted to? The fit's rows (and the rows its
#' `na.action` removed) are its rows, and every variable the fit holds
#' has the fit's values on them. `update()` reads a variable the fit
#' does not hold from such data only; from changed data every other
#' variable would be read again too, and the refit would be another
#' model.
#'
#' @noRd
refit_data_matches <- function(fit, data) {
  rd <- fit_refit_data(fit, "data")
  if (is.null(rd) || is.null(data) || is.environment(data) ||
        !is.list(data) || !length(data)) {
    return(FALSE)
  }
  n <- if (is.data.frame(data)) nrow(data) else NROW(data[[1L]])
  rn <- if (is.data.frame(data)) rownames(data) else {
    as.character(seq_len(n))
  }
  dropped <- unclass(fit[["frame"]][["na_action"]])
  if (length(dropped) && (any(dropped < 1L) || any(dropped > n))) {
    return(FALSE)
  }
  if (nrow(rd) + length(dropped) != n ||
        !setequal(c(rownames(rd), rn[dropped]), rn)) {
    return(FALSE)
  }
  # every variable the fit read from its data by name: one missing here
  # would be read from wherever R finds it, and silently
  if (!all(fit[["frame"]][["data_vars"]] %in% names(data))) {
    return(FALSE)
  }
  rows <- match(rownames(rd), rn)
  norm <- function(x) {
    if (is.factor(x)) return(as.character(x))
    d <- dim(x)
    attributes(x) <- NULL
    if (!is.null(d)) dim(x) <- d
    x
  }
  for (v in intersect(names(rd), names(data))) {
    x <- data[[v]]
    x <- if (is.null(dim(x))) x[rows] else x[rows, , drop = FALSE]
    if (!identical(norm(x), norm(rd[[v]]))) return(FALSE)
  }
  TRUE
}

#' The data argument of a call as it is spelled, for a message: a data
#' frame passed by value has no spelling worth printing.
#'
#' @noRd
data_label <- function(x) {
  if (is.name(x) || is.call(x)) substr(deparse1(x), 1L, 50L) else "a value"
}
