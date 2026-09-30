# brms's row-selection addition terms, subset() and index(), and the
# mi(x, idx = ) predictor that uses them.
#
# In brms each response of a multivariate model is fitted on its own
# rows: `y | subset(s) ~ x` keeps the rows where `s` is TRUE, so two
# responses can have different numbers of observations. A predictor
# that reads another response's values, mi(x), then cannot line the two
# up by row number, and brms asks for the rows to be matched by value:
# `x | index(id)` names each of x's rows, and `mi(x, idx = ref)` picks,
# for each row of the using response, the row of x whose `id` equals
# `ref`.
#
# frmtmb follows that design. One model frame holds every row, which is
# what keeps the NA handling brms's; each response then takes its own
# rows of it before anything is built, so its design matrices, smooth
# bases, grouping levels and addition-term values are computed on those
# rows alone, as brms computes them on the subsetted data.

#' Whether any response of a model carries `subset()`.
#'
#' @noRd
spec_has_subset <- function(spec) {
  any(vapply(spec$responses, function(r) !is.null(r$aterms[["subset"]]),
             TRUE))
}

#' Rows of a model frame, keeping the `terms` attribute that the design
#' builders read.
#'
#' @noRd
frame_rows <- function(mf, rows) {
  if (is.null(rows)) return(mf)
  out <- mf[rows, , drop = FALSE]
  attr(out, "terms") <- attr(mf, "terms")
  attr(out, "na.action") <- attr(mf, "na.action")
  out
}

#' An addition term's value on the model frame: the frame column when
#' the term is a variable of it, else the expression evaluated there.
#'
#' @noRd
aterm_raw_value <- function(resp, nm, mf) {
  ex <- resp$aterms[[nm]]
  v <- mf[[deparse1(ex)]]
  if (is.null(v)) v <- eval(ex, mf, resp$formula_env)
  v
}

#' One response's `subset()` on a frame, with brms's checks.
#'
#' @noRd
subset_eval <- function(resp, mf) {
  v <- aterm_raw_value(resp, "subset", mf)
  lab <- paste0("subset(", deparse1(resp$aterms[["subset"]]), ")")
  if (length(v) != nrow(mf)) {
    frm_stop(lab, " of response '", resp$resp_name, "': the length of ",
             "'subset' does not match the rows of 'data'. It takes a ",
             "logical variable with one value per row", call. = FALSE)
  }
  # as.logical() reads "TRUE" and "FALSE", and brms refuses a character
  # subset variable ("invalid argument type"), so frmtmb does too
  if (is.character(v)) {
    frm_stop(lab, " of response '", resp$resp_name, "': a character ",
             "subset variable is refused, as brms refuses it. Give a ",
             "logical variable, or a number or factor that as.logical() ",
             "reads", call. = FALSE)
  }
  v <- as.logical(v)
  if (anyNA(v)) {
    frm_stop(lab, " of response '", resp$resp_name, "': subset variables ",
             "may not contain NAs, or values that are not TRUE or FALSE",
             call. = FALSE)
  }
  v
}

#' The variables of the combined model frame that one response uses,
#' named the way `model.frame()` names its columns.
#'
#' @noRd
frame_rhs_columns <- function(rhs) {
  if (is.null(rhs)) return(character(0))
  vars <- as.list(attr(stats::terms(stats::as.formula(call("~", rhs))),
                       "variables"))[-1L]
  vapply(vars, function(x) {
    paste(deparse(x, width.cutoff = 500L,
                  backtick = !is.symbol(x) && is.language(x)),
          collapse = " ")
  }, "")
}

#' The rows to drop from a multivariate frame with `subset()`, brms's
#' `na_omit()`: an NA drops its row unless every response that uses the
#' variable leaves that row out through its own `subset()`. `exempt`
#' are the columns whose NAs are read rather than dropped (`mi()`
#' responses).
#'
#' @noRd
subset_na_rows <- function(spec, mf, resp_cols, exempt) {
  subs <- lapply(spec$responses, function(r) {
    if (!is.null(r$aterms[["subset"]])) subset_eval(r, mf)
  })
  bad <- rep(FALSE, nrow(mf))
  for (cn in setdiff(names(mf), exempt)) {
    v <- mf[[cn]]
    isna <- if (is.matrix(v)) rowSums(is.na(v)) > 0 else is.na(v)
    if (!any(isna)) next
    keep <- rep(TRUE, nrow(mf))
    used <- FALSE
    for (r in names(spec$responses)) {
      if (!cn %in% resp_cols[[r]]) next
      used <- TRUE
      keep <- keep & (if (is.null(subs[[r]])) FALSE else !subs[[r]])
    }
    if (!used) keep <- FALSE
    bad <- bad | (isna & !keep)
  }
  bad
}

#' The rows of the model frame each response uses: a named list with an
#' integer vector for every response that carries `subset()`.
#'
#' @noRd
subset_rows_of <- function(spec, mf) {
  out <- list()
  for (r in spec$responses) {
    if (is.null(r$aterms[["subset"]])) next
    rows <- which(subset_eval(r, mf))
    if (!length(rows)) {
      frm_stop("subset(", deparse1(r$aterms[["subset"]]), ") of response '",
               r$resp_name, "' leaves it no rows. Make sure the variables ",
               "have no NAs on the rows the response is meant to use, and ",
               "that the subset variable is TRUE on at least one of them",
               call. = FALSE)
    }
    out[[r$resp_name]] <- rows
  }
  out
}

#' Models `subset()` cannot join, refused before the frame is built.
#'
#' @noRd
subset_check_spec <- function(spec, na.action) {
  if (!spec_has_subset(spec)) return(invisible(NULL))
  if (identical(na.action, stats::na.exclude)) {
    frm_stop("subset() cannot be combined with na.action = na.exclude: ",
             "a response fitted on some rows has no place to pad the rows ",
             "it leaves out. Use na.omit (the default)", call. = FALSE)
  }
  if (isTRUE(spec$rescor)) {
    frm_stop("subset() cannot be combined with rescor = TRUE: the residual ",
             "correlation pairs the responses row by row, so every ",
             "response needs the same rows. Use set_rescor(FALSE)",
             call. = FALSE)
  }
  has_me <- any(vapply(spec$responses, function(r) {
    any(vapply(r$dpars, function(dp) length(dp[["meterms"]] %||% list()) > 0L,
               TRUE))
  }, TRUE))
  if (has_me) {
    frm_stop("Argument 'subset' is not supported when using 'me' terms, as ",
             "in brms: a noise-free value is shared by every response on ",
             "the same row", call. = FALSE)
  }
  if (length(spec$responses) > 1L) {
    for (r in spec$responses) {
      if (is.null(r$aterms[["subset"]])) next
      if (!is.null(r$autocor)) {
        frm_stop("subset() on response '", r$resp_name, "' cannot be ",
                 "combined with a residual correlation term in a ",
                 "multivariate model yet: the time and grouping structure ",
                 "has not been checked on a response's own rows",
                 call. = FALSE)
      }
      if (!is.null(fam_structure(r$family))) {
        frm_stop("subset() on response '", r$resp_name, "' cannot be ",
                 "combined with the structured family '",
                 r$family[["family"]], "' in a multivariate model: its ",
                 "sequences and groups are read over the whole frame",
                 call. = FALSE)
      }
    }
  }
  invisible(NULL)
}

#' One response's `index()` on its own rows, with brms's duplicate
#' check. The raw values are kept: `mi(x, idx = )` matches by value, as
#' brms's `match()` does.
#'
#' @noRd
index_eval <- function(resp, mf) {
  v <- aterm_raw_value(resp, "index", mf)
  if (length(v) == 1L && nrow(mf) > 1L) {
    frm_stop("index(", deparse1(resp$aterms[["index"]]), ") of response '",
             resp$resp_name, "' is one value; it takes a variable that ",
             "names each row", call. = FALSE)
  }
  if (anyDuplicated(v)) {
    frm_stop("Index of response '", resp$resp_name, "' contains duplicated ",
             "values. index(", deparse1(resp$aterms[["index"]]), ") names ",
             "each row of the response once", call. = FALSE)
  }
  if (is.factor(v)) v <- as.character(v)
  v
}

#' For a predictor `mi(x)` of response `resp`, the row of `x` that each
#' of `resp`'s rows reads: `NULL` when the two share their rows, else
#' brms's `idxl`, `match(idx, index(x))`. Refused as brms refuses where
#' the rows cannot be lined up.
#'
#' @noRd
mi_idx_rows <- function(ent, vn, resp, tgt, mf, index_vals, sub_rows) {
  if (is.null(ent$idx)) {
    if (!is.null(sub_rows[[resp$resp_name]])) {
      frm_stop("mi() terms in subsetted formulas require the 'idx' ",
               "argument to be specified: response '", resp$resp_name,
               "' uses subset(), so its rows are not the rows of '", vn,
               "'. Write mi(", vn, ", idx = <variable>) and give '", vn,
               "' an index() term", call. = FALSE)
    }
    if (!is.null(sub_rows[[vn]])) {
      frm_stop("mi() terms of subsetted variables require the 'idx' ",
               "argument to be specified: response '", vn, "' uses ",
               "subset(), so its rows are not the rows of '",
               resp$resp_name, "'. Write mi(", vn,
               ", idx = <variable>) and give '", vn, "' an index() term",
               call. = FALSE)
    }
    return(NULL)
  }
  if (is.null(tgt$aterms[["index"]])) {
    frm_stop("Response '", vn, "' needs to have an 'index' addition term ",
             "to compare with 'idx': write bf(", vn, " | mi() + ",
             "index(<variable>) ~ ...). See ?mi for examples",
             call. = FALSE)
  }
  iv <- eval(ent$idx, mf, resp$formula_env)
  if (is.factor(iv)) iv <- as.character(iv)
  m <- match(iv, index_vals[[vn]])
  if (anyNA(m)) {
    frm_stop("Could not match all indices in response '", vn, "': ",
             sum(is.na(m)), " value(s) of ", deparse1(ent$idx), " on the ",
             "rows of '", resp$resp_name, "' are not among the values of ",
             "index(", deparse1(tgt$aterms[["index"]]), ") on the rows of '",
             vn, "'", call. = FALSE)
  }
  m
}

#' brms's rule for the post-fit methods of a multivariate model with
#' `subset()`: the responses have different rows, so they are asked for
#' one at a time.
#'
#' @noRd
subset_resp_check <- function(object, resp, what) {
  if (!length(object$frame[["subset_rows"]])) return(invisible(NULL))
  if (length(resp) != 1L) {
    frm_stop(what, ": argument 'resp' must be a single variable name for ",
             "models using addition argument 'subset'. The responses are ",
             "fitted on different rows, so their answers do not stack; ",
             "ask for one of ",
             paste0("'", names(object$spec$responses), "'", collapse = ", "),
             call. = FALSE)
  }
  invisible(NULL)
}

#' The rows of `newdata` a response with `subset()` predicts: those
#' where its subset is TRUE, as brms keeps them. `newdata` must hold
#' the subset variable, as brms requires.
#'
#' @noRd
subset_newdata <- function(object, resp, newdata) {
  if (is.null(newdata) || length(resp) != 1L) return(newdata)
  rspec <- object$spec$responses[[resp]]
  ex <- rspec$aterms[["subset"]]
  if (is.null(ex)) return(newdata)
  v <- tryCatch(eval(ex, newdata, rspec$formula_env),
                error = function(e) NULL)
  if (is.null(v) || length(v) != nrow(newdata)) {
    miss <- setdiff(all.vars(ex), names(newdata))
    frm_stop("subset(", deparse1(ex), ") of response '", resp, "' could ",
             "not be evaluated on newdata",
             if (length(miss)) {
               paste0(": newdata has no column ",
                      paste(miss, collapse = ", "))
             },
             ". brms predicts a subsetted response on the rows of newdata ",
             "where the subset is TRUE, so newdata needs the variable",
             call. = FALSE)
  }
  # the fit's rule: brms refuses a character subset variable
  if (is.character(v)) {
    frm_stop("subset(", deparse1(ex), ") of response '", resp, "' on ",
             "newdata: a character subset variable is refused, as brms ",
             "refuses it. Give a logical variable, or a number or factor ",
             "that as.logical() reads", call. = FALSE)
  }
  v <- as.logical(v)
  if (anyNA(v)) {
    frm_stop("subset(", deparse1(ex), ") of response '", resp, "' is NA ",
             "or not TRUE or FALSE on some rows of newdata", call. = FALSE)
  }
  newdata[v, , drop = FALSE]
}

#' The values of `x` a predictor `mi(x, idx = )` reads on `newdata`:
#' for each row, `x` on the row of `newdata` whose `index()` equals the
#' row's `idx`, among the rows `x`'s own `subset()` keeps there.
#'
#' @noRd
mi_idx_newdata <- function(fit, mt, newdata, env) {
  tgt <- fit$spec$responses[[mt$var]]
  xd <- subset_newdata(fit, mt$var, newdata)
  iv <- tryCatch(eval(tgt$aterms[["index"]], xd, tgt$formula_env),
                 error = function(e) NULL)
  ref <- tryCatch(eval(mt$idx_expr, newdata, env), error = function(e) NULL)
  if (is.null(iv) || is.null(ref)) {
    frm_stop("mi(", mt$var, ", idx = ", deparse1(mt$idx_expr), "): ",
             "newdata needs the variables of idx and of index(",
             deparse1(tgt$aterms[["index"]]), ")", call. = FALSE)
  }
  if (is.factor(iv)) iv <- as.character(iv)
  if (is.factor(ref)) ref <- as.character(ref)
  m <- match(ref, iv)
  if (anyNA(m)) {
    frm_stop("Could not match all indices in response '", mt$var, "' ",
             "on newdata: ", sum(is.na(m)), " value(s) of ",
             deparse1(mt$idx_expr), " are not among the values of index(",
             deparse1(tgt$aterms[["index"]]), ") on the rows of newdata ",
             "that '", mt$var, "' uses", call. = FALSE)
  }
  xd[[mt$var]][m]
}
