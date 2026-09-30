#' Build a conditions data frame for conditional_effects()
#'
#' Makes one row per combination of values of `vars`, which is the
#' data frame [conditional_effects()] takes as `conditions` to draw one
#' panel per row. This is brms's `make_conditions()`, and its output is
#' identical to brms's on the same input.
#'
#' A factor, character or logical variable contributes each of its
#' levels. A numeric variable contributes three values: its mean minus
#' one standard deviation, its mean, and its mean plus one standard
#' deviation. The last variable of `vars` varies fastest.
#'
#' @param x A data frame, or a fitted model. From a `frmtmb_fit` the
#'   model data is used; from any other list that has a `data` element,
#'   that element is used, as brms does.
#' @param vars Character vector of the names of the variables to
#'   condition on.
#' @param ... Passed to the labeler that makes the `cond__` column:
#'   `digits` (default 2) rounds a numeric value in the label, `sep`
#'   (default `" & "`) separates the variables, and `incl_vars`
#'   (default `TRUE`) puts the variable name before each value.
#' @return A data frame with one column per variable of `vars`, in the
#'   order of `vars`, and a column `cond__` that labels each row, for
#'   example `"zBase = -1 & zAge = 1"`. As in brms, `cond__` is the
#'   variable itself when there is one variable and `incl_vars = FALSE`.
#' @seealso [conditional_effects()], [conditional_smooths()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(60), z = rnorm(60),
#'                  f = factor(rep(c("a", "b"), 30)))
#' dd$y <- rnorm(60, dd$x + dd$z * (dd$f == "b"))
#' fit <- frm(bf(y ~ x + z * f), family = gaussian(), data = dd)
#' conds <- make_conditions(fit, c("z", "f"))
#' conds
#' ce <- conditional_effects(fit, effects = "x", conditions = conds)
#' @export
make_conditions <- function(x, vars, ...) {
  vars <- rev(as.character(vars))
  if (inherits(x, "frmtmb_fit")) {
    x <- x$frame[["data_frame"]]
  } else if (!is.data.frame(x) && "data" %in% names(x)) {
    x <- x[["data"]]
  }
  x <- as.data.frame(x)
  miss <- setdiff(vars, names(x))
  if (length(miss)) {
    frm_stop("make_conditions(): ", paste0("`", miss, "`", collapse = ", "),
             " is not a column of the data. It has: ",
             paste(names(x), collapse = ", "), call. = FALSE)
  }
  out <- stats::setNames(vector("list", length(vars)), vars)
  for (v in vars) {
    tmp <- x[[v]]
    out[[v]] <- if (brms_like_factor(tmp)) {
      levels(as.factor(tmp))
    } else {
      mean(tmp, na.rm = TRUE) + (-1:1) * stats::sd(tmp, na.rm = TRUE)
    }
  }
  out <- rev(expand.grid(out))
  out$cond__ <- brms_rows2labels(out, ...)
  out
}

#' brms's `is_like_factor()`: what a conditions labeler rounds and what
#' it prints as a level. Logical and character count as factors there.
#'
#' @noRd
brms_like_factor <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
}

#' One label per row, as brms's `rows2labels()` makes them, because the
#' labels become panel titles and a ported script compares them as
#' strings.
#'
#' @noRd
brms_rows2labels <- function(x, digits = 2, sep = " & ", incl_vars = TRUE,
                             ...) {
  x <- as.data.frame(x)
  check_flag(incl_vars, "incl_vars")
  out <- x
  for (i in seq_along(out)) {
    if (!brms_like_factor(out[[i]])) {
      out[[i]] <- round(out[[i]], digits)
    }
    if (incl_vars) {
      out[[i]] <- paste0(names(out)[i], " = ", out[[i]])
    }
  }
  Reduce(function(a, b) paste(a, b, sep = sep), out)
}

#' Update the addition terms of a formula
#'
#' Replaces or adds the addition terms of a model formula, the terms
#' after `|` on the left-hand side, such as `trials()`, `weights()`,
#' `se()` or `cens()`. This is brms's `update_adterms()`, and for a
#' plain formula its output is identical to brms's on the same input.
#'
#' With `action = "update"` an addition term of `adform` replaces the
#' term of the same name in `formula`, and the other terms of `formula`
#' stay. With `action = "replace"` the addition terms of `formula` are
#' all removed first, so only those of `adform` remain. A term name
#' spelled with brms's legacy `resp_` prefix counts as the name without
#' it.
#'
#' @param formula A two-sided formula, or a formula made by [bf()]. For
#'   a [bf()] formula the response formula is updated and every other
#'   part (distributional formulas, constants, `nl`, the family) stays.
#' @param adform A one-sided formula of addition terms, such as
#'   `~ trials(10)` or `~ weights(w) + se(sei)`.
#' @param action `"update"` (default) or `"replace"`.
#' @return An object of the class of `formula`: a formula with the same
#'   attributes and environment, or the updated [bf()] formula.
#' @seealso [bf()]
#' @examples
#' update_adterms(y | trials(size) ~ x, ~ trials(10))
#' update_adterms(y | trials(size) ~ x, ~ weights(w))
#' update_adterms(y | trials(size) ~ x, ~ weights(w), action = "replace")
#' update_adterms(bf(y ~ x, sigma ~ x), ~ weights(w))
#' @export
update_adterms <- function(formula, adform, action = c("update", "replace")) {
  action <- frm_match_arg(action)
  if (inherits(formula, "frmtmb_mvformula")) {
    frm_stop("update_adterms() updates one response formula. Call it on ",
             "each response's bf() before combining them with mvbf()",
             call. = FALSE)
  }
  if (inherits(formula, "frmtmb_formula")) {
    formula[["formula"]] <- update_adterms(formula[["formula"]], adform,
                                           action)
    return(formula)
  }
  if (!inherits(formula, "formula")) {
    frm_stop("update_adterms(): `formula` must be a formula or a bf() ",
             "formula, not ", arg_desc(formula), call. = FALSE)
  }
  if (!inherits(adform, "formula")) {
    frm_stop("update_adterms(): `adform` must be a one-sided formula of ",
             "addition terms, such as ~ trials(10), not ", arg_desc(adform),
             call. = FALSE)
  }
  if (length(formula) != 3L) {
    frm_stop("update_adterms() needs a formula with a response: ",
             "addition terms are on the left-hand side, and ",
             deparse1(formula), " has none", call. = FALSE)
  }
  # the string route is brms's own, kept so that the output is the one a
  # brms script compares against, term order and spacing included
  str_formula <- gsub("[ \t\r\n]+", " ",
                      gsub("[\t\r\n]+", " ",
                           Reduce(paste, deparse(formula)), perl = TRUE),
                      perl = TRUE)
  old_ad <- regmatches(str_formula,
                       gregexpr("(?<=\\|)[^~]*(?=~)", str_formula,
                                perl = TRUE))[[1L]]
  new_ad_terms <- attr(stats::terms(adform), "term.labels")
  if (action == "update" && length(old_ad)) {
    old_ad <- stats::formula(paste("~", old_ad))
    old_ad_terms <- attr(stats::terms(old_ad), "term.labels")
    ad_name <- function(tl) {
      sub("^resp_", "", unlist(regmatches(tl, gregexpr("^[^\\(]+", tl))))
    }
    keep <- !ad_name(old_ad_terms) %in% ad_name(new_ad_terms)
    new_ad_terms <- c(old_ad_terms[keep], new_ad_terms)
  }
  if (length(new_ad_terms)) {
    new_ad_terms <- paste("|", paste(new_ad_terms, collapse = "+"))
  }
  resp <- gsub("\\|.+", "", paste(deparse(formula[[2L]]), collapse = ""))
  out <- stats::formula(paste(resp, paste(new_ad_terms, collapse = ""),
                              "~1"))
  out[[3L]] <- formula[[3L]]
  attributes(out) <- attributes(formula)
  out
}
