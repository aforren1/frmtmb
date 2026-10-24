# Fix round 1: a wrapper is found at any depth, and re_formula is checked
# the same way (A7), on the tree named as the argument.
root <- commandArgs(trailingOnly = TRUE)[1]
rep_file <- function(p, old, new) {
  p <- file.path(root, p)
  s <- paste(readLines(p, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, s, fixed = TRUE)[[1]])
  if (regexpr(old, s, fixed = TRUE) < 0 || n != 1) {
    stop(p, ": count ", n, ": ", substr(old, 1, 70))
  }
  writeLines(sub(old, new, s, fixed = TRUE), p)
}
rep_file("R/parse.R", 'check_bar_wrapper <- function(tm, known) {
  if (!is.call(tm) || !is.name(tm[[1L]])) return(invisible(NULL))
  hd <- as.character(tm[[1L]])
  if (hd %in% c("(", "|", "||", ":", "*", "/", "-", known)) {
    return(invisible(NULL))
  }
  is_bar <- function(e) {
    while (is.call(e) && identical(e[[1L]], as.name("("))) e <- e[[2L]]
    is.call(e) && is.name(e[[1L]]) &&
      as.character(e[[1L]]) %in% c("|", "||")
  }
  if (!any(vapply(as.list(tm)[-1L], is_bar, TRUE))) return(invisible(NULL))
  if (hd %in% c("us_t", "diag_t")) {', 'check_bar_wrapper <- function(tm, known, where = NULL) {
  hd <- bar_wrapper_head(tm, known)
  if (is.null(hd)) return(invisible(NULL))
  if (!is.null(where)) {
    frm_stop(where, ": ", hd, "() is not a covariance structure, so ",
             deparse1(tm), " names no group-level term of the model. Name ",
             "a term as the model writes it, (1 | g) or diag(x | g)",
             call. = FALSE)
  }
  if (hd %in% c("us_t", "diag_t")) {')
rep_file("R/parse.R", '#\' Is `tm` a call to a registered field special, as a whole term?', '#\' The head of the first call between the top of `tm` and a bar term
#\' that names no structure, or `NULL`. The path to every bar is walked
#\' through `(` and through the known wrappers, at any depth, so
#\' `foo(bar(1 | g))` and `diag(foo(1 | g))` are found as `foo(1 | g)` is;
#\' a bar\'s own arguments (gr(), mm(), mo()) are not on the path. The
#\' crossing operators keep their own refusal in parse_linpred().
#\'
#\' @noRd
bar_wrapper_head <- function(tm, known) {
  is_bar <- function(e) {
    is.call(e) && is.name(e[[1L]]) && as.character(e[[1L]]) %in% c("|", "||")
  }
  has_bar <- function(e) {
    if (!is.call(e)) return(FALSE)
    if (is_bar(e)) return(TRUE)
    any(vapply(as.list(e)[-1L], function(a) {
      !identical(a, quote(expr = )) && has_bar(a)
    }, NA))
  }
  rec <- function(e) {
    if (!is.call(e) || is_bar(e)) return(NULL)
    hd <- if (is.name(e[[1L]])) as.character(e[[1L]]) else deparse1(e[[1L]])
    if (hd %in% c(":", "*", "/", "-")) return(NULL)
    if (hd %in% c("(", "+", known)) {
      for (a in as.list(e)[-1L]) {
        if (identical(a, quote(expr = ))) next
        h <- rec(a)
        if (!is.null(h)) return(h)
      }
      return(NULL)
    }
    if (has_bar(e)) hd else NULL
  }
  rec(tm)
}

#\' Is `tm` a call to a registered field special, as a whole term?')
rep_file("R/re-formula.R", '    } else if (op %in% c(cs_names(), "rr") &&
                 length(e) >= 2L && is.call(e[[2L]]) &&
                 deparse1(e[[2L]][[1L]]) %in% c("|", "||")) {
      # a covariance wrapper, diag(x | g), which brms does not have:
      # the term it wraps names the columns, as the fit\'s own term does
      walk(e[[2L]])
    }', '    } else if (op %in% c(cs_names(), "rr")) {
      # a covariance wrapper, diag(x | g), which brms does not have:
      # the term it wraps names the columns, as the fit\'s own term does
      for (a in as.list(e)[-1L]) {
        if (is.call(a)) walk(a)
      }
    } else {
      # any other call around a bar term names no term of the model; it
      # was ignored, and the prediction kept a term it did not name
      check_bar_wrapper(e, c(cs_names(), "rr"), where = "re_formula")
    }')
cat("fix-a7 applied to", root, "\n")
