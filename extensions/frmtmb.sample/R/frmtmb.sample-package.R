#' @keywords internal
#' @import frmtmb
#' @importFrom stats coef family formula getCall nobs predict update
#' @importFrom graphics pairs
"_PACKAGE"

# The stats and graphics generics above are imported by name because
# this package registers `frmtmb_draws` METHODS on them. A method is
# registered against the generic the namespace can see, so without the
# import the namespace fails to load with "object 'family' not found",
# which is what R CMD check reports rather than anything about the
# method itself.

# frmtmb's whole export surface is imported rather than qualified call
# by call. The two packages are developed and released together from one
# repository and frmtmb.sample is written against a version bound
# (`Imports: frmtmb (>= 0.46.0)`), so the ambiguity `@import` is usually
# warned about - not knowing which package a name came from - does not
# arise here: every unqualified name in this package is either defined
# in it or is frmtmb's.
#
# The internals this package needs are frmtmb EXPORTS as of 0.46.0, not
# `:::` reaches. dev/draws-extraction.md is the inventory that made them
# exports, one documented contract each, and `?frmtmb::"frmtmb-sampling-api"`
# is the contract as shipped. Nothing here reaches into frmtmb's
# namespace past its exports; that was the point of taking the inventory
# before moving the code.

# ---- helpers this package keeps its own copy of ----------------------
#
# Context-free utilities that frmtmb also has. They are duplicated
# rather than exported from there on purpose: each is a few lines that
# know nothing about frmtmb, and making three argument-check messages
# into stable cross-package API would cost more than the duplication
# does. See the (c) column of dev/draws-extraction.md.

#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' A single number that is a whole count, at least `min`.
#'
#' @noRd
check_count <- function(x, arg, min = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != as.integer(x) ||
        x < min) {
    stop(arg, " must be a single whole number", if (min > 0L) {
      paste0(" of at least ", min)
    }, "; got ", paste(format(x), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' A single number strictly inside (0, 1).
#'
#' @noRd
check_probability <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    stop(arg, " must be a single number strictly between 0 and 1; got ",
         paste(format(x), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' A named list, with an example of the spelling wanted.
#'
#' @noRd
check_named_list <- function(x, arg, example) {
  ok <- is.list(x) && (!length(x) ||
                         (!is.null(names(x)) && all(nzchar(names(x)))))
  if (!ok) {
    stop(arg, " must be a named list, as in ", example, call. = FALSE)
  }
  invisible(TRUE)
}

#' Pointwise percentile of a draws matrix (draws in rows), NA where a
#' column has no usable draw rather than an error out of `quantile()`.
#'
#' @noRd
ce_pctl <- function(m, p) {
  apply(m, 2, function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) NA_real_ else unname(stats::quantile(v, p))
  })
}
