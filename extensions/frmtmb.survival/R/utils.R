# Argument checks for this package's constructors. frmtmb has the same
# checks and does not export them, and a constructor's argument check is
# not a seam worth a public promise, so these are the package's own,
# with core's wording so that a refused argument reads the same in both.

`%||%` <- function(x, y) if (is.null(x)) y else x

#' A short description of a refused value: its class and, for a scalar,
#' the value itself.
#'
#' @noRd
sv_arg_desc <- function(x) {
  if (is.null(x)) return("NULL")
  cls <- class(x)[1L]
  if (is.atomic(x) && length(x) == 1L) {
    if (is.na(x)) return(paste0(cls, " NA"))
    return(paste0(cls, " ", encodeString(format(x), quote = "\"")))
  }
  art <- if (substr(cls, 1L, 1L) %in% c("a", "e", "i", "o", "u")) "an" else
    "a"
  paste0(art, " ", cls, " of length ", length(x))
}

#' @noRd
sv_check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    frm_stop("`", arg, "` must be TRUE or FALSE, not ", sv_arg_desc(x),
             call. = FALSE)
  }
  invisible(x)
}

#' @noRd
sv_check_count <- function(x, arg, min = 0L) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x == round(x) && x >= min
  if (!ok) {
    frm_stop("`", arg, "` must be a single whole number of at least ", min,
             ", not ", sv_arg_desc(x), call. = FALSE)
  }
  invisible(as.integer(x))
}

#' The link a constructor was handed, as brms and `stats::family()` read
#' it: an UNQUOTED link name is taken as its name, so `cox(identity)`
#' does not evaluate base R's identity function.
#'
#' @noRd
sv_link_value <- function(slink, link, choices, default) {
  if (is.character(slink) && length(slink) == 1L && !is.na(slink)) {
    return(slink)
  }
  text <- if (is.name(slink)) as.character(slink) else deparse1(slink)
  if (text %in% choices) return(text)
  if (is.name(slink)) link <- tryCatch(link, error = function(e) text)
  if (is.null(link) || identical(link, NA)) return(default)
  if (is.character(link) || is.list(link)) return(link)
  text
}

#' A link name checked against the set brms allows for the family, or a
#' custom link object passed through to `frmtmb_family()`.
#'
#' @noRd
sv_check_link <- function(value, choices, family) {
  if (is.list(value)) return(value)
  allowed <- paste0("'", choices, "'", collapse = ", ")
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    frm_stop(family, "(link =) takes a single link name, not ",
             sv_arg_desc(value), ". Supported links are: ", allowed,
             ". See ?frmtmb::`frmtmb-links`", call. = FALSE)
  }
  if (!value %in% choices) {
    frm_stop("'", value, "' is not a supported link for family '", family,
             "'. Supported links are: ", allowed,
             ". See ?frmtmb::`frmtmb-links`", call. = FALSE)
  }
  value
}
