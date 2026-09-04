# Argument plumbing, cloned rather than imported.
#
# These three are frmtmb internals that hmm.R and lca.R reached while
# they lived inside it. They are not extension API and should not
# become it: a general TRUE/FALSE validator on frmtmb's public surface
# would be there forever to serve one call here. Copies keep the
# message a user sees identical to the character, and message
# uniqueness is a per-package property, so a template repeated from
# frmtmb resolves unambiguously to this tree.
#
# See dev/out-of-tree-inventory.md for the whole resolution table.

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' How a rejected argument value is described back to the user.
#'
#' @noRd
arg_desc <- function(x) {
  if (is.null(x)) return("NULL")
  cls <- class(x)[1L]
  if (is.atomic(x) && length(x) == 1L) {
    if (is.na(x)) return(paste0(cls, " NA"))
    return(paste0(cls, " ", encodeString(format(x), quote = "\"")))
  }
  art <- if (substr(cls, 1L, 1L) %in% c("a", "e", "i", "o", "u")) "an" else "a"
  paste0(art, " ", cls, " of length ", length(x))
}

#' The one refusal every TRUE/FALSE argument in the package shares.
#'
#' @noRd
check_flag <- function(x, arg, what = NULL) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE, not ", arg_desc(x),
         if (is.null(what)) "" else paste0(". ", what), call. = FALSE)
  }
  invisible(x)
}
