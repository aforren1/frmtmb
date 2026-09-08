#' Load hook for the optional downstream integrations. It registers the
#' `frmtmb_fit` class with emmeans when that suggested package is
#' installed, and adds the class to the `marginaleffects` model
#' whitelist, because marginaleffects only dispatches to methods for
#' classes on that list.
#'
#' @noRd
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("emmeans", quietly = TRUE)) {
    emmeans::.emm_register("frmtmb_fit", pkgname)
  }
  # marginaleffects gates on a class whitelist; the methods themselves
  # are registered via delayed S3 registration
  cls <- getOption("marginaleffects_model_classes", NULL)
  if (!"frmtmb_fit" %in% cls) {
    options(marginaleffects_model_classes = c(cls, "frmtmb_fit"))
  }
  invisible()
}

#' Session state for the notices frmtmb gives once.
#'
#' A convention this package applies silently, and which a user porting
#' a model from another package would not expect, is announced ONCE per
#' session. A warning on every fit is noise by the tenth model, and no
#' signal at all is how a ported model changes its answer quietly. The
#' state is a package-level environment, so it resets with the session
#' and never reaches a saved fit.
#'
#' @noRd
frmtmb_notice_state <- new.env(parent = emptyenv())

#' Give a notice the first time its `tag` comes up in a session.
#'
#' Suppressible two ways: `suppressMessages()`, because it is a
#' message and nothing else, and `options(frmtmb.notices = FALSE)` for
#' a session that wants none of them. The tag is marked as given before
#' the message is raised, so a suppressed notice is spent rather than
#' saved for the next fit. The package's own test suite sets the option
#' off in `setup.R`; the tests that exercise a notice turn it back on.
#'
#' @noRd
notify_once <- function(tag, ...) {
  if (!isTRUE(getOption("frmtmb.notices", TRUE))) return(invisible(FALSE))
  if (exists(tag, envir = frmtmb_notice_state, inherits = FALSE)) {
    return(invisible(FALSE))
  }
  assign(tag, TRUE, envir = frmtmb_notice_state)
  message(...)
  invisible(TRUE)
}
