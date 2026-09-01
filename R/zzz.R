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
