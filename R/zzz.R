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
