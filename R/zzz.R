.onLoad <- function(libname, pkgname) {
  if (requireNamespace("emmeans", quietly = TRUE)) {
    emmeans::.emm_register("frmtmb_fit", pkgname)
  }
  invisible()
}
