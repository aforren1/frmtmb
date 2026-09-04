#' Load hook.
#'
#' The one thing this package must tell frmtmb about itself: a frame
#' check that refuses a dynamics input which is not constant inside a
#' solve group. Registering at load time is what keeps frmtmb free of any
#' mention of ODEs.
#'
#' @noRd
.onLoad <- function(libname, pkgname) {
  frmtmb::frmtmb_register_frame_check(check_ode_constancy)
  invisible()
}
