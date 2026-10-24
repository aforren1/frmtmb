#' @keywords internal
#'
#' @section How this package reaches frmtmb:
#' Through exported functions only. The structures are registered with
#' [frmtmb::frmtmb_register_covstruct()] and their rows in
#' [frmtmb::frm_compat()] with [frmtmb::frmtmb_register_compat()], both
#' from `.onLoad()`. The package reads no frmtmb internal: no
#' `frmtmb:::`, no `getFromNamespace()`, no `asNamespace("frmtmb")`.
#'
#' Without this package, a formula that writes `car()` or `spde()` is
#' refused by name, and the refusal says to load it.
"_PACKAGE"

# frmtmb is a Depends, so that one library(frmtmb.spatial) call gives a
# user the formula grammar and the frm() the field is fitted with. The
# seams this package builds on are imported by name as well, because a
# namespace that is loaded and not attached reaches nothing through the
# search path.
#' @importFrom frmtmb compat_rule_builder frm_message frm_stop
#'   frmtmb_register_compat frmtmb_register_covstruct
#' @importFrom RTMB dgmrf dmvnorm
#' @importFrom stats rnorm
NULL
