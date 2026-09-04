#' @keywords internal
"_PACKAGE"

# frmtmb is a Depends, so that one library(frmtmb.ddm) call gives a user
# the formula grammar and the frm() the family is passed to. The compat
# seam is imported by name as well, because a namespace that is loaded
# and not attached reaches nothing through the search path.
#' @importFrom frmtmb custom_family frmtmb_register_compat
#'   compat_rule_builder
NULL
