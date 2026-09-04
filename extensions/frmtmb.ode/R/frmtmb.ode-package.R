#' @keywords internal
"_PACKAGE"

# frmtmb is a Depends, so that one library(frmtmb.ode) call gives a user
# the formula grammar frm_ode() is written inside. The two seams this
# package builds on are imported by name as well, because a namespace
# that is loaded and not attached reaches nothing through the search
# path.
#' @importFrom frmtmb frmtmb_register_frame_check frmtmb_ad_overload
NULL
