# What this package tells frmtmb about itself at load time.
#
# Only the compatibility matrix. A pair rule belongs to whichever
# package makes the pair possible, and every rule naming hmm() or lca()
# against a core feature only becomes true when this package is loaded.
#
# Nothing else needs registering. check_lca_structure() runs through
# lca()'s own frmtmb_structure(check_frame =) slot, which frmtmb calls
# for whatever family the model names;
# frmtmb_register_frame_check() is for a check that must run on models
# the registering package's family does NOT appear in, which is not
# the case for either family here.
#
# Registering from .onLoad() rather than at top level is what a
# contributor outside frmtmb must do: by then every namespace is
# sealed, so the collation-order question that governs frmtmb's own
# in-package contributors does not arise.

#' @noRd
.onLoad <- function(libname, pkgname) {
  frmtmb_register_compat(features = c(hmm = "structure"),
                         rules = hmm_compat_rules)
  frmtmb_register_compat(features = c(lca = "structure"),
                         rules = lca_compat_rules)
  invisible()
}
