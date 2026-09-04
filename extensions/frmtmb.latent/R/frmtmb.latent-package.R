#' @keywords internal
"_PACKAGE"

# frmtmb is a Depends, so that one library(frmtmb.latent) call gives a
# user the frm() grammar that hmm() and lca() are families inside. The
# accessors are imported by name as well, because a namespace that is
# loaded and not attached reaches nothing through the search path, and
# these two files call every one of them unqualified.
#' @importFrom frmtmb as_frmtmb_family compat_rule_builder dpar_linpred
#' @importFrom frmtmb eval_dpars frame_block_of frmtmb_family
#' @importFrom frmtmb frmtmb_register_compat frmtmb_structure
#' @importFrom frmtmb latent_probs mixture_multimodal_refusals
#' @importFrom frmtmb mixture_posterior response_mean single_response
#' @importFrom frmtmb structure_supports_all
NULL
