# frmtmb's generics, re-exported.
#
# frmtmb keeps every generic it still has a method for, and this package
# registers its `frmtmb_draws` methods on those rather than defining a
# second generic of the same name. Re-exporting them means
# `library(frmtmb.sample); loo(ds)` works on its own: one generic per
# name in the session, whichever of the two packages is attached.
#
# The generics with no frmtmb method left - log_lik(), psis(),
# posterior_epred() and the rest - are DEFINED in this package, not
# re-exported, and are not on this page.

#' @importFrom frmtmb loo
#' @export
frmtmb::loo

#' @importFrom frmtmb waic
#' @export
frmtmb::waic

#' @importFrom frmtmb LOO
#' @export
frmtmb::LOO

#' @importFrom frmtmb WAIC
#' @export
frmtmb::WAIC

#' @importFrom frmtmb loo_compare
#' @export
frmtmb::loo_compare

#' @importFrom frmtmb bayes_R2
#' @export
frmtmb::bayes_R2

#' @importFrom frmtmb posterior_summary
#' @export
frmtmb::posterior_summary

#' @importFrom frmtmb as_draws
#' @export
frmtmb::as_draws

#' @importFrom frmtmb as_draws_matrix
#' @export
frmtmb::as_draws_matrix

#' @importFrom frmtmb as_draws_array
#' @export
frmtmb::as_draws_array

#' @importFrom frmtmb as_draws_df
#' @export
frmtmb::as_draws_df

#' @importFrom frmtmb as_draws_list
#' @export
frmtmb::as_draws_list

#' @importFrom frmtmb as_draws_rvars
#' @export
frmtmb::as_draws_rvars

#' @importFrom frmtmb ndraws
#' @export
frmtmb::ndraws

#' @importFrom frmtmb nchains
#' @export
frmtmb::nchains

#' @importFrom frmtmb niterations
#' @export
frmtmb::niterations

#' @importFrom frmtmb nvariables
#' @export
frmtmb::nvariables

#' @importFrom frmtmb hypothesis
#' @export
frmtmb::hypothesis

#' @importFrom frmtmb prior_summary
#' @export
frmtmb::prior_summary

#' @importFrom frmtmb variables
#' @export
frmtmb::variables

#' @importFrom frmtmb pp_check
#' @export
frmtmb::pp_check

#' @importFrom frmtmb conditional_effects
#' @export
frmtmb::conditional_effects

#' @importFrom frmtmb fixef
#' @export
frmtmb::fixef

#' @importFrom frmtmb ranef
#' @export
frmtmb::ranef

#' @importFrom frmtmb VarCorr
#' @export
frmtmb::VarCorr

#' @importFrom frmtmb ngrps
#' @export
frmtmb::ngrps

#' @importFrom frmtmb frm
#' @export
frmtmb::frm

#' @importFrom frmtmb bf
#' @export
frmtmb::bf

#' @importFrom frmtmb set_prior
#' @export
frmtmb::set_prior
