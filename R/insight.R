# Registered methods for the insight package (the accessor layer under
# performance, parameters, and the rest of easystats). The default
# methods already work off our standard accessors; these fill the gaps
# the defaults cannot see: the random-effect formula split (which gates
# is_mixed_model and r2_nakagawa) and clean parameter extraction.

#' @exportS3Method insight::find_formula
find_formula.frmtmb_fit <- function(x, verbose = TRUE, ...) {
  rspec <- x$spec$responses[[1L]]
  dp <- rspec$dpars$mu %||% rspec$dpars[[1L]]
  out <- list(
    conditional = stats::as.formula(
      call("~", rspec$resp_expr, reformulas::RHSForm(dp$fixed)),
      env = rspec$formula_env
    )
  )
  bars <- lapply(dp$re %||% list(), function(rt) {
    stats::as.formula(call("~", rt$bar), env = rspec$formula_env)
  })
  if (length(bars)) {
    out$random <- if (length(bars) == 1L) bars[[1L]] else unname(bars)
  }
  structure(out, class = c("insight_formula", "list"))
}

#' @exportS3Method insight::find_random
find_random.frmtmb_fit <- function(x, split_nested = FALSE,
                                   flatten = FALSE, ...) {
  bks <- Filter(function(bk) bk$covstruct != "smooth",
                x$frame$re_blocks)
  grps <- unique(vapply(bks, `[[`, "", "group_name"))
  if (!length(grps)) return(NULL)
  if (flatten) grps else list(random = grps)
}

#' @exportS3Method insight::get_parameters
get_parameters.frmtmb_fit <- function(x, ...) {
  bd <- x$estimates$betad
  if (length(fx <- x$frame$betad_fixed_idx)) bd <- bd[-fx]
  est <- c(x$estimates$beta, bd)
  data.frame(
    Parameter = estimated_coef_names(x),
    Estimate = unname(est),
    Component = c(rep("conditional", length(x$estimates$beta)),
                  rep("dispersion", length(bd)))
  )
}

#' @exportS3Method insight::get_varcov
get_varcov.frmtmb_fit <- function(x, ...) {
  vcov(x)
}

#' @exportS3Method insight::find_statistic
find_statistic.frmtmb_fit <- function(x, ...) {
  "z-statistic"
}

#' @exportS3Method insight::link_inverse
link_inverse.frmtmb_fit <- function(x, ...) {
  rspec <- x$spec$responses[[1L]]
  dp <- rspec$dpars$mu %||% rspec$dpars[[1L]]
  dp$link$linkinv
}

#' @exportS3Method insight::link_function
link_function.frmtmb_fit <- function(x, ...) {
  rspec <- x$spec$responses[[1L]]
  dp <- rspec$dpars$mu %||% rspec$dpars[[1L]]
  dp$link$linkfun
}