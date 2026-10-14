p <- "extensions/frmtmb.sample/R/sample.R"
x <- paste(readLines(p), collapse = "\n")
r <- function(old, new) {
  if (!grepl(old, x, fixed = TRUE)) stop("not found: ", substr(old, 1, 60))
  x <<- sub(old, new, x, fixed = TRUE)
}
# 1. the predicted / scalar tell, from the formula as brms reads it
r("#' Whether `lp` is the `xi` of a `gen_extreme_value()` response with a
#' predictor of its own, which is where brms's default on `xi` is a
#' prior on an intercept frmtmb also has. A scalar `xi` is not: brms
#' samples it as `tmp_xi`, scaled into the support the data leave, and
#' frmtmb has no such parameter.
#'
#' @noRd
gev_xi_predicted <- function(rspec, lp) {
  identical(rspec$family[[\"family\"]], \"gen_extreme_value\") &&
    identical(lp[[\"dpar\"]], \"xi\") &&
    (ncol(lp[[\"X\"]]) > 1L || !is.null(lp[[\"Z\"]]))
}",
"#' Whether `lp` is the `xi` of a `gen_extreme_value()` response with a
#' formula of its own, `xi ~ 1` included: brms then samples `xi` on its
#' `log1p` link with `normal(0, 4)` on the intercept. Without a formula
#' brms samples `tmp_xi` with `normal(0, 2.5)` (`gev_xi_scalar()`). The
#' two spellings have the same design, one intercept column, so the tell
#' is the formula, which `dpar_has_predictor()` reads; a design test
#' gave `xi ~ 1` the scalar's note and no default (review fix 6).
#'
#' @noRd
gev_xi_predicted <- function(spec, rspec, lp) {
  identical(rspec$family[[\"family\"]], \"gen_extreme_value\") &&
    identical(lp[[\"dpar\"]], \"xi\") && is.null(lp[[\"constant\"]]) &&
    dpar_has_predictor(spec, lp)
}

#' Whether `lp` is a `gen_extreme_value()` `xi` without a formula.
#'
#' @noRd
gev_xi_scalar <- function(spec, rspec, lp) {
  identical(rspec$family[[\"family\"]], \"gen_extreme_value\") &&
    identical(lp[[\"dpar\"]], \"xi\") && is.null(lp[[\"constant\"]]) &&
    is.null(lp[[\"equate\"]]) && !dpar_has_predictor(spec, lp)
}")
r("    } else if (gev_xi_predicted(rspec, lp)) {",
  "    } else if (gev_xi_scalar(fit$spec, rspec, lp)) {
      # brms's default on a scalar xi is normal(0, 2.5) on tmp_xi, and
      # frm_sample() samples that xi through tmp_xi as brms does
      # (gev_xi_chart_frame()), where this density is read on tmp_xi.
      # Where the chart cannot be used it is a density on xi itself
      add(natural_dpar_prior(\"normal(0, 2.5)\", \"xi\", resp = rs))
    } else if (gev_xi_predicted(fit$spec, rspec, lp)) {")
# 2. the note: only where the chart is blocked
r("    xi_lp <- Filter(function(lp) {
      identical(lp[[\"resp\"]], rspec$resp_name) &&
        identical(lp[[\"dpar\"]], \"xi\") && is.null(lp[[\"constant\"]]) &&
        identical(rspec$family[[\"family\"]], \"gen_extreme_value\") &&
        !gev_xi_predicted(rspec, lp)
    }, fit$frame[[\"linpreds\"]])
    # a prior of the user's on xi is the one the note asks for
    xi_own <- any(vapply(unclass(user_pl %||% list()), function(s) {
      identical(s$class, \"xi\") ||
        (identical(s$class, \"Intercept\") && identical(s$dpar, \"xi\"))
    }, NA))
    if (length(xi_lp) && !xi_own) {
      notes <- c(notes, paste0(
        \"no default for xi\", of, \": brms puts normal(0, 2.5) on tmp_xi, \",
        \"its xi scaled into the support the data leave, which frmtmb does \",
        \"not have. A flat xi is improper on its log1p scale, so give it \",
        \"a proper prior, such as set_prior(\\\"normal(0, 2.5)\\\", \",
        \"class = \\\"xi\\\")\"))
    }",
"    for (tm in gev_xi_chart_terms(fit)) {
      if (!identical(tm$resp, rspec$resp_name) || !length(tm$blocked)) next
      notes <- c(notes, gev_xi_blocked_text(tm, of))
    }")
writeLines(strsplit(x, "\n")[[1]], p)
