# One-off edit of R/conditional-effects.R for item 2 (offset variables are
# not default displays). Kept as the record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/conditional-effects.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("#' the only place the name survives.
#'
#' @noRd
ce_lp_vars <- function(lp) {
  v <- if (is.null(lp[[\"terms\"]])) {
    character(0)
  } else {
    all.vars(stats::delete.response(lp[[\"terms\"]]))
  }",
"#' the only place the name survives. `offsets = FALSE` leaves out a
#' variable that only an `offset()` term reads, as brms's default
#' displays do.
#'
#' @noRd
ce_lp_vars <- function(lp, offsets = TRUE) {
  v <- if (is.null(lp[[\"terms\"]])) {
    character(0)
  } else {
    tt <- stats::delete.response(lp[[\"terms\"]])
    tv <- as.list(attr(tt, \"variables\"))[-1L]
    if (!offsets && length(attr(tt, \"offset\"))) tv <- tv[-attr(tt, \"offset\")]
    unique(unlist(lapply(tv, all.vars)))
  }")
rep1("ce_plot_vars <- function(x, rspec, lp, resp, seen = character(0)) {
  v <- ce_lp_vars(lp)",
"ce_plot_vars <- function(x, rspec, lp, resp, seen = character(0),
                         offsets = TRUE) {
  v <- ce_lp_vars(lp, offsets)")
rep1("        v <- c(v, ce_plot_vars(x, rspec, lpn, resp, c(seen, np)))",
"        v <- c(v, ce_plot_vars(x, rspec, lpn, resp, c(seen, np),
                               offsets))")
rep1("ce_plot_vars_any <- function(x, rspec, resp) {
  v <- character(0)
  for (dp in names(rspec$dpars)) {
    lpn <- x$frame[[\"linpreds\"]][[linpred_key(resp, dp)]]
    if (!is.null(lpn)) v <- c(v, ce_plot_vars(x, rspec, lpn, resp))",
"ce_plot_vars_any <- function(x, rspec, resp, offsets = TRUE) {
  v <- character(0)
  for (dp in names(rspec$dpars)) {
    lpn <- x$frame[[\"linpreds\"]][[linpred_key(resp, dp)]]
    if (!is.null(lpn)) {
      v <- c(v, ce_plot_vars(x, rspec, lpn, resp, offsets = offsets))
    }")
rep1("  vars <- ce_plot_vars(x, rspec, lp, resp)
  # a model whose SELECTED predictor has no terms still has covariates
  # somewhere (bf(y ~ 1, theta1 ~ x)): naming the one predictor the
  # search looked at was a refusal to draw a model that has something
  # to draw
  if (is.null(effects) && !length(intersect(vars, names(base)))) {
    vars <- ce_plot_vars_any(x, rspec, resp)
  }",
"  # an offset's variable is held at its mean like any other, but it is
  # not a default display: brms draws the terms, and an offset is not
  # one
  vars <- ce_plot_vars(x, rspec, lp, resp, offsets = !is.null(effects))
  # a model whose SELECTED predictor has no terms still has covariates
  # somewhere (bf(y ~ 1, theta1 ~ x)): naming the one predictor the
  # search looked at was a refusal to draw a model that has something
  # to draw
  if (is.null(effects) && !length(intersect(vars, names(base)))) {
    vars <- ce_plot_vars_any(x, rspec, resp, offsets = FALSE)
  }")
writeLines(x, p)
