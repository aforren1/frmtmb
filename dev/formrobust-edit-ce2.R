# One-off edit of R/conditional-effects.R: rename the offsets flag to
# `rsv` and leave out brms's reserved column of ones as well. Kept as
# the record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/conditional-effects.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("#' the only place the name survives. `offsets = FALSE` leaves out a
#' variable that only an `offset()` term reads, as brms's default
#' displays do.",
"#' the only place the name survives. `rsv = FALSE` leaves out what
#' brms's default displays leave out: a variable that only an `offset()`
#' term reads, and the column of ones of brms's deprecated
#' `0 + intercept` (`rsv_vars()` in brms).")
rep1("ce_lp_vars <- function(lp, offsets = TRUE) {", "ce_lp_vars <- function(lp, rsv = TRUE) {")
rep1("    if (!offsets && length(attr(tt, \"offset\"))) tv <- tv[-attr(tt, \"offset\")]
    unique(unlist(lapply(tv, all.vars)))
  }",
"    if (!rsv && length(attr(tt, \"offset\"))) tv <- tv[-attr(tt, \"offset\")]
    v <- unique(unlist(lapply(tv, all.vars)))
    if (!rsv && isTRUE(lp[[\"rsv_lower\"]])) v <- setdiff(v, \"intercept\")
    v
  }")
rep1("                         offsets = TRUE) {
  v <- ce_lp_vars(lp, offsets)", "                         rsv = TRUE) {
  v <- ce_lp_vars(lp, rsv)")
rep1("                               offsets))", "                               rsv))")
rep1("ce_plot_vars_any <- function(x, rspec, resp, offsets = TRUE) {", "ce_plot_vars_any <- function(x, rspec, resp, rsv = TRUE) {")
rep1("      v <- c(v, ce_plot_vars(x, rspec, lpn, resp, offsets = offsets))", "      v <- c(v, ce_plot_vars(x, rspec, lpn, resp, rsv = rsv))")
rep1("  vars <- ce_plot_vars(x, rspec, lp, resp, offsets = !is.null(effects))", "  vars <- ce_plot_vars(x, rspec, lp, resp, rsv = !is.null(effects))")
rep1("    vars <- ce_plot_vars_any(x, rspec, resp, offsets = FALSE)", "    vars <- ce_plot_vars_any(x, rspec, resp, rsv = FALSE)")
rep1("  # an offset's variable is held at its mean like any other, but it is
  # not a default display: brms draws the terms, and an offset is not
  # one", "  # an offset's variable is held at its mean like any other, but it is
  # not a default display: brms draws the terms, and an offset is not
  # one, nor is its deprecated reserved `intercept`")
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
close(con)
