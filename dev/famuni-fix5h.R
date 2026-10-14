p <- "extensions/frmtmb.sample/R/sample.R"
x <- paste(readLines(p), collapse = "\n")
old <- "mo_sample_entries <- function(entries, fit) {
  terms <- Filter(function(tm) !length(tm$blocked), gev_xi_chart_terms(fit))"
new <- "mo_sample_entries <- function(entries, fit = NULL) {
  # without the fit there is no GEV chart to read a prior in
  terms <- if (!is.null(fit)) {
    Filter(function(tm) !length(tm$blocked), gev_xi_chart_terms(fit))
  }"
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], p)
