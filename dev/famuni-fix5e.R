p <- "extensions/frmtmb.sample/R/sample.R"
x <- paste(readLines(p), collapse = "\n")
old <- "  if (mo_soft) m <- mo_draws_to_chart(m, fit, laplace)
  if (gev_chart && !identical(obj, fit$obj)) {
    m <- gev_draws_to_xi(m, obj, fit, gev_terms, laplace)
  }"
new <- "  # first, while every column is still in the sampled coordinates the
  # objective reads
  if (gev_chart && !identical(obj, fit$obj)) {
    m <- gev_draws_to_xi(m, obj, fit, gev_terms, laplace)
  }
  if (mo_soft) m <- mo_draws_to_chart(m, fit, laplace)"
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], p)
