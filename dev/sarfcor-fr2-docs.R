# Fix round 2, (1): docs, compat notes, vignettes and NEWS say that
# newdata covers every row of the term's matrix, the dropped ones too.
# Part A. Rscript <this> <root>   (core NEWS.md is shared: run on both)
root <- commandArgs(TRUE)[1]
edit <- function(p, pairs) {
  p <- file.path(root, p)
  if (!file.exists(p)) return(invisible())
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}
note <- paste0(
  "newdata with one row per row of the term's matrix (the rows of the ",
  "data, the ones the fit dropped for missing values among them) is ",
  "read in their order, as brms reads it, so a counterfactual (x + 1 ",
  "everywhere) works and a row with a missing predictor gives NA; other ",
  "newdata is refused, as brms refuses it without a matrix over the new ",
  "rows (its newdata2).")
edit("R/autocor-seam.R", list(c(
  "#' fitted rows, so [fitted()], [predict()] and [simulate()] take
#' `newdata` with the same number of rows as the fit, read as the
#' fitted rows in their order (a counterfactual on the same units, as
#' brms takes it), and refuse other `newdata`, as brms refuses it
#' without a matrix over the new rows in its `newdata2`.",
  "#' rows of the data, so [fitted()], [predict()] and [simulate()] take
#' `newdata` with one row per row of the term's matrix, read in their
#' order (a counterfactual on the same units, as brms takes it). When
#' the fit dropped rows for missing values, `newdata` covers them too:
#' the term is built again over every row, its parameters held at
#' their natural values (through `bounds`, for a bounded one), and a
#' row with a missing predictor gives `NA`. Other `newdata` is refused,
#' as brms refuses it without a matrix over the new rows in its
#' `newdata2`.")))
edit("R/autocor.R", list(c(
  "#' [fitted()], [predict()] and [simulate()] take `newdata` with as many
#' rows as the fit, read as the fitted rows in their order, as brms
#' takes it, so a counterfactual such as `x + 1` everywhere works. Other
#' `newdata` is refused, because `M` covers the fitted rows only; brms
#' refuses it too unless a matrix over the new rows comes in its
#' `newdata2`. [frm_linpred()] gives `mu` on any new data. Pearson residuals divide by the marginal
#' sd, \\eqn{\\sigma \\sqrt{M_{ii}}}.",
  "#' [fitted()], [predict()] and [simulate()] take `newdata` with one row
#' per row of `M`, read in their order, as brms takes it, so a
#' counterfactual such as `x + 1` everywhere works; when the fit dropped
#' rows for missing values, `newdata` covers them too, and a row with a
#' missing predictor gives `NA`. Other `newdata` is refused; brms
#' refuses it too unless a matrix over the new rows comes in its
#' `newdata2`. [frm_linpred()] gives `mu` on any new data. Pearson
#' residuals divide by the marginal sd, \\eqn{\\sigma \\sqrt{M_{ii}}}.")))
edit("extensions/frmtmb.spatial/R/sar.R", list(c(
  "#' [fitted()], [predict()] and [simulate()] take `newdata` with as many
#' rows as the fit, read as the fitted locations in their order, as
#' brms takes it: a counterfactual such as `x + 1` everywhere gives the
#' lag model's spillover, \\eqn{A^{-1} \\eta_{new}}. Other `newdata` is
#' refused, because `M` covers the fitted rows only; brms refuses it
#' too unless a matrix over the new rows comes in its `newdata2`.",
  "#' [fitted()], [predict()] and [simulate()] take `newdata` with one row
#' per row of `M`, read as the locations in their order, as brms takes
#' it: a counterfactual such as `x + 1` everywhere gives the lag model's
#' spillover, \\eqn{A^{-1} \\eta_{new}}. When the fit dropped rows for
#' missing values, `newdata` covers them too, and A is then over every
#' row of `M`, at the fitted rho; a missing predictor gives `NA`, which
#' the lag form spreads to every row, as in brms. Other `newdata` is
#' refused; brms refuses it too unless a matrix over the new rows comes
#' in its `newdata2`.")))
for (k in c("compat.R", "zzz.R")) {
  p <- if (k == "compat.R") "R/compat.R" else
    "extensions/frmtmb.spatial/R/zzz.R"
  f <- file.path(root, p)
  x <- readLines(f)
  old <- paste0("newdata with as many rows as the fit is read as the ",
                "fitted rows in their order, as brms reads it, so a ",
                "counterfactual (x + 1 everywhere) works; other newdata is ",
                "refused, as brms refuses it without a matrix over the new ",
                "rows (its newdata2).")
  hit <- grepl(old, x, fixed = TRUE)
  stopifnot(sum(hit) == 3L)
  x <- gsub(old, note, x, fixed = TRUE)
  writeLines(x, f)
}
edit("NEWS.md", list(c(
  "They also take `newdata`
  with as many rows as the fit, read as the fitted rows in their order,
  as brms reads it, so a counterfactual (`x + 1` everywhere) works;",
  "They also take `newdata`
  with one row per row of `M`, read in their order, as brms reads it,
  so a counterfactual (`x + 1` everywhere) works; when the fit dropped
  rows for missing values, `newdata` covers them too, and a row with a
  missing predictor gives `NA`, as in brms;")))
edit("extensions/frmtmb.spatial/NEWS.md", list(c(
  "- `fitted()`, `predict()` and `simulate()` take `newdata` with as many
  rows as the fit, read as the fitted locations in their order, as
  brms reads it.",
  "- `fitted()`, `predict()` and `simulate()` take `newdata` with one row
  per row of `M`, read as the locations in their order, as brms reads
  it; when the fit dropped rows for missing values, `newdata` covers
  them too, and `I - rho M` is then over every row of `M`.")))
edit("extensions/frmtmb.spatial/vignettes/spatial.Rmd", list(c(
  "`simulate()` take `newdata` with as many rows as the fit, read as the
fitted locations in their order, as brms takes it,",
  "`simulate()` take `newdata` with one row per row of `M`, the rows a
fit dropped for missing values among them, read as the locations in
their order, as brms takes it,")))
cat("edited", root, "\n")
