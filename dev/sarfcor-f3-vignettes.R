# Fix round 1, F3: the vignettes say what newdata on the fitted rows
# does. Part A only: Rscript <this> <root>
root <- commandArgs(TRUE)[1]
edit <- function(p, pairs) {
  p <- file.path(root, p)
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}
edit("vignettes/brms-migration.Rmd", list(
  c("a predicted `sigma` gives `D M D`, as in brms; no `newdata`, as brms without `newdata2` |",
    "a predicted `sigma` gives `D M D`, as in brms; `newdata` only on the fitted rows, as brms without `newdata2` |")))
edit("extensions/frmtmb.spatial/vignettes/spatial.Rmd", list(
  c("`se()` and `mi()` are refused, and so is `newdata` in `fitted()`,
`predict()` and `simulate()`: `M` describes the fitted rows only.
brms refuses the same calls unless its `newdata2` holds a matrix over
the new rows.",
    "`se()` and `mi()` are refused. `fitted()`, `predict()` and
`simulate()` take `newdata` with as many rows as the fit, read as the
fitted locations in their order, as brms takes it, so a counterfactual
such as `x + 1` everywhere gives the lag model's spillover. Other
`newdata` is refused, because `M` describes the fitted rows only; brms
refuses it too unless its `newdata2` holds a matrix over the new
rows.")))
cat("edited", root, "\n")
