# Fix round 2: the shared NEWS files. Core NEWS's newdata sentence on
# the stage-A copy (the worktree got it from dev/sarfcor-fr2-docs.R),
# and frmtmb.sample's rounded announcement on both.
# Rscript <this> <worktree> <stageA>
args <- commandArgs(TRUE)
edit <- function(root, p, pairs) {
  p <- file.path(root, p)
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}
edit(args[2], "NEWS.md", list(c(
  "They also take `newdata`
  with as many rows as the fit, read as the fitted rows in their order,
  as brms reads it, so a counterfactual (`x + 1` everywhere) works;",
  "They also take `newdata`
  with one row per row of `M`, read in their order, as brms reads it,
  so a counterfactual (`x + 1` everywhere) works; when the fit dropped
  rows for missing values, `newdata` covers them too, and a row with a
  missing predictor gives `NA`, as in brms;")))
for (root in args) {
  edit(root, "extensions/frmtmb.sample/NEWS.md", list(c(
    "`rho` over its interval, which `frm_sample()` announces as
  `uniform(lower, upper)`.",
    "`rho` over its interval, which `frm_sample()` announces as
  `uniform(lower, upper)`, the ends shown to 6 digits.")))
}
cat("ok\n")
