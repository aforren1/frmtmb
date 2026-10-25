# Fix round 2: the comments at the three newdata sites. Part A:
# Rscript <this> <root>
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
edit("R/predict.R", list(c(
  "      # newdata of the fitted rows, as brms takes it: the term's mean
      # over them, at newdata's mu (a counterfactual on the same units)",
  "      # newdata over the rows of the term's matrix, as brms takes it:
      # the term's mean over them, at newdata's mu (a counterfactual on
      # the same units)")))
edit("R/simulate-newdata.R", list(c(
  "    # newdata of the fitted rows keeps the fitted block, as brms keeps
    # its matrix for newdata of as many rows",
  "    # newdata over the rows of the term's matrix, as brms reads it: the
    # fitted block, or the block over every row of the matrix when the
    # fit dropped rows")))
edit("R/predict-brms.R", list(c(
  "    # a registry residual term (fcor(), sar()) takes newdata of the
    # fitted rows, as brms does, and refuses any other",
  "    # a registry residual term (fcor(), sar()) takes newdata over the
    # rows of its matrix, as brms does, and refuses any other")))
cat("edited", root, "\n")
