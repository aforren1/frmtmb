# Fix round 1: car(M) without gr, brms's plainest spelling, is refused
# with brms's own reason (B1's probes), in whichever file of the tree
# named as the argument holds car_term_parse() (core's R/covstruct.R at
# the AB stage, frmtmb.spatial's R/car.R after the move).
root <- commandArgs(trailingOnly = TRUE)[1]
cand <- file.path(root, c("R/covstruct.R",
                          "extensions/frmtmb.spatial/R/car.R"))
old <- '  if (is.null(args[["M"]]) || is.null(args[["gr"]])) {
    frm_stop("car() needs an adjacency matrix and a grouping variable: ",
             "car(M, gr = g, type = \\"icar\\")", call. = FALSE)
  }
  # brms\'s gr = NA default (one location per observation) is deprecated
  # there and never supported here: the field\'s levels come from the
  # factor, so the call has to name one
  if (identical(args[["gr"]], NA) || identical(args[["gr"]], quote(NA))) {'
new <- '  if (is.null(args[["M"]])) {
    frm_stop("car() needs an adjacency matrix and a grouping variable: ",
             "car(M, gr = g, type = \\"icar\\")", call. = FALSE)
  }
  # brms\'s gr = NA default (one location per observation) is deprecated
  # there ("Using CAR terms without a grouping factor is deprecated") and
  # never supported here: the field\'s levels come from the factor, so the
  # call has to name one. car(M) is that default, written by omission
  if (is.null(args[["gr"]]) || identical(args[["gr"]], NA) ||
        identical(args[["gr"]], quote(NA))) {'
done <- 0L
for (p in cand[file.exists(cand)]) {
  s <- paste(readLines(p, warn = FALSE), collapse = "\n")
  if (!grepl(old, s, fixed = TRUE)) next
  writeLines(sub(old, new, s, fixed = TRUE), p)
  done <- done + 1L
}
stopifnot(done == 1L)
cat("fix-car applied to", root, "\n")
