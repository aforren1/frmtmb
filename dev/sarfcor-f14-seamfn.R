# Fix round 1, F14: seam_fn() takes a function of `...`. Run on the
# worktree and the stage-A copy: Rscript <this> <root>
root <- commandArgs(TRUE)[1]
p <- file.path(root, "R/covstruct-seam.R")
txt <- paste(readLines(p), collapse = "\n")
old <- "  if (!is.function(x) || length(formals(x)) < argn) {
    frm_stop(what, \" must be a function of \", argn, \" argument\","
new <- "  # a function of `...` takes any number of arguments
  if (!is.function(x) || (length(formals(x)) < argn &&
                            !\"...\" %in% names(formals(x)))) {
    frm_stop(what, \" must be a function of \", argn, \" argument\","
stopifnot(lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE))) == 1L)
txt <- sub(old, new, txt, fixed = TRUE)
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
cat("edited", root, "\n")
