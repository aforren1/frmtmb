# Fix round 1, A8: the how-to names the exported registry helper, on the
# tree named as the argument.
root <- commandArgs(trailingOnly = TRUE)[1]
p <- file.path(root, "R/structure.R")
s <- paste(readLines(p, warn = FALSE), collapse = "\n")
old <- "#'     registry, in a test that calls `local_frmtmb_registries()` first."
new <- "#'     registry, in a test that calls [frmtmb_local_registries()]\n#'     first, so that every registry is restored when the test ends."
stopifnot(grepl(old, s, fixed = TRUE))
writeLines(sub(old, new, s, fixed = TRUE), p)
p <- file.path(root, "_pkgdown.yml")
s <- readLines(p, warn = FALSE)
if (!any(grepl("frmtmb_local_registries", s))) {
  i <- grep("^  - frmtmb_register_covstruct$", s)
  stopifnot(length(i) == 1L)
  s <- append(s, "  - frmtmb_local_registries", i)
  writeLines(s, p)
}
cat("fix-a8 applied to", root, "\n")
