## REVIEW claim 7: RENDER the three changed Rd files and print the new
## sections, so the text is checked as the reader sees it rather than as
## the source spells it.
root <- "C:/Users/adf44/source/r/frmtmb-wt-thresrefit"
for (f in c("frm.Rd", "frm_bootstrap.Rd", "influence.frmtmb_fit.Rd")) {
  p <- file.path(root, "man", f)
  cat("\n================ ", f, " ================\n")
  tmp <- tempfile(fileext = ".txt")
  tools::Rd2txt(p, out = tmp)
  writeLines(readLines(tmp, warn = FALSE))
}
