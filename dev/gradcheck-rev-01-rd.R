## Reviewer: render the two edited manual pages from the WORKTREE source
## (the lane library's help was built before the last roxygenise, so the
## worktree Rd is the authoritative copy) and grep the rendered text.
WT <- "C:/Users/adf44/source/r/frmtmb-wt-gradcheck"
for (f in c("man/diagnose.Rd", "man/frmtmb_control.Rd")) {
  cat("\n======================== ", f, "\n")
  out <- tempfile()
  tools::Rd2txt(file.path(WT, f), out = out, options = list(width = 78))
  txt <- readLines(out, warn = FALSE)
  writeLines(txt)
}
