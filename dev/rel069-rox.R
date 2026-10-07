args <- commandArgs(trailingOnly = TRUE)
REL <- "C:/Users/adf44/source/r/rellib-r7"
.libPaths(c(REL, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-release"
for (p in args) {
  path <- if (p == "frmtmb") root else file.path(root, "extensions", p)
  cat("== roxygenise", p, "\n")
  roxygen2::roxygenise(path)
}
cat("ROX DONE\n")
