# roxygenise core or extensions of the release tree against the library
# the install script names (REL071_LIB), so an extension's Rd sees the
# core it will be installed with.
args <- commandArgs(trailingOnly = TRUE)
LIB <- Sys.getenv("REL071_LIB")
stopifnot(nzchar(LIB))
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-release"
for (p in args) {
  path <- if (p == "frmtmb") root else file.path(root, "extensions", p)
  cat("== roxygenise", p, "\n")
  roxygen2::roxygenise(path)
}
cat("ROX DONE\n")
