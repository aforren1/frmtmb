# roxygenise core or extensions of the release tree against the library
# dev/rel073-install.sh names (REL073_LIB), so an extension's Rd sees the
# core it will be installed with. Adapted from dev/rel072-rox.R.
args <- commandArgs(trailingOnly = TRUE)
LIB <- Sys.getenv("REL073_LIB")
stopifnot(nzchar(LIB))
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-release"
for (p in args) {
  path <- if (p == "frmtmb") root else file.path(root, "extensions", p)
  cat("== roxygenise", p, "\n")
  roxygen2::roxygenise(path)
}
cat("ROX DONE\n")
