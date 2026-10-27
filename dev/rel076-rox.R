# roxygenise core or extensions of the release tree against the library
# dev/rel076-install.sh names (REL076_LIB), so an extension's Rd sees the
# core it will be installed with. Adapted from dev/rel075-rox.R.
args <- commandArgs(trailingOnly = TRUE)
LIB <- Sys.getenv("REL076_LIB")
stopifnot(nzchar(LIB))
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- Sys.getenv("REL076_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-rel076")
for (p in args) {
  path <- if (p == "frmtmb") root else file.path(root, "extensions", p)
  cat("== roxygenise", p, "\n")
  roxygen2::roxygenise(path)
}
cat("ROX DONE\n")
