# roxygenise core or extensions of the release tree against the library
# dev/rel075-install.sh names (REL075_LIB), so an extension's Rd sees the
# core it will be installed with. Adapted from dev/rel074-rox.R.
args <- commandArgs(trailingOnly = TRUE)
LIB <- Sys.getenv("REL075_LIB")
stopifnot(nzchar(LIB))
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- Sys.getenv("REL075_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-rel075")
for (p in args) {
  path <- if (p == "frmtmb") root else file.path(root, "extensions", p)
  cat("== roxygenise", p, "\n")
  roxygen2::roxygenise(path)
}
cat("ROX DONE\n")
