# roxygenise core or extensions of the release tree against the library
# dev/rel077-install.sh names (REL077_LIB), so an extension's Rd sees the
# core it will be installed with. Adapted from dev/rel076-rox.R.
args <- commandArgs(trailingOnly = TRUE)
LIB <- Sys.getenv("REL077_LIB")
stopifnot(nzchar(LIB))
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- Sys.getenv("REL077_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-rel077")
for (p in args) {
  path <- if (p == "frmtmb") root else file.path(root, "extensions", p)
  cat("== roxygenise", p, "\n")
  roxygen2::roxygenise(path)
}
cat("ROX DONE\n")
