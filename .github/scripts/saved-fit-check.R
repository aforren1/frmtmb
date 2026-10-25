# Runs frmtmb's test-extension-home-mixed.R where every extension is
# installed, and fails on any failure, error or skip. Called by
# pkgdown.yaml, the one job that installs every extension.
#
# Two of the file's blocks save a fit in a child process, a cox() and
# car() fit with frmtmb.survival and frmtmb.spatial loaded and a car()
# and sar() fit with frmtmb.spatial loaded, then read it in child
# processes where those packages are absent. Its helper
# (helper-absent-package.R) builds the "not installed" stack by
# dropping every library that holds them; on a
# runner they sit in the one library that also holds frmtmb's
# dependencies, so that arm would skip. They are moved to a library of
# their own, beside that one so that the move is a rename, for the
# duration of the test, and moved back whatever happens.
#
#   Rscript .github/scripts/saved-fit-check.R
# Run from the repository root, with frmtmb and every extension
# installed.
moved <- c("frmtmb.survival", "frmtmb.spatial")
home <- vapply(moved, function(p) dirname(find.package(p)), "")
ext <- file.path(dirname(home[[1L]]), "frmtmb-saved-fit-ext-lib")
dir.create(ext, showWarnings = FALSE)
back <- function() {
  for (p in moved) {
    if (dir.exists(file.path(ext, p))) {
      file.rename(file.path(ext, p), file.path(home[[p]], p))
    }
  }
  # only when empty: a package left in it is evidence for the error
  if (!length(list.files(ext, all.files = TRUE, no.. = TRUE))) {
    unlink(ext, recursive = TRUE)
  }
}
r <- tryCatch({
  for (p in moved) {
    if (!file.rename(file.path(home[[p]], p), file.path(ext, p))) {
      stop("could not move ", p, " out of ", home[[p]])
    }
  }
  .libPaths(c(ext, .libPaths()))
  Sys.setenv(NOT_CRAN = "true")
  as.data.frame(testthat::test_file(
    "tests/testthat/test-extension-home-mixed.R", package = "frmtmb",
    load_package = "installed", reporter = "summary"))
}, finally = back())
left <- vapply(moved, function(p) dir.exists(file.path(home[[p]], p)), NA)
print(r[, c("test", "nb", "failed", "skipped", "error")])
bad <- any(r$skipped) || any(r$error) || sum(r$failed) > 0 ||
  sum(r$nb) == 0 || !all(left)
if (bad) {
  stop("test-extension-home-mixed.R did not pass with every extension ",
       "installed", if (!all(left)) " (a package was not moved back)")
}
