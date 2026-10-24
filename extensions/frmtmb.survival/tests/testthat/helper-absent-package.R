# Library stacks for a child process that reads a saved fit where an
# extension package cannot be used. frmtmb.survival and frmtmb.spatial
# keep identical copies of this file (tests/testthat/
# helper-absent-package.R); the master is frmtmb's
# dev/rel074-helper-absent-package.R.
#
# Two stacks, because R CMD check can build only one of them. It gives
# its tests one temporary library that holds the package under check
# and every dependency, then R's own library. So a stack that drops
# every library holding the package drops frmtmb's dependencies too,
# and frmtmb cannot load there; the test skips that arm and says why. A
# stack that puts an unloadable stub of the package in front of every
# library can be built anywhere, so the refusal stays covered under R
# CMD check.

# Run `code` in a fresh R process whose library paths are `libs`, and
# return its output lines (stdout and stderr).
#
# The child inherits R_LIBS_SITE, and `.libPaths(new)` appends the site
# libraries it names. covr sets R_LIBS_SITE to its whole library stack,
# the library holding the package included, so a plain `.libPaths(libs)`
# put the package back. `include.site` exists from R 4.2.0 on.
child_with_libs <- function(libs, code) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(sprintf(paste0("if ('include.site' %%in%% ",
                              "names(formals(.libPaths))) ",
                              ".libPaths(%s, include.site = FALSE) ",
                              "else .libPaths(%s)"),
                       paste0(deparse(libs), collapse = ""),
                       paste0(deparse(libs), collapse = "")),
               code), f)
  suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                           c("--vanilla", shQuote(f)), stdout = TRUE,
                           stderr = TRUE))
}

# Installed but unloadable: a `pkg` that is not a valid installed
# package (a DESCRIPTION and nothing else), in front of every library
# this session reads. Both stacks live in the session temporary
# directory, which R removes when the session ends; neither package
# declares withr, so no deferred cleanup.
unloadable_package_libs <- function(pkg) {
  broken <- tempfile("broken-lib")
  dir.create(file.path(broken, pkg), recursive = TRUE)
  writeLines(c(paste("Package:", pkg), "Version: 0.0.0"),
             file.path(broken, pkg, "DESCRIPTION"))
  c(broken, .libPaths())
}

# Not installed: this session's frmtmb, copied, in front of every
# library that does not hold `pkg`; the copy masks any other frmtmb
# there. Skips the calling test when frmtmb cannot load in that stack,
# which is the case under R CMD check (see the top of this file).
absent_package_libs <- function(pkg) {
  lib <- tempfile("absent-lib")
  dir.create(lib)
  file.copy(find.package("frmtmb"), lib, recursive = TRUE)
  rest <- Filter(function(l) !dir.exists(file.path(l, pkg)), .libPaths())
  libs <- c(lib, rest)
  probe <- child_with_libs(libs, c("suppressMessages(library(frmtmb))",
                                   "cat('FRMTMB LOADED\\n')"))
  testthat::skip_if_not(any(probe == "FRMTMB LOADED"),
                        paste("frmtmb's dependencies are only in a",
                              "library that also holds", pkg))
  libs
}
