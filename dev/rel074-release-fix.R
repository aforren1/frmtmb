# Release fix of the 0.74.0 consolidation, sourced by
# dev/rel074-resolve.R (which defines rd(), wr() and edit_file()).
#
# frmtmb.survival's test-saved-fits.R built only a "not installed"
# library stack, by dropping every library that holds the package. R CMD
# check gives its tests one temporary library holding the package and
# every dependency, so there the child could not load frmtmb and 7
# expectations failed (Status: 1 ERROR; dev/round-20261023.md, "R CMD
# check --as-cran"; dev/rel074-log/survival-check-probe.txt).
# frmtmb.spatial's test-guardrail.R already had both stacks and the
# probe. Both move into one helper that both packages copy
# (dev/rel074-helper-absent-package.R); frmtmb.survival gains the
# "installed but unloadable" arm, which R CMD check can build, and only
# the "not installed" arm may skip. Approved by the coordinator as a
# release fix, not a lane change.
hp <- rd("dev/rel074-helper-absent-package.R")
for (p in c("frmtmb.survival", "frmtmb.spatial")) {
  wr(hp, file.path("extensions", p, "tests/testthat",
                   "helper-absent-package.R"))
}

ts <- "extensions/frmtmb.survival/tests/testthat/test-saved-fits.R"
edit_file(ts, r"---[# this package's name. The child process sees a library holding a copy
# of the installed frmtmb and none of the libraries that hold this
# package.

saved_child <- function(lib, code) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(sprintf(".libPaths(%s)", paste0(deparse(lib),
                                                collapse = "")),
               code), f)
  suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                           c("--vanilla", shQuote(f)), stdout = TRUE,
                           stderr = TRUE))
}
]---", r"---[# this package's name. The child process reads the fits in two library
# stacks (helper-absent-package.R): this package installed but
# unloadable, which every setup can build, and this package not
# installed, which R CMD check cannot build and so skips.
]---")
edit_file(ts, r"---[  # a library with frmtmb alone, ahead of every library that does not
  # hold this package
  lib <- tempfile("lib")
  dir.create(lib)
  on.exit(unlink(lib, recursive = TRUE), add = TRUE)
  file.copy(system.file(package = "frmtmb"), lib, recursive = TRUE)
  rest <- Filter(function(l) {
    !dir.exists(file.path(l, "frmtmb.survival"))
  }, .libPaths())
  out <- saved_child(c(lib, rest), c(
    "suppressMessages(library(frmtmb))",
    "cat('SURVIVAL ', nzchar(system.file(package = 'frmtmb.survival')),",
    "    '\\n')",]---", r"---[  every_method <- function(libs) {
  out <- child_with_libs(libs, c(
    "suppressMessages(library(frmtmb))",
    "cat('HAVE', requireNamespace('frmtmb.survival', quietly = TRUE),",
    "    '\\n')",]---")
edit_file(ts, r"---[  expect_true(any(out == "CHILDOK"), info = paste(out, collapse = "\n"))
  skip_if(any(grepl("^SURVIVAL +TRUE", out)),
          "this package is visible to the child; nothing to check")]---",
r"---[  expect_true(any(out == "CHILDOK"), info = paste(out, collapse = "\n"))
  # the stack is the one meant: the child cannot use this package
  expect_true(any(out == "HAVE FALSE "), info = paste(out, collapse = "\n"))]---")
edit_file(ts, r"---[  expect_equal(ll, as.numeric(logLik(fc)), tolerance = 1e-12)
})]---", r"---[  expect_equal(ll, as.numeric(logLik(fc)), tolerance = 1e-12)
  }
  # Installed but unloadable: built under R CMD check too, so the
  # refusal is covered there
  every_method(unloadable_package_libs("frmtmb.survival"))
  # Not installed: skipped where frmtmb cannot load without the library
  # that holds this package (R CMD check)
  every_method(absent_package_libs("frmtmb.survival"))
})]---")
# the body of every_method() is the old test body, one level deeper
x <- rd(ts)
a <- grep("^  every_method <- function[(]libs[)] [{]$", x)
z <- grep("^  # Installed but unloadable: built under R CMD check too",
          x) - 1L
stopifnot(length(a) == 1L, length(z) == 1L, x[z] == "  }", z > a)
x[(a + 1L):(z - 1L)] <- ifelse(nzchar(x[(a + 1L):(z - 1L)]),
                               paste0("  ", x[(a + 1L):(z - 1L)]), "")
wr(x, ts)
cat("indented the body of every_method() in", ts, "\n")

tg <- "extensions/frmtmb.spatial/tests/testthat/test-guardrail.R"
edit_file(tg, r"---[# The same, with the child's library paths given.
spatial_child_lib <- function(libs, code) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(sprintf(".libPaths(%s)", paste0(deparse(libs), collapse = "")),
               "suppressMessages(library(frmtmb))", code), f)
  suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                           c("--vanilla", shQuote(f)), stdout = TRUE,
                           stderr = TRUE))
}

]---", "")
edit_file(tg, r"---[    out <- spatial_child_lib(libs, c(
      sprintf(]---", r"---[    out <- child_with_libs(libs, c(
      "suppressMessages(library(frmtmb))",
      sprintf(]---")
edit_file(tg, r"---[  # Installed but unloadable: a frmtmb.spatial that is not a valid
  # installed package, in front of every library this session reads.
  broken <- tempfile("spatial-broken-lib")
  dir.create(file.path(broken, "frmtmb.spatial"), recursive = TRUE)
  on.exit(unlink(broken, recursive = TRUE), add = TRUE)
  writeLines(c("Package: frmtmb.spatial", "Version: 0.0.0"),
             file.path(broken, "frmtmb.spatial", "DESCRIPTION"))
  every_method(c(broken, .libPaths()))
  # Not installed: this session's frmtmb, copied, in front of every
  # library without frmtmb.spatial; the copy masks any other frmtmb there.
  # R CMD check puts frmtmb's dependencies in the library that holds
  # frmtmb.spatial, so there the arm cannot be built and is skipped.
  lib <- tempfile("spatial-absent-lib")
  dir.create(lib)
  on.exit(unlink(lib, recursive = TRUE), add = TRUE)
  file.copy(find.package("frmtmb"), lib, recursive = TRUE)
  rest <- Filter(function(l) {
    !dir.exists(file.path(l, "frmtmb.spatial"))
  }, .libPaths())
  probe <- spatial_child_lib(c(lib, rest), "cat('FRMTMB LOADED\\n')")
  skip_if_not(any(probe == "FRMTMB LOADED"),
              paste("frmtmb's dependencies are only in a library that",
                    "also holds frmtmb.spatial"))
  every_method(c(lib, rest))]---",
r"---[  # Installed but unloadable, then not installed
  # (helper-absent-package.R); only the second may skip, under R CMD
  # check, which cannot build it
  every_method(unloadable_package_libs("frmtmb.spatial"))
  every_method(absent_package_libs("frmtmb.spatial"))]---")
# the fixes the release review asked for
source("dev/rel074-release-review-fix.R", local = TRUE)
cat("RELEASE FIX APPLIED\n")
