# codemeta.json has to describe the DESCRIPTION beside it.
#
# WHY THIS EXISTS. Nothing regenerates codemeta.json on its own, so it
# lagged behind the package at release after release: 0.50.0 while the
# package was at 0.64.0 (without `ordinal`, a Suggests by then), and
# 0.79.1 while the package was at 0.82.0. The release flow regenerates
# it with codemetar::write_codemeta() (dev/release/checklist.md), and
# this file is what makes a forgotten step fail instead of ship.
#
# The test reads the SOURCE tree, so it runs where the release runs its
# test files one per process and skips under R CMD check and covr, which
# hold no codemeta.json. A regex reads the file because a JSON parser is
# not a dependency of this package.

codemeta_root <- function() {
  root <- testthat::test_path("..", "..")
  desc <- file.path(root, "DESCRIPTION")
  ok <- file.exists(desc) && file.exists(file.path(root, "codemeta.json")) &&
    identical(unname(read.dcf(desc)[1L, "Package"]), "frmtmb")
  if (ok) root else NULL
}

test_that("codemeta.json carries the version of DESCRIPTION", {
  root <- codemeta_root()
  if (is.null(root)) skip("not run from the source tree")
  desc <- read.dcf(file.path(root, "DESCRIPTION"))
  cm <- readLines(file.path(root, "codemeta.json"), warn = FALSE)
  # the package's own version is the one at the top level, two spaces in;
  # a dependency's version constraint is nested deeper
  top <- grep('^  "version": "', cm, value = TRUE)
  expect_length(top, 1L)
  have <- sub('^  "version": "([^"]*)".*$', "\\1", top)
  expect_identical(
    have, unname(desc[1L, "Version"]),
    info = paste("regenerate it: codemetar::write_codemeta() at the",
                 "repository root"))
})

test_that("codemeta.json lists every dependency of DESCRIPTION", {
  root <- codemeta_root()
  if (is.null(root)) skip("not run from the source tree")
  desc <- read.dcf(file.path(root, "DESCRIPTION"))
  fields <- intersect(c("Depends", "Imports", "Suggests"), colnames(desc))
  want <- trimws(sub("[(].*", "", unlist(strsplit(
    paste(desc[1L, fields], collapse = ","), ","))))
  want <- sort(unique(want[nzchar(want)]))
  cm <- readLines(file.path(root, "codemeta.json"), warn = FALSE)
  id <- sub('^ *"identifier": "([^"]*)".*$', "\\1",
            grep('^ *"identifier": "', cm, value = TRUE))
  have <- sort(unique(setdiff(id, "frmtmb")))
  # the guard has to see the dependencies at all, or the comparison
  # below passes on a file it failed to read
  expect_gt(length(want), 10L)
  expect_identical(have, want,
                   info = paste("regenerate it: codemetar::write_codemeta()",
                                "at the repository root"))
})
