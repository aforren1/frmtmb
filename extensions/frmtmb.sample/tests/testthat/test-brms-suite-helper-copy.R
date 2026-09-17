# frmtmb.sample's helper-brms-suite.R is a copy of core's, written by
# dev/brmsport-gen.R. A drifted copy would run the sample half of the
# brms suite under different rules from the core half, so every object
# the two define must be the same. Core's file exists only in the source
# tree, which is where the gated tier runs; an installed-package check
# has no core tests to compare against.

skip_unless_brms_suite()

test_that("the sample helper is identical to core's, object by object", {
  core <- testthat::test_path("..", "..", "..", "..", "tests", "testthat",
                              "helper-brms-suite.R")
  copy <- testthat::test_path("helper-brms-suite.R")
  skip_if_not(file.exists(core), "core's helper is not in this tree")
  a <- new.env()
  b <- new.env()
  sys.source(core, envir = a, keep.source = FALSE)
  sys.source(copy, envir = b, keep.source = FALSE)
  expect_identical(sort(ls(a, all.names = TRUE)),
                   sort(ls(b, all.names = TRUE)))
  for (nm in ls(a, all.names = TRUE)) {
    x <- get(nm, envir = a)
    y <- get(nm, envir = b)
    if (is.function(x)) {
      expect_identical(deparse(x), deparse(y), info = nm)
    } else if (!is.environment(x)) {
      expect_identical(x, y, info = nm)
    }
  }
})
