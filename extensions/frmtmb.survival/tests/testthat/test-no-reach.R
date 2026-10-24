# Guardrail 2 of frmtmb's dev/core-boundary.md: this package uses only
# frmtmb's exported functions. A reach into the namespace would tie it
# to an internal that core may change without notice, which is what the
# seams it uses (simplex_pars =, aterm_data =, post$curve_check) exist
# to avoid. Checked on the package sources and on its tests, which are
# where a reach is easiest to write.

test_that("no source or test file of this package reaches into frmtmb", {
  root <- testthat::test_path("..", "..")
  desc <- file.path(root, "DESCRIPTION")
  is_src <- file.exists(desc) && dir.exists(file.path(root, "R")) &&
    !dir.exists(file.path(root, "Meta")) &&
    identical(unname(read.dcf(desc, fields = "Package")[1L, 1L]),
              "frmtmb.survival")
  skip_if_not(is_src,
              "package sources are not available (installed-package run)")
  files <- c(list.files(file.path(root, "R"), pattern = "[.][Rr]$",
                        full.names = TRUE),
             list.files(file.path(root, "tests"), pattern = "[.][Rr]$",
                        full.names = TRUE, recursive = TRUE),
             list.files(file.path(root, "vignettes"), pattern = "[.]Rmd$",
                        full.names = TRUE))
  expect_gt(length(files), 10L)
  # spelled in pieces so that this file does not match itself
  pats <- c(paste0("frmtmb", ":::"), paste0("getFrom", "Namespace"),
            paste0("asNamespace", "[(]['\"]frmtmb['\"][)]"))
  hits <- character(0)
  for (f in files) {
    x <- readLines(f, warn = FALSE)
    for (p in pats) {
      i <- grep(p, x)
      if (length(i)) {
        hits <- c(hits, paste0(basename(f), ":", i, ": ", trimws(x[i])))
      }
    }
  }
  expect_identical(hits, character(0),
                   info = paste(hits, collapse = "\n"))
  # the guard present: each pattern matches the reach it describes
  probe <- c(paste0("frmtmb", ":::family_ctor"),
             paste0("utils::getFrom", "Namespace('x', 'frmtmb')"),
             paste0("asNamespace", "(\"frmtmb\")"))
  expect_true(all(vapply(seq_along(pats), function(k) {
    grepl(pats[[k]], probe[[k]])
  }, NA)))
})
