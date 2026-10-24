# The contributed compatibility rows, with every extension loaded.
#
# test-compat.R checks core's rows in whatever session it is run in,
# and it must not depend on which extensions an earlier file loaded
# (covr and R CMD check run every file in one session). What the
# contributed rows promise can only be checked with every owner
# loaded, so this file loads every one in a CHILD process and leaves
# this session as it found it. helper-compat-contract.R says what is
# checked.
#
# It skips unless all of them are installed: under R CMD check, which
# hides packages the DESCRIPTION does not name, and on the core
# coverage job, which removes them. The release gate's one-session runs
# (dev/release/run-session.R) and a local run of this suite reach it.

compat_extensions <- c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent",
                       "frmtmb.learn", "frmtmb.ode", "frmtmb.sample",
                       "frmtmb.survival",
                       "frmtmb.spatial")

compat_run_child <- function(code) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  head <- sprintf(".libPaths(%s)",
                  paste0(deparse(.libPaths()), collapse = ""))
  writeLines(c(head, code), f)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(f)), stdout = TRUE, stderr = TRUE))
  out
}

test_that("every contributed rule resolves, wins a pair and names a
           function, with every extension loaded", {
  skip_on_cran()
  # system.file() and not requireNamespace(): this session must not load
  # them
  have <- vapply(compat_extensions, function(p) {
    nzchar(system.file(package = p))
  }, TRUE)
  skip_if_not(all(have), paste("not installed:",
                               toString(compat_extensions[!have])))
  helper <- normalizePath(testthat::test_path("helper-compat-contract.R"),
                          "/")
  out <- compat_run_child(c(
    "suppressMessages(library(frmtmb))",
    sprintf("for (p in %s) loadNamespace(p)",
            paste0(deparse(compat_extensions), collapse = "")),
    sprintf("source(%s)", deparse(helper)),
    "v <- compat_contract_violations()",
    "for (x in v) cat('VIOLATION ', x, '\\n', sep = '')",
    "cat('CHILDOK\\n')"))
  expect_true(any(out == "CHILDOK"), info = paste(out, collapse = "\n"))
  expect_equal(sub("^VIOLATION ", "", grep("^VIOLATION ", out,
                                           value = TRUE)),
               character(0))
})
