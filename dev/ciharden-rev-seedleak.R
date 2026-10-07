# Reviewer: does a seeded ported file leave the caller's RNG state as
# it found it? Runs the file between two draws from a known seed and
# compares with the same two draws around a file that touches no RNG.
.libPaths(c("C:/Users/adf44/source/r/wt-ciharden-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ciharden/tests/testthat"
ctl <- tempfile(fileext = ".R")
writeLines('test_that("control", expect_true(TRUE))', ctl)
around <- function(f) {
  set.seed(99)
  r <- testthat::test_file(f, package = "frmtmb",
                           env = testthat::test_env("frmtmb"),
                           reporter = "silent")
  df <- as.data.frame(r)
  list(next_draw = runif(1), passed = sum(df$passed),
       failed = sum(df$failed), skipped = sum(df$skipped))
}
c0 <- around(ctl)
for (f in c("test-brms-suite-standata.R", "test-brms-suite-priors.R")) {
  a <- around(file.path(wt, f))
  cat(sprintf("%s: next draw after file %.15f, control %.15f, same %s;",
              f, a$next_draw, c0$next_draw,
              identical(a$next_draw, c0$next_draw)),
      "pass", a$passed, "fail", a$failed, "skip", a$skipped, "\n")
}
