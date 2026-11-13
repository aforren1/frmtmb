# Lane arspat, fix round (review MINOR-8): which test of a core test file
# raises frmtmb_autocor_gap without catching it. One file per call.
#   Rscript dev/arspat-gapnoise.R <test file> ...
.libPaths(c("C:/Users/adf44/source/r/wt-arspat-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
Sys.setenv(NOT_CRAN = "true")
for (f in commandArgs(trailingOnly = TRUE)) {
  where <- character(0)
  withCallingHandlers(
    suppressWarnings(testthat::test_file(
      file.path("tests/testthat", f), package = "frmtmb",
      env = testthat::test_env("frmtmb"), reporter = "silent")),
    frmtmb_autocor_gap = function(m) {
      # the innermost test_that() on the stack names the test
      cl <- sys.calls()
      desc <- NA_character_
      for (c1 in rev(cl)) {
        if (is.call(c1) && identical(c1[[1L]], as.name("test_that"))) {
          desc <- as.character(c1[[2L]])
          break
        }
      }
      where <<- c(where, desc)
    })
  cat(sprintf("%s: %d uncaught or caught gap messages\n", f, length(where)))
  if (length(where)) print(table(where))
}
