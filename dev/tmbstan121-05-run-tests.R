# One test file, one R process, in one arm of this lane. Derived from
# dev/release/run-tests.R, which produced dev/suite-baseline.tsv, so the
# counts are comparable; only the library paths and the version
# assertion differ. Usage:
#   Rscript tmbstan121-05-run-tests.R <arm> <package> <path-to-test-file>
ROOT <- "C:/Users/adf44/source/r/tmbstan121-lib"
PIN <- "C:/Users/adf44/source/r/pinlib"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
p <- a[2]
f <- a[3]
want <- switch(arm,
  A = list(paths = c(file.path(ROOT, "A"), file.path(ROOT, "common"), PIN,
                     USER), sh = "2.32.10"),
  B = list(paths = c(file.path(ROOT, "Bsrc"), file.path(ROOT, "common"),
                     USER), sh = "2.39.1"),
  # The control: CRAN's tmbstan 1.2.0 built against 2.39.1, the defect.
  C = list(paths = c(file.path(ROOT, "C120"), file.path(ROOT, "common"),
                     USER), sh = "2.39.1", tv = "1.2.0"),
  stop("unknown arm ", arm))
.libPaths(want$paths)

# A wrong StanHeaders or a tmbstan resolved from the user library would
# make this arm measure a different configuration without saying so.
stopifnot(
  identical(format(utils::packageVersion("StanHeaders")), want$sh),
  identical(format(utils::packageVersion("tmbstan")),
            if (is.null(want$tv)) "1.2.1" else want$tv),
  identical(normalizePath(dirname(find.package("tmbstan"))),
            normalizePath(want$paths[1])),
  identical(normalizePath(dirname(find.package("frmtmb.sample"))),
            normalizePath(file.path(ROOT, "common"))))

suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))

res <- tryCatch(
  test_file(f, package = p, env = testthat::test_env(p),
            reporter = "silent"),
  error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })

if (!is.null(res)) {
  r <- as.data.frame(res)
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      "\n", sep = "")
  for (blk in res) {
    for (x in blk$results) {
      if (inherits(x, c("expectation_skip", "expectation_failure",
                        "expectation_error"))) {
        cat("DETAIL ", basename(f), " [", class(x)[1], "] ", blk$test,
            " :: ", gsub("[\r\n]+", " ", substr(conditionMessage(x), 1,
                                                 300)), "\n", sep = "")
      }
    }
  }
}
