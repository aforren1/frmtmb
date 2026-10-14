# Lane remeasure: does a SETUP line that reads a stale name pass the
# staleness on? brms_port() refuses an assertion that reads a name whose
# last assignment failed, but brms_setup() runs its line regardless and
# marks what it assigns as fresh. Constructed with the harness itself,
# in record mode, as the generated bin-2 files use it.
#
#   Rscript dev/remeasure-bin2-stale-setup.R \
#     > dev/remeasure-out/probe/bin2-stale-setup.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
`%||%` <- function(a, b) if (is.null(a)) b else a
source("tests/testthat/helper-brms-suite.R")
rec <- normalizePath("dev/remeasure-out/probe/bin2-stale-setup.tsv",
                     winslash = "/", mustWork = FALSE)
unlink(rec)
Sys.setenv(FRMTMB_BRMSPORT_RECORD = rec)
test_that("stale through a setup line", {
  brms_setup("x:1", a <- 1)
  brms_setup("x:2", a <- stop("the line brms's block would stop at"))
  # reads the stale a directly: refused
  brms_port("x:3", "pass", "", expect_equal(a, 1))
  # reads it through a setup line
  brms_setup("x:4", b <- a + 1)
  brms_port("x:5", "pass", "", expect_equal(b, 2))
  # the not-run form: a setup line reaching brms's namespace
  brms_setup("x:6", prep <- list(data = list(Y = 1:3)))
  brms_setup_not_run("x:7", "prep")
  brms_setup("x:8", n <- length(prep$data$Y))
  brms_port("x:9", "pass", "", expect_equal(n, 3))
})
x <- utils::read.delim(rec, header = FALSE, quote = "",
                       colClasses = "character")
for (k in seq_len(nrow(x))) {
  cat(x[k, 1], x[k, 3], "held=", x[k, 5], "msg=", substr(x[k, 7], 1, 90),
      "\n")
}
