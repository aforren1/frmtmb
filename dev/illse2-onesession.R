# Lane illse2: core's whole suite in ONE R session, every file in
# alphabetical order as covr and R CMD check run it, with every
# extension installed in the library path (lane-env.md, "Test
# isolation"). Prints the pass/fail line and every failing test.
#   Rscript dev/illse2-onesession.R
.libPaths(c("C:/Users/adf44/source/r/wt-illse2-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
res <- testthat::test_dir(
  "C:/Users/adf44/source/r/frmtmb-wt-illse2/tests/testthat",
  package = "frmtmb", load_package = "none", stop_on_failure = FALSE,
  reporter = "silent")
df <- as.data.frame(res)
cat(sprintf("ONE SESSION: files %d tests %d pass %d fail %d error %d skip %d warn %d\n",
            length(unique(df$file)), nrow(df), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
bad <- df[df$failed > 0 | df$error, c("file", "test", "failed", "error")]
if (nrow(bad)) print(bad, row.names = FALSE)
