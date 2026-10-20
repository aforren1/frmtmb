# Core's whole suite in ONE session, alphabetical, as R CMD check and
# covr run it. Log: dev/drop1raw-out/onesession.txt
.libPaths(c("C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
r <- testthat::test_dir("tests/testthat", package = "frmtmb",
                        load_package = "none", stop_on_failure = FALSE,
                        reporter = "silent")
df <- as.data.frame(r)
for (i in which(df$failed > 0 | df$error)) {
  cat("FAILED:", df$file[i], "::", df$test[i], "\n")
}
cat(sprintf("ONESESSION files=%d pass=%d fail=%d error=%d skip=%d warn=%d\n",
            length(unique(df$file)), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
