# Lane adhess: one package's whole suite in ONE R session, every file in
# alphabetical order as covr and R CMD check run it, with every
# extension installed in the library path (lane-env.md, "Test
# isolation"). Prints the pass/fail line and every failing test.
#   Rscript dev/adhess-onesession.R <package>
p <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-adhess-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-adhess"
dir <- if (p == "frmtmb") file.path(wt, "tests/testthat") else
  file.path(wt, "extensions", p, "tests/testthat")
cat(p, "from", find.package(p), "| frmtmb from", find.package("frmtmb"), "\n")
res <- testthat::test_dir(dir, package = p, load_package = "none",
                          stop_on_failure = FALSE, reporter = "silent")
df <- as.data.frame(res)
cat(sprintf("ONE SESSION %s: files %d tests %d pass %d fail %d error %d skip %d warn %d\n",
            p, length(unique(df$file)), nrow(df), sum(df$passed),
            sum(df$failed), sum(df$error), sum(df$skipped), sum(df$warning)))
bad <- df[df$failed > 0 | df$error, c("file", "test", "failed", "error")]
if (nrow(bad)) print(bad, row.names = FALSE)
