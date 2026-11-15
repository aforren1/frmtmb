# Short check, item 8: one test file per process on the lane library.
#   NOT_CRAN=true Rscript dev/reviews/arspat2-short-test.R <pkg> <file>
a <- commandArgs(TRUE)
p <- a[1]
f <- a[2]
.libPaths(c("C:/Users/adf44/source/r/wt-arspat2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", dirname(system.file(package = "frmtmb")), "|", p, "from",
    dirname(system.file(package = p)), "\n")
t0 <- proc.time()[["elapsed"]]
r <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                         reporter = "silent")
df <- as.data.frame(r)
cat(sprintf("RESULT %s: pass=%d fail=%d error=%d skip=%d warn=%d (%.0f s)\n",
            basename(f), sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning), proc.time()[["elapsed"]] - t0))
bad <- df[df$failed > 0 | df$error | df$warning > 0, "test"]
if (length(bad)) cat("BAD:", paste(bad, collapse = " | "), "\n")
