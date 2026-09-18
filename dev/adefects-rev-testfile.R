LIB <- Sys.getenv("REV_LIB")
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE")))
  Sys.setenv(FRMTMB_STAN_CACHE =
    "C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/stan-cache")
suppressPackageStartupMessages({ library(testthat) })
a <- commandArgs(trailingOnly = TRUE)
pkg <- a[1]; f <- a[2]
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
r <- as.data.frame(test_file(f, reporter = "silent",
                             package = pkg)) 
cat(sprintf("RESULT %s %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            pkg, basename(f), sum(r$passed), sum(r$failed),
            sum(r$error), sum(r$skipped), sum(r$warning)))
for (i in seq_len(nrow(r))) {
  if (r$failed[i] > 0 || r$error[i]) {
    cat(sprintf("  BLOCK fail=%-3d err=%-5s %s\n", r$failed[i],
                r$error[i], r$test[i]))
  }
}
