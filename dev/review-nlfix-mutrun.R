# One test file under one mutant. Rscript review-nlfix-mutrun.R <mutant> <file> <pkg>
a <- commandArgs(trailingOnly = TRUE)
Sys.setenv(REV_MUT = a[1])
lib <- Sys.getenv("REV_TEST_LIB", "C:/Users/adf44/source/r/rev-nlfix-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-nlfix/dev/stan-cache", NOT_CRAN = "true")
p <- a[3]
suppressPackageStartupMessages({ library(testthat); library(p, character.only = TRUE) })
r <- as.data.frame(test_file(a[2], package = p, env = test_env(p), reporter = "silent", stop_on_failure = FALSE))
cat(sprintf("MUT %-16s %-26s pass=%d fail=%d err=%d skip=%d lib=%s\n", a[1], basename(a[2]), sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped), dirname(find.package("frmtmb"))))
bad <- r[r$failed > 0 | r$error, ]
if (nrow(bad)) cat(paste0("   BAD: ", unique(bad$test)), sep = "\n")
