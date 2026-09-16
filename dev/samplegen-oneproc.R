# The whole frmtmb.sample suite in ONE process, the way R CMD check runs
# it: an earlier file's loaded namespaces are still there for a later
# one. frmtmb's lane found a defect only this run could expose.
#   Rscript dev/samplegen-oneproc.R
.libPaths(c("C:/Users/adf44/source/r/samplegen-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages({ library(testthat); library(frmtmb.sample) })
message("frmtmb.sample from: ", find.package("frmtmb.sample"))
res <- testthat::test_dir(
  "C:/Users/adf44/source/r/frmtmb-wt-samplegen/extensions/frmtmb.sample/tests/testthat",
  reporter = "silent", package = "frmtmb.sample", stop_on_failure = FALSE,
  load_package = "none")
df <- as.data.frame(res)
tot <- function(nm) sum(df[[nm]], na.rm = TRUE)
cat("FILES  ", length(unique(df$file)), "\n")
cat("BLOCKS ", nrow(df), "\n")
cat("ASSERT ", tot("nb"), "\n")
cat("FAIL   ", tot("failed"), "\n")
cat("ERROR  ", tot("error"), "\n")
cat("SKIP   ", tot("skipped"), "\n")
cat("brms loaded at the end:", isNamespaceLoaded("brms"), "\n")
for (i in which(df$failed > 0 | df$error)) cat("BAD    ", df$file[i], ":", df$test[i], "\n")
