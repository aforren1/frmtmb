## The whole suite in ONE process, which is what the CI workflow's
## comment counts and what R CMD check runs. The standing one-file-per
## -process rule is about isolating a correctness run; this is a
## measurement of the number that comment carries.
##
## Run: Rscript dev/coh-suite-onefile.R
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "false")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
setwd("extensions/frmtmb.coupling/tests/testthat")
t0 <- Sys.time()
res <- testthat::test_dir(".", reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(res)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
per <- stats::aggregate(cbind(passed, failed, error, skipped) ~ file, d,
                        sum)
print(per)
cat(sprintf("\nTOTAL pass %d fail %d error %d skip %d in %.0f s\n",
            sum(d$passed), sum(d$failed), sum(d$error), sum(d$skipped),
            secs))
cat(sprintf(paste0("less test-message-uniqueness.R, which skips itself",
                   " inside R CMD check: %d\n"),
            sum(d$passed) - sum(per$passed[per$file ==
                                             "test-message-uniqueness.R"])))
