## One test file, one R process: the standing rule for a correctness run.
##
## Usage: Rscript dev/famlink-run1.R <pkgdir> <filter> [base]
##   pkgdir  "." for core, or extensions/frmtmb.sample
##   filter  testthat filter, i.e. the file name without test- and .R
##   base    any third argument runs against the SHARED REFERENCE
##           library alone, the base commit aa9227e, so a test can be seen
##           failing before the change is installed.
args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]
filt <- args[2]
base <- length(args) >= 3L
own <- if (base) character() else "C:/Users/adf44/source/r/famlink-lib"
.libPaths(c(own,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n")
cat("version    :", as.character(utils::packageVersion("frmtmb")), "\n")
if (!identical(pkgdir, ".")) {
  pkg <- basename(pkgdir)
  suppressMessages(library(pkg, character.only = TRUE))
  cat(pkg, "from:", dirname(system.file(package = pkg)), "\n")
}
this_pkg <- if (identical(pkgdir, ".")) "frmtmb" else basename(pkgdir)
setwd(file.path(pkgdir, "tests", "testthat"))
t0 <- Sys.time()
# `package =` is what puts the package NAMESPACE on the test
# environment's search path. Without it every test that calls an
# unexported helper by its bare name errors with "could not find
# function", and a runner that only sums `failed` prints a green line
# for the file anyway.
res <- testthat::test_dir(".", filter = filt, package = this_pkg,
                          reporter = testthat::SummaryReporter$new(max_reports = Inf),
                          stop_on_failure = FALSE)
d <- as.data.frame(res)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("\nFAMLINK %s pass %d fail %d error %d skip %d in %.0f s\n",
            filt, sum(d$passed), sum(d$failed), sum(d$error),
            sum(d$skipped), secs))
