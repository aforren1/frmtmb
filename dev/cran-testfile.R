# Lane cran: one core test file, one R process, as CRAN would run it
# (NOT_CRAN unset, so skip_on_cran() skips). Prints the whole process's
# CPU and elapsed time, the file's, a fixed control's, and one TEST line
# per test_that() block with its own CPU and elapsed time.
#
#   Rscript dev/cran-testfile.R <path-to-test-file>
#
# The control is the same arithmetic in every process: its spread across
# the run is the machine-load factor that every other time here carries
# (dev/lane-rules.md, "An examples-timing NOTE on this box measures LOAD").
local({
t_proc0 <- proc.time()
.libPaths(c("C:/Users/adf44/source/r/wt-cran-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(!nzchar(Sys.getenv("NOT_CRAN")))
ctl <- system.time({
  s <- 0
  for (i in seq_len(4e6)) s <- s + sqrt(i)
})
suppressMessages(library(testthat))
t_lib0 <- proc.time()
suppressMessages(library(frmtmb))
t_lib <- proc.time() - t_lib0
f <- commandArgs(trailingOnly = TRUE)[1]
cat("lib: frmtmb ", as.character(utils::packageVersion("frmtmb")), " ",
    dirname(find.package("frmtmb")), "\n", sep = "")
res <- NULL
t0 <- proc.time()
r <- tryCatch({
  res <- test_file(f, package = "frmtmb",
                   env = testthat::test_env("frmtmb"), reporter = "silent")
  as.data.frame(res)
}, error = function(e) {
  cat("LOADERROR ", conditionMessage(e), "\n", sep = "")
  NULL
})
t1 <- proc.time() - t0
if (!is.null(res)) {
  for (t in res) {
    for (x in t$results) {
      k <- if (inherits(x, "expectation_failure")) "FAIL" else
        if (inherits(x, "expectation_error")) "ERROR" else
          if (inherits(x, "expectation_skip")) "SKIP" else
            if (inherits(x, "expectation_warning")) "WARN" else NA
      if (!is.na(k)) {
        cat("DETAIL ", k, " [", t$test, "] ",
            gsub("\n", " | ", conditionMessage(x)), "\n", sep = "")
      }
    }
  }
}
if (!is.null(r)) {
  for (j in seq_len(nrow(r))) {
    cat("TEST\t", basename(f), "\t", gsub("[\t\n]", " ", r$test[j]), "\t",
        r$nb[j], "\t", r$skipped[j], "\t",
        sprintf("%.2f", r$user[j] + r$system[j]), "\t",
        sprintf("%.2f", r$real[j]), "\n", sep = "")
  }
}
tp <- proc.time() - t_proc0
cpu <- function(t) sum(t[c(1, 2)])
cat("RESULT\t", basename(f),
    "\tpass=", if (is.null(r)) NA else sum(r$passed),
    "\tfail=", if (is.null(r)) NA else sum(r$failed),
    "\terr=", if (is.null(r)) NA else sum(r$error),
    "\tskipblocks=", if (is.null(r)) NA else sum(r$skipped),
    "\tblocks=", if (is.null(r)) NA else nrow(r),
    "\twarn=", if (is.null(r)) NA else sum(r$warning),
    "\tfile_cpu=", sprintf("%.2f", cpu(t1)),
    "\tfile_real=", sprintf("%.2f", t1[[3]]),
    "\tproc_cpu=", sprintf("%.2f", cpu(tp)),
    "\tproc_real=", sprintf("%.2f", tp[[3]]),
    "\tlib_real=", sprintf("%.2f", t_lib[[3]]),
    "\tctl_real=", sprintf("%.2f", ctl[[3]]), "\n", sep = "")
})
