# Run this package's test suite one FILE at a time, resumably.
#
#   Rscript dev/run-suite.R [ledger] [filter]
#
# This script exists for the reason dev/fuzz/run-fuzz.R exists in
# frmtmb: a full run of this suite is minutes of NUTS, and testthat's
# reporter swallows progress on a long job. It exists RESUMABLY for a
# reason of its own: every file here samples, so a run that is
# interrupted - a timeout, a crash, a laptop lid - has cost real
# sampling time, and repeating the files that already passed is the
# expensive way to find out nothing changed.
#
# The ledger is a tab-separated append-only log, one line per file.
# A file recorded PASS is skipped on the next run. Delete the ledger to
# force a full run; pass a filter (a regular expression over file
# names) to run a subset.
#
# Environment, as in frmtmb (see its CONTRIBUTING.md):
#   NOT_CRAN=true              the heavy reference tier; set here
#   FRMTMB_SAMPLER_GATES=false turns the chain-agreement gates off, for
#                              a machine whose seeded Stan chains do
#                              not match this one's. Structural and
#                              exactness assertions never take it.
#   FRMTMB_BRMS_FIT_TESTS=true lets the brms-oracle blocks compile Stan

args <- commandArgs(trailingOnly = TRUE)
ledger <- if (length(args) >= 1L) args[[1L]] else "dev/suite-ledger.tsv"
filt <- if (length(args) >= 2L) args[[2L]] else NULL

if (!dir.exists("R") || !file.exists("DESCRIPTION")) {
  stop("run this from the package root: Rscript dev/run-suite.R")
}

Sys.setenv(NOT_CRAN = "true")
suppressMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.sample)
})

tdir <- file.path("tests", "testthat")
files <- sort(basename(list.files(tdir, pattern = "^test-.*[.]R$")))
if (!is.null(filt)) files <- grep(filt, files, value = TRUE)
if (!length(files)) stop("no test files matched")

passed <- character(0)
if (file.exists(ledger)) {
  for (l in readLines(ledger, warn = FALSE)) {
    p <- strsplit(l, "\t", fixed = TRUE)[[1L]]
    if (length(p) >= 2L && identical(p[[2L]], "PASS")) {
      passed <- c(passed, p[[1L]])
    }
  }
}

ledger_dir <- dirname(ledger)
if (nzchar(ledger_dir) && !dir.exists(ledger_dir)) {
  dir.create(ledger_dir, recursive = TRUE)
}

owd <- setwd(tdir)
on.exit(setwd(owd), add = TRUE)
rel <- file.path(owd, ledger)

bad <- character(0)
for (f in files) {
  if (f %in% passed) {
    cat(sprintf("%-34s  skipped (ledger says PASS)\n", f))
    next
  }
  t0 <- Sys.time()
  res <- tryCatch(
    as.data.frame(test_file(f, reporter = "silent",
                            package = "frmtmb.sample")),
    error = function(e) e)
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (inherits(res, "error")) {
    status <- "ERROR"
    detail <- conditionMessage(res)
  } else {
    nbad <- sum(res$failed) + sum(res$error)
    status <- if (nbad > 0L) "FAIL" else "PASS"
    detail <- sprintf("pass=%d fail=%d skip=%d", sum(res$passed), nbad,
                      sum(res$skipped))
  }
  if (!identical(status, "PASS")) bad <- c(bad, f)
  cat(sprintf("%-34s  %-5s %7.1fs  %s\n", f, status, el, detail))
  cat(paste(f, status, el, detail, sep = "\t"), "\n", sep = "",
      file = rel, append = TRUE)
}

cat("\nledger: ", ledger, "\n", sep = "")
if (length(bad)) {
  cat("not passing: ", paste(bad, collapse = ", "), "\n", sep = "")
  quit(status = 1L)
}
cat("all files passed\n")
