# The new tests against the bugs they are for, one mutant per process.
# Usage: Rscript dev/refuse-f2-mutants.R <name>. Reports the failures.
name <- commandArgs(trailingOnly = TRUE)[1L]
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
ns <- asNamespace("frmtmb")
sub_fn <- function(fn, from, to) {
  txt <- deparse(get(fn, ns))
  i <- grep(from, txt, fixed = TRUE)
  stopifnot(length(i) >= 1L)
  txt[i] <- gsub(from, to, txt[i], fixed = TRUE)
  f <- eval(parse(text = txt)); environment(f) <- ns
  assignInNamespace(fn, f, "frmtmb")
}
file <- "tests/testthat/test-mi-rows-read.R"
switch(name,
  # one pass of the peel: the old classification
  one_pass = sub_fn("mi_rows_read", "identical(lengths(read), before)", "TRUE"),
  # a weight of 0 reads
  weight_reads = sub_fn("mi_rows_live", "rep_len(w > 0, n)", "rep(TRUE, n)"),
  # a live missing row counts however it has left
  gone_reads = sub_fn("mi_rows_live", "!is.null(m)", "FALSE"),
  # never warns
  no_warn = sub_fn("mi_laplace_note", "abs(total) > mi_laplace_limit", "FALSE"),
  # empty miss reaches the chart
  empty_miss = sub_fn("mi_values", "length(mm_$rows)", "TRUE"),
  # a reader reads the placeholder at a row that left
  placeholder = sub_fn("mi_values", "length(mm_$unread)", "FALSE"),
  stop("no such mutant"))
res <- as.data.frame(test_file(file, reporter = "silent", stop_on_failure = FALSE,
                               env = test_env("frmtmb"), package = "frmtmb"))
cat(sprintf("MUT %-14s pass=%d fail=%d error=%d: %s\n", name, sum(res$passed),
            sum(res$failed), sum(res$error),
            paste(unique(substr(res$test[res$failed > 0L | res$error], 1, 50)),
                  collapse = " | ")))
