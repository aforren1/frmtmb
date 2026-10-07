# The ciharden review's script on the 0.69.0 release library (rellib-r7)
# and the release tree's test-perf.R (the merged-tree rerun it asks).
# Reviewer, punch round 1: does test-perf.R's guard test fail when its
# instrument is broken, and skip (not pass) where memory profiling is
# unavailable? Runs text variants of the lane's test-perf.R.
#   Rscript dev/rel069-perfguard.R <variant>
# variants: asis, noprofmem, bytesNA, bytesZero, noloop
.libPaths(c("C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
v <- commandArgs(TRUE)[1]
src <- readLines(
  "C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat/test-perf.R")
sub1 <- function(pat, rep) {
  hit <- grepl(pat, src, fixed = TRUE)
  if (!any(hit)) stop("pattern not found: ", pat)
  src <<- gsub(pat, rep, src, fixed = TRUE)
}
switch(v,
  asis = NULL,
  # R built without memory profiling
  noprofmem = sub1('isTRUE(capabilities("profmem"))', "FALSE"),
  # the bytes half of the census broken: it reports nothing
  bytesNA = sub1("c(nodes = nrow(RTMB::GetTape(obj)$data.frame()), bytes = by)",
                 "c(nodes = nrow(RTMB::GetTape(obj)$data.frame()), bytes = NA_real_)"),
  # broken the other way: a constant
  bytesZero = sub1("c(nodes = nrow(RTMB::GetTape(obj)$data.frame()), bytes = by)",
                   "c(nodes = nrow(RTMB::GetTape(obj)$data.frame()), bytes = 1)"),
  # the guard's mutation absent: the mock adds no loop
  noloop = sub1("for (i in seq_len(n)) m[i] <- b * 2", "NULL"),
  stop("variant"))
cat("capabilities(profmem) here:", capabilities("profmem"), "\n")
tf <- tempfile("test-perf-", fileext = ".R")
writeLines(src, tf)
r <- test_file(tf, package = "frmtmb", env = test_env("frmtmb"),
               reporter = "silent")
for (t in r) {
  k <- vapply(t$results, function(x) class(x)[1], "")
  msg <- vapply(t$results, function(x) {
    if (inherits(x, "expectation_success")) "" else
      substr(gsub("\n", " ", conditionMessage(x)), 1, 110)
  }, "")
  cat(sprintf("%-8s | %-55s | %s\n", v, t$test,
              paste(sub("expectation_", "", k), collapse = ",")))
  for (m in msg[nzchar(msg)]) cat("           ", m, "\n")
}
