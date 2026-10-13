# Lane optrap: run one test file against a library, one process, and
# print one line per test (no reporter cap on failures).
#   Rscript dev/optrap-testfile.R <lib or ""> <pkg> <file>
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkg <- a[2]
suppressMessages({library(testthat); library(pkg, character.only = TRUE)})
# OPTRAP_ROUND0=1: the round-0 unit rule patched in (punch round 1, B1)
if (nzchar(Sys.getenv("OPTRAP_ROUND0"))) source("dev/optrap-round0-unit.R")
# OPTRAP_ROUND1=1: punch round 1's unit spread and Newton steps patched
# in (punch round 2)
if (nzchar(Sys.getenv("OPTRAP_ROUND1"))) source("dev/optrap-p2-round1.R")
# OPTRAP_ROUND2=1: punch round 2's newton_polish() patched in (the
# consolidation fixes)
if (nzchar(Sys.getenv("OPTRAP_ROUND2"))) source("dev/optrap-c3-round2.R")
# OPTRAP_GATE=<G>: newton_polish()'s gain gate at G instead of 0.1 (the
# consolidation fixes; 0.25 was the integration's candidate)
if (nzchar(Sys.getenv("OPTRAP_GATE"))) local({
  ns <- asNamespace("frmtmb")
  b <- deparse(get("newton_polish", ns))
  i <- grep("0.1 * gain", b, fixed = TRUE)
  stopifnot(length(i) == 1L)
  b[i] <- sub("0.1 * gain", paste0(Sys.getenv("OPTRAP_GATE"), " * gain"),
              b[i], fixed = TRUE)
  f <- eval(parse(text = b))
  environment(f) <- ns
  utils::assignInNamespace("newton_polish", f, "frmtmb")
  cat("newton_polish: gain gate", Sys.getenv("OPTRAP_GATE"), "\n")
})
# OPTRAP_NOPOLISH=1: newton_polish() taken out (punch round 1, m10)
if (nzchar(Sys.getenv("OPTRAP_NOPOLISH"))) {
  utils::assignInNamespace("newton_polish", function(obj, opt, ...) opt,
                           "frmtmb")
  cat("newton_polish: taken out\n")
}
cat("lib:", pkg, dirname(find.package(pkg)), "frmtmb",
    dirname(find.package("frmtmb")), "\n")
res <- testthat::test_file(a[3], package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = "silent")
r <- as.data.frame(res)
for (i in seq_len(nrow(r))) {
  cat(sprintf("TEST %s | nb=%d failed=%d error=%s skipped=%s warning=%d\n",
              r$test[i], r$nb[i], r$failed[i], r$error[i], r$skipped[i],
              r$warning[i]))
}
for (t in res) for (x in t$results) {
  if (inherits(x, c("expectation_failure", "expectation_error"))) {
    cat("  FAIL [", t$test, "]", substr(conditionMessage(x), 1, 200), "\n")
  }
  if (inherits(x, "expectation_warning")) {
    cat("  WARN [", t$test, "]", substr(conditionMessage(x), 1, 120), "\n")
  }
}
cat(sprintf("RESULT pass=%d fail=%d error=%d skip=%d warn=%d\n",
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
            sum(r$warning)))
