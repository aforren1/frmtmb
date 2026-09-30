# One test file, this lane's library first (or the base build with
# --base), one line of counts and every failing expectation's text.
#   Rscript dev/postfit2-runtest.R <pkg> <file> [--base] [--gated]
# --gated sets the gated tiers' switches here, in R: exported from a
# driver under xargs they did not reach R, and a whole gated run
# skipped (dev/postfit2-gated-run1/).
args <- commandArgs(TRUE)
pkg <- args[1]
file <- args[2]
base_only <- "--base" %in% args
if ("--gated" %in% args) {
  Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true",
             FRMTMB_DRMTMB_FIT_TESTS = "true", FRMTMB_FUZZ = "true")
  # the same driver lost TMP, and a fresh Stan compile then wrote its
  # temporary files to C:\WINDOWS and failed
  if (!nzchar(Sys.getenv("TMP"))) {
    Sys.setenv(TMP = "C:/Users/adf44/AppData/Local/Temp",
               TEMP = "C:/Users/adf44/AppData/Local/Temp")
  }
}
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (base_only) libs <- libs[-1]
# --r0: the lane's round-0 build (the tarballs its final R CMD check
# built, installed into a private library of their own), the build the
# punch-round blockers were found on
if ("--r0" %in% args) {
  libs[1] <- "C:/Users/adf44/source/r/wt-postfit2-r0lib"
}
# --r1: the round-1 build (punch round 1's installed packages, copied
# before round 2 was installed), the build punch round 2 was found on
if ("--r1" %in% args) {
  libs[1] <- "C:/Users/adf44/source/r/wt-postfit2-r1lib"
}
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-postfit2/dev/stan-cache")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2"
dir <- if (pkg == "frmtmb") file.path(wt, "tests/testthat") else
  file.path(wt, "extensions", pkg, "tests/testthat")
suppressMessages(library(testthat))
# --attach: attach frmtmb and the package, as each package's own
# tests/testthat.R does; an extension's tests call frm() unqualified
if ("--attach" %in% args) {
  suppressMessages(library(frmtmb))
  suppressMessages(library(pkg, character.only = TRUE))
}
cat("libs:", paste(.libPaths(), collapse = " | "), "\n")
cat(pkg, "version", format(packageVersion(pkg)), "from",
    find.package(pkg), "\n")
r <- testthat::test_file(file.path(dir, file), package = pkg,
                         env = testthat::test_env(pkg),
                         reporter = testthat::ListReporter$new(),
                         stop_on_failure = FALSE)
df <- as.data.frame(r)
cat(sprintf("RESULT %s %s: tests=%d failed=%d error=%d skipped=%d warning=%d passed=%d\n",
            pkg, file, sum(df$nb), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning), sum(df$passed)))
for (res in r) {
  for (e in res$results) {
    if (inherits(e, c("expectation_failure", "expectation_error",
                      "expectation_warning"))) {
      cat("----", class(e)[1], "in", res$test, "\n")
      cat(substr(conditionMessage(e), 1, 1500), "\n")
    }
  }
}
sk <- unique(unlist(lapply(r, function(res) {
  vapply(Filter(function(e) inherits(e, "expectation_skip"), res$results),
         conditionMessage, "")
})))
if (length(sk)) cat("SKIP reasons:", paste(head(sk, 8), collapse = " | "), "\n")
