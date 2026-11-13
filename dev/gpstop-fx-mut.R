# Lane gpstop, fix round: test-gp-start.R on the lane build with a mutant,
# to see the new tests fail where the fix is undone. Args: mutant
#   regular   gp_relin_regular() returns TRUE (the build before the fix)
#   template  relin also on a fitted template (refit() runs it)
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-gpstop-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/stan-cache")
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
put <- function(name, f) {
  environment(f) <- ns
  unlockBinding(name, ns); assign(name, f, envir = ns); lockBinding(name, ns)
}
if (a[1] == "regular") put("gp_relin_regular", function(...) TRUE)
if (a[1] == "template") {
  f <- get("fit_assembled", ns)
  b <- deparse(body(f))
  i <- grep("relin = gp_body && announce_start", b, fixed = TRUE)
  stopifnot(length(i) == 1L,
            grepl("(gp_cold || prefit_seeded)", b[i + 1L], fixed = TRUE))
  b[i + 1L] <- sub("(gp_cold || prefit_seeded)", "TRUE", b[i + 1L],
                   fixed = TRUE)
  b[i] <- sub("announce_start && ", "", b[i], fixed = TRUE)
  body(f) <- parse(text = b)[[1]]
  put("fit_assembled", f)
}
cat("MUTANT", a[1], "\n")
r <- as.data.frame(testthat::test_file(
  "C:/Users/adf44/source/r/frmtmb-wt-gpstop/tests/testthat/test-gp-start.R",
  package = "frmtmb", env = testthat::test_env("frmtmb"),
  reporter = "summary"))
cat(sprintf("RESULT mutant %s pass=%d fail=%d err=%d skip=%d warn=%d\n", a[1],
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
            sum(r$warning)))
