# Review of lane names: one ported brms-suite file in tier mode (verdicts
# asserted both ways), on the lane or base build, every extension loaded.
# Usage: Rscript dev/reviews/names-review-suite.R <base|lane> <pkg> <file>
a <- commandArgs(TRUE)
arm <- a[1]; p <- a[2]; f <- a[3]
lib <- c(if (arm == "lane") "C:/Users/adf44/source/r/wt-names-lib",
         "C:/Users/adf44/source/r/rellib-r20",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = p,
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
for (e in c("frmtmb.sample", "frmtmb.coupling", "frmtmb.eam",
            "frmtmb.latent", "frmtmb.learn", "frmtmb.ode", "frmtmb.spatial",
            "frmtmb.survival")) {
  suppressMessages(library(e, character.only = TRUE))
}
cat("lib:", p, format(packageVersion(p)), find.package(p), "| frmtmb",
    find.package("frmtmb"), "\n")
tdir <- if (p == "frmtmb") "tests/testthat" else
  file.path("extensions", p, "tests/testthat")
if (!file.exists(file.path(tdir, basename(f)))) {
  # a bin-2 file: staged beside the package's helpers
  stage <- file.path("dev/reviews/names-review-stage", arm, p,
                     sub("[.]R$", "", basename(f)))
  unlink(stage, recursive = TRUE)
  dir.create(stage, recursive = TRUE)
  file.copy(list.files(tdir, "^helper.*[.]R$", full.names = TRUE), stage)
  file.copy("dev/brmsport-bin2-adapters.R",
            file.path(stage, "helper-brms-bin2-adapters.R"))
  file.copy(file.path("dev/brmsport-bin2", p, basename(f)), stage)
  f <- file.path(stage, basename(f))
} else {
  f <- file.path(tdir, basename(f))
}
tr <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent", stop_on_failure = FALSE)
for (t in tr) {
  for (e in t$results) {
    if (inherits(e, c("expectation_failure", "expectation_error"))) {
      cat("FAILURE [", t$test, "]: ", substr(conditionMessage(e), 1, 600),
          "\n", sep = "")
    }
  }
}
r <- as.data.frame(tr)
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d\n", basename(f),
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped)))
