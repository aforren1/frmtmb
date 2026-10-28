# The tests that build their own library stack, run under covr as
# .github/workflows/test-coverage.yaml calls it (install_path set, so
# covr's R_LIBS_SITE carries the whole stack): lane-env.md, last
# section. covr installs the package from the worktree into
# install_path and runs the one test file there.
#   Rscript dev/rel077-covr.R <package dir> <test file name>
# Log: dev/rel077-log/covr-<package>.txt. Adapted from dev/rel076-covr.R.
a <- commandArgs(trailingOnly = TRUE)
lib <- "C:/Users/adf44/source/r/rellib-r15"
.libPaths(c(lib,
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
pkg <- read.dcf(file.path(a[1], "DESCRIPTION"), "Package")[1, 1]
ip <- file.path(tempdir(), "covr-install")
code <- sprintf(paste0(
  "library(testthat); library(%s); ",
  "r <- as.data.frame(testthat::test_file(file.path(%s, 'tests', ",
  "'testthat', %s), package = %s, reporter = 'silent', ",
  "env = testthat::test_env(%s))); ",
  "cat('COVR RESULT', %s, 'pass', sum(r$passed), ",
  "'fail', sum(r$failed), 'err', sum(r$error), ",
  "'skip', sum(r$skipped), '\\n')"),
  pkg, deparse(normalizePath(a[1], winslash = "/")), deparse(a[2]),
  deparse(pkg), deparse(pkg), deparse(a[2]))
cov <- covr::package_coverage(a[1], type = "none", code = code,
                              quiet = FALSE, clean = FALSE,
                              install_path = ip)
# covr runs `code` with R CMD BATCH, R_LIBS_SITE set to its whole
# stack; the output is the .Rout it leaves in install_path
for (f in list.files(ip, "[.]Rout$", recursive = TRUE, full.names = TRUE)) {
  x <- readLines(f, warn = FALSE)
  cat("== ", f, "\n", sep = "")
  cat(grep("COVR RESULT|Error|FAIL", x, value = TRUE), sep = "\n")
}
print(cov)
