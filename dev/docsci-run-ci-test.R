# The one test file this lane's change can reach: it is the only test
# that reads .github/workflows. One file, one process.
pkgload::load_all(".", quiet = TRUE)
res <- testthat::test_file(
  "tests/testthat/test-ci-siblings.R",
  package = "frmtmb",
  reporter = testthat::MultiReporter$new(list(
    testthat::ProgressReporter$new(),
    testthat::CheckReporter$new()
  ))
)
df <- as.data.frame(res)
cat("PASS=", sum(df$passed), " FAIL=", sum(df$failed),
    " ERROR=", sum(df$error), " SKIP=", sum(df$skipped), "\n", sep = "")
if (any(df$skipped)) {
  cat("skipped: ", paste(df$test[df$skipped], collapse = "; "), "\n",
      sep = "")
}
