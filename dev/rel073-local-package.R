# Brief item 6: how long withr::local_package() keeps a package attached
# when a test file calls it at its top level. Two scratch directories of
# two files each (a attaches emmeans, b reports whether it is still
# attached), run with test_dir() in one session, as R CMD check and
# covr run a suite: once with the default environment (the file's own),
# once with teardown_env().
#   Rscript dev/rel073-local-package.R   (log dev/rel073-log/local-package.txt)
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(testthat))
cat("testthat", format(packageVersion("testthat")), "withr",
    format(packageVersion("withr")), "\n")
d <- tempfile("lp")
for (arm in c("file", "teardown")) {
  td <- file.path(d, arm)
  dir.create(td, recursive = TRUE)
  writeLines(c(if (arm == "file") 'withr::local_package("emmeans")' else
               'withr::local_package("emmeans", .local_envir = teardown_env())',
               'test_that("a", expect_true("package:emmeans" %in% search()))'),
             file.path(td, "test-a.R"))
  writeLines(c('test_that("b", {',
               '  cat("in the next file, emmeans attached:",',
               '      "package:emmeans" %in% search(), "\\n")',
               '  expect_true(TRUE)',
               '})'), file.path(td, "test-b.R"))
  cat("== arm", arm, "\n")
  test_dir(td, reporter = "silent", stop_on_failure = FALSE)
  cat("after the run, emmeans attached:", "package:emmeans" %in% search(),
      "\n")
  if ("package:emmeans" %in% search()) detach("package:emmeans")
}
cat("END\n")
