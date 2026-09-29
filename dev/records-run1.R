LIB <- "C:/Users/adf44/source/r/wt-records-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
suppressMessages(library(p, character.only = TRUE))
r <- tryCatch(
  as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })
if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      "\n", sep = "")
  bad <- r[r$failed > 0 | r$error, , drop = FALSE]
  if (nrow(bad)) print(bad[, c("test", "failed", "error")])
}
