# One test file, one R process, on the lane's library stack.
# A copy of dev/release/run-tests.R with the library swapped: the lane
# build first, then the round's reference build for the extensions this
# lane did not change, then the user library.
#
#   Rscript dev/adefects-run.R <package> <path-to-test-file>
LIB <- c("C:/Users/adf44/source/r/adefects-lib",
         "C:/Users/adf44/source/r/rellib-r3",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(LIB)

# The StanHeaders pin is gone; anything that compiles a fresh Stan
# program needs the user Makevars flag instead, and R only reads
# <HOME>/.R/Makevars.win, where HOME depends on the launcher.
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}

# A guard that fails CLOSED. Without NOT_CRAN every skip_on_cran()
# block is skipped, and a file made entirely of them reports
# pass=0 fail=0 err=0 skip=0, which is indistinguishable from a file
# that ran and asserted nothing. Refusing beats reporting zero.
if (!nzchar(Sys.getenv("NOT_CRAN"))) {
  stop("NOT_CRAN is unset: every skip_on_cran() block would be skipped ",
       "and the file would report zero expectations. Export NOT_CRAN=true.")
}

suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
suppressMessages(library(p, character.only = TRUE))

r <- tryCatch(
  as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", p, " ", basename(f), " LOADERROR ", conditionMessage(e),
        "\n", sep = "")
    NULL
  })

if (!is.null(r)) {
  cat("RESULT ", p, " ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      "\n", sep = "")
  bad <- r[r$failed > 0 | r$error > 0, , drop = FALSE]
  for (i in seq_len(nrow(bad))) {
    cat("  BLOCK ", bad$test[i], " fail=", bad$failed[i], " err=",
        bad$error[i], "\n", sep = "")
  }
}
