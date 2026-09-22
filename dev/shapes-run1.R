# One test file with the FULL failure text, for triage.
#   Rscript dev/shapes-run1.R <package> <path-to-test-file>
.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
suppressMessages(library(p, character.only = TRUE))
# the summary reporter prints every failure's own text, which the
# silent one keeps in an object this testthat version does not hand back
invisible(test_file(a[2], package = p, env = testthat::test_env(p),
                    reporter = "summary"))
