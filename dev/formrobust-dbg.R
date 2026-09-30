# Run one test file with a verbose reporter, for debugging failures.
# Usage: Rscript dev/formrobust-dbg.R <pkg> <file>
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(a[1], character.only = TRUE))
invisible(test_file(a[2], package = a[1], env = test_env(a[1]),
                    reporter = ProgressReporter$new(show_praise = FALSE)))
