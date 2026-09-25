# Lane wt-mvprior: one test file with the check reporter, to read what a
# failure says. Usage: MVPRIOR_ARM=base|lane Rscript dev/mvprior-onefile.R <pkg> <file>
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
suppressMessages(library(testthat))
testthat::set_max_fails(Inf)
a <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(a[1], character.only = TRUE))
cat("frmtmb from", find.package("frmtmb"), "\n")
invisible(test_file(a[2], package = a[1], env = testthat::test_env(a[1]),
                    reporter = "progress"))
