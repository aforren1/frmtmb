# Run one test file of one package against a given library stack and
# print its counts, so a lane's "before" (rellib-r18) and "after" (its
# own library first) arms differ only in the libraries.
#
#   Rscript dev/ci-runfile.R <test file> <pkg> <lib1;lib2;...>
args <- commandArgs(trailingOnly = TRUE)
f <- normalizePath(args[[1L]], winslash = "/")
p <- args[[2L]]
.libPaths(c(strsplit(args[[3L]], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", p, as.character(packageVersion(p)), find.package(p), "\n")
# La_library() and the BLAS entry are empty on Windows, so R.home() is
# what tells an OpenBLAS copy of R from the installed one
cat("R.home:", R.home(), "\n")
cat("La_library:", La_library(), "\n")
cat("BLAS:", extSoftVersion()[["BLAS"]], "\n")
res <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                           reporter = "check", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warning=%d\n",
            basename(f), sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning)))
