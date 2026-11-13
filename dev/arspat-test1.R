# Lane arspat: run ONE test file of one package on the lane build.
#   Rscript dev/arspat-test1.R <package> <test file path>
# ARSPAT_LIB="" runs it on the base build (rellib-r20) alone.
args <- commandArgs(trailingOnly = TRUE)
p <- args[1]
f <- args[2]
lib <- Sys.getenv("ARSPAT_LIB", "C:/Users/adf44/source/r/wt-arspat-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-arspat/dev/stan-cache")
}
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", dirname(find.package("frmtmb")), "|",
    p, "from", dirname(find.package(p)), "\n")
r <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                         reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(r)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning)))
