# Run one test file (optionally one test_that() block by its description)
# against the lane build, and print a RESULT line with the counts.
#   Rscript dev/fams2-run-test.R <pkg> <file> [desc]
# <pkg> is frmtmb or frmtmb.sample; <file> is relative to its tests dir.
args <- commandArgs(TRUE)
pkg <- args[1]
file <- args[2]
desc <- if (length(args) >= 3) args[3] else NULL
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (identical(Sys.getenv("FAMS2_BASE"), "true")) {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r3",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
}
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/stan-cache")
options(frmtmb.brms_lp_report = TRUE)
root <- "C:/Users/adf44/source/r/frmtmb-wt-fams2"
tdir <- if (pkg == "frmtmb") file.path(root, "tests/testthat") else
  file.path(root, "extensions", pkg, "tests/testthat")
suppressPackageStartupMessages(library(testthat))
cat("lib:", find.package(pkg), "\n")
res <- testthat::test_file(file.path(tdir, file), package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = "summary", desc = desc,
                           load_package = "installed")
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s%s pass=%d fail=%d err=%d skip=%d warn=%d blocks=%d\n",
            pkg, file, if (is.null(desc)) "" else paste0(" [", desc, "]"),
            sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning), nrow(df)))
