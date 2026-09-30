# Lane sampfix: run ONE test file against one build and print its counts.
#
#   Rscript dev/sampfix-runtest.R <lane|ref> <frmtmb|frmtmb.sample> <file>
#
# lane = the private library first, then rellib-r3; ref = rellib-r3.
# The file is read from this worktree in both arms, so a ref run is the
# new test against the unfixed code.

a <- commandArgs(trailingOnly = TRUE)
arm <- a[1L]
pkg <- a[2L]
file <- a[3L]
LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
WT <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix"
dir <- if (pkg == "frmtmb") file.path(WT, "tests/testthat") else
  file.path(WT, "extensions", pkg, "tests/testthat")
suppressMessages(library(testthat))
suppressMessages(library(pkg, character.only = TRUE))
cat("ARM ", arm, " ", pkg, " from ", dirname(find.package(pkg)), "\n",
    sep = "")
res <- testthat::test_file(file.path(dir, file), package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = testthat::ProgressReporter$new(
                             show_praise = FALSE),
                           stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s: tests %d pass %d fail %d error %d skip %d warn %d\n",
            arm, file, nrow(df), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
