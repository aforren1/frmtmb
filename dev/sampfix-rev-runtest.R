# Reviewer, lane sampfix: run ONE test file against one build, print its
# counts and every warning that escaped (dev/warnleak-scan.R's rule,
# with the library chosen by the arm).
#
#   Rscript dev/sampfix-rev-runtest.R <lane|ref> <pkg> <file>
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1L]; pkg <- a[2L]; file <- a[3L]
LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-rev-stan-cache")
WT <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix"
dir <- if (pkg == "frmtmb") file.path(WT, "tests/testthat") else
  file.path(WT, "extensions", pkg, "tests/testthat")
suppressMessages(library(testthat))
suppressMessages(library(pkg, character.only = TRUE))
cat("ARM ", arm, " ", pkg, " from ", dirname(find.package(pkg)), "\n",
    sep = "")
r <- testthat::test_file(file.path(dir, file), package = pkg,
                         env = testthat::test_env(pkg),
                         reporter = "list", stop_on_failure = FALSE)
df <- as.data.frame(r)
cat(sprintf("RESULT %s %s: tests %d pass %d fail %d error %d skip %d warn %d\n",
            arm, file, nrow(df), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
for (x in r) for (e in x$results) {
  cls <- class(e)[1L]
  if (cls %in% c("expectation_success", "expectation_skip")) next
  sr <- e$srcref
  loc <- if (is.null(sr)) "?" else paste0(basename(file), ":", sr[1])
  msg <- gsub("[\r\n\t]+", " ", conditionMessage(e))
  cat(sprintf("%s\t%s\t%s\t%s\n", toupper(sub("expectation_", "", cls)),
              loc, x$test, substr(msg, 1, 300)))
}
