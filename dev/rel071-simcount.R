# Count the statistical assertions of tests/testthat/test-simulate-density.R
# on a build, for the count its header states (SIM_ALPHA's trade). Each
# call of expect_lt() is classified by its label, which the file's
# helpers set: "chi-square (df", " mean", " variance", "P(Y = " and the
# logistic-normal correlation check ("lncor").
#
#   Rscript dev/rel071-simcount.R <lib> [test file] > dev/rel071-log/simcount.txt
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
labs <- character(0)
tr <- quote(labs <<- c(labs, if (is.null(label)) NA_character_ else label))
# the helpers call testthat::expect_lt(); the two correlation checks
# call a bare expect_lt(), which resolves through the attached
# package:testthat and not the namespace (the release review, 6c), so
# both bindings are traced
trace("expect_lt", where = asNamespace("testthat"), print = FALSE,
      tracer = tr)
trace("expect_lt", where = as.environment("package:testthat"),
      print = FALSE, tracer = tr)
f <- commandArgs(trailingOnly = TRUE)[2]
if (is.na(f)) f <- "C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat/test-simulate-density.R"
res <- as.data.frame(test_file(f, package = "frmtmb",
                               env = test_env("frmtmb"), reporter = "silent"))
untrace("expect_lt", where = asNamespace("testthat"))
untrace("expect_lt", where = as.environment("package:testthat"))
cat("file: pass", sum(res$nb) - sum(res$failed), "fail", sum(res$failed),
    "error", sum(res$error), "skip", sum(res$skipped), "\n")
cls <- ifelse(grepl("chi-square \\(df", labs), "chi-square",
       ifelse(grepl(" mean$", labs), "mean",
       ifelse(grepl(" variance$", labs), "variance",
       ifelse(grepl("P\\(Y = ", labs), "P(Y = k)",
       ifelse(grepl("lncor$", labs), "correlation", "other")))))
print(table(cls, useNA = "ifany"))
cat("statistical (all but other):", sum(cls != "other"), "\n")
a <- 1e-6
n <- sum(cls != "other")
cat(sprintf("P(any false alarm) at 1e-6: %.4f, at 1e-3: %.4f\n",
            1 - (1 - a)^n, 1 - (1 - 1e-3)^n))
if (any(cls == "other")) print(unique(labs[cls == "other"]))
