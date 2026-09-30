# Reviewer (lane ceplot): one test file with the summary reporter, which
# prints every failing expectation and its message.
#   Rscript dev/ceplot-rev-debug1.R <file> [package] [lane|base]
a <- commandArgs(trailingOnly = TRUE)
f <- a[1]
p <- if (length(a) >= 2) a[2] else "frmtmb"
arm <- if (length(a) >= 3) a[3] else "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib", libs)
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
library(testthat)
suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
cat("ARM", arm, "|", p, "from", find.package(p), "\n")
invisible(test_file(f, package = p, env = test_env(p),
                    reporter = SummaryReporter$new(max_reports = Inf)))
