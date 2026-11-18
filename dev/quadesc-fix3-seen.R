# Lane quadesc, fix round 1: the test of finding 5 (an abandoned restart
# carries evals and nonfinite_trials) seen failing. The lane's
# escape_stationary() is parsed from R/fit.R without the two lines that
# set the counts and assigned into the installed lane namespace; then
# test-quad-escape-record.R runs.
#   Rscript dev/quadesc-fix3-seen.R
.libPaths(c("C:/Users/adf44/source/r/wt-quadesc-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
src <- readLines("R/fit.R")
i0 <- grep("^escape_stationary <- function", src)
i1 <- i0 + which(src[(i0 + 1):length(src)] == "}")[1]
body <- src[i0:i1]
drop <- grep("b[$]evals <- tl[$]evals|b[$]nonfinite_trials <- tl[$]n",
             body)
cat("lines dropped:", length(drop), "\n")
body <- body[-drop]
# the if () { } around the trials line is left empty
eval(parse(text = body))
environment(escape_stationary) <- ns
unlockBinding("escape_stationary", ns)
assign("escape_stationary", escape_stationary, envir = ns)
r <- as.data.frame(test_file("tests/testthat/test-quad-escape-record.R",
                             package = "frmtmb",
                             env = testthat::test_env("frmtmb"),
                             reporter = "silent"))
print(r[, c("test", "nb", "failed", "error")])
