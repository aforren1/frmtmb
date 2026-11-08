# Review of lane loglik, short check: does a binding uniform(0, 0.2) on a
# slope with no lb/ub fail on the base build too? Arg 1: before | after.
# Seed 703, the data of review-loglik-short.R section 2.
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "after"
libs <- c("C:/Users/adf44/source/r/wt-loglik-lib",
          "C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "before") libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
cat("arm:", arm, " lib:", dirname(find.package("frmtmb")), "\n")
set.seed(703)
d3 <- data.frame(x = rnorm(120), g = factor(rep(1:20, each = 6)))
d3$y <- 1 + 0.5 * d3$x + rnorm(20, 0, 0.7)[d3$g] + rnorm(120)
r <- tryCatch(suppressWarnings(frm(y ~ x + (1 | g), data = d3,
              prior = set_prior("uniform(0, 0.2)", class = "b", coef = "x"))),
              error = function(e) conditionMessage(e))
cat(if (is.character(r)) paste("ERROR", substr(r, 1, 160)) else
  paste("b_x", r$estimates$beta[["x"]]), "\n")
