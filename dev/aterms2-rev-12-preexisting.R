# Reviewer: are two failures seen on the lane build older than the lane?
# emmeans() on an offset(log(time)) poisson model, and an addition term
# given an expression (weights(wt * 2)). Per arm:
#   Rscript dev/aterms2-rev-12-preexisting.R <base|lane>
# Log: dev/aterms2-rev-log-12-preexisting.txt
arm <- commandArgs(TRUE)[[1]]
libs <- c("C:/Users/adf44/source/r/wt-aterms2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
set.seed(606)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4), wt = runif(n, 0.5, 2))
d$y <- rpois(n, exp(0.2 + 0.5 * d$x) * d$time)
r <- function(e) tryCatch({suppressWarnings(suppressMessages(e)); "ok"},
                          error = function(e) conditionMessage(e))
fo <- suppressWarnings(frm(y ~ x + offset(log(time)), data = d, family = poisson()))
cat(arm, "emmeans offset:", r(summary(emmeans::emmeans(fo, ~ 1, type = "response"))), "\n")
cat(arm, "weights(wt * 2):", r(frm(y | weights(wt * 2) ~ x, data = d, family = poisson())), "\n")
