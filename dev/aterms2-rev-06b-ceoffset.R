# Reviewer: conditional_effects() on an offset(log(time)) poisson model,
# per arm, to see whether its failure predates the lane.
# Log: dev/aterms2-rev-log-06b-ceoffset.txt
arm <- commandArgs(TRUE)[[1]]
libs <- c("C:/Users/adf44/source/r/wt-aterms2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
set.seed(606)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4))
d$y <- rpois(n, exp(0.2 + 0.5 * d$x) * d$time)
fo <- suppressWarnings(frm(y ~ x + offset(log(time)), data = d,
                           family = poisson()))
r <- tryCatch({conditional_effects(fo); "ok"},
              error = function(e) conditionMessage(e))
cat(arm, find.package("frmtmb"), ":", r, "\n")
