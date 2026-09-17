# Reviewer, lane wt-priorform: is the gp() fit difference between ref
# and lane a lane change or run-to-run noise? Fits the same gp model
# three times in ONE process and prints full-precision values.
#   Rscript dev/priorform-rev-gp.R ref|lane
mode <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (mode == "lane") "C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(20260916)
n <- 180
d <- data.frame(x = runif(n, 0, 3))
d$y <- sin(d$x) + rnorm(n, 0, .3)
for (i in 1:3) {
  f <- suppressWarnings(suppressMessages(frm(bf(y ~ gp(x)) + gaussian(), data = d)))
  cat(mode, i, sprintf("%.17g", f$opt$objective),
      paste(sprintf("%.17g", f$opt$par), collapse = " "), "\n")
}
