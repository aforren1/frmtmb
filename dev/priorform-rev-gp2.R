# Reviewer, lane wt-priorform: does the gp() fit depend on the RNG state
# at the call, and does the state reaching it differ between builds?
# Part A's data, then gp fitted after set.seed(1), set.seed(2), and after
# the 19 models part A fits before it (their RNG use is what differs).
#   Rscript dev/priorform-rev-gp2.R ref|lane
mode <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (mode == "lane") "C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(20260916)
n <- 180
d <- data.frame(x = runif(n, 0, 3))
d$y <- 1 + 0.5 * d$x + rnorm(n)
for (s in 1:2) {
  set.seed(s)
  f <- suppressWarnings(suppressMessages(frm(bf(y ~ gp(x)) + gaussian(), data = d)))
  cat(mode, "seed", s, sprintf("%.17g", f$opt$objective), "\n")
}
