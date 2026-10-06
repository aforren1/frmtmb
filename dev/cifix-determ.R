# Is the fit of dev/cifix-scan2.R at n = 140 deterministic? Fit the same
# data set in one process, print the estimates to 17 digits.
# Usage: Rscript dev/cifix-determ.R <lib or "base"> <n> <seed> [prefit]
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
n <- as.integer(args[2])
s <- as.integer(args[3])
if (length(args) > 3) {
  # fit the earlier data sets of the scan first, in this process
  for (nn in c(100L, 120L)) for (ss in seq_len(s)) {
    set.seed(ss)
    d <- data.frame(x = sort(stats::runif(nn, 0, 6)),
                    fac = factor(rep(c("A", "B"), length.out = nn)))
    d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
      stats::rnorm(nn, 0, 0.3)
    invisible(suppressWarnings(
      frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)), data = d)))
  }
}
set.seed(s)
d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
  stats::rnorm(n, 0, 0.3)
fit <- suppressWarnings(
  frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)), data = d))
cat(sprintf("%.17g", fit$opt$par), "\n")
cat("evals:", fit$opt$evals, " objective:", sprintf("%.17g", fit$opt$objective),
    "\n")
