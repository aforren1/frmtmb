# Fix round 1: the cost of a truncated poisson under importance
# sampling (test-importance.R, "trunc() and cens() are corrected"),
# which hit std::bad_alloc once in a one-session run: time and the R
# process's peak memory for the fit, on the given library.
#   Rscript dev/cdfs-fix1-imp.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", find.package("frmtmb"), "\n")
set.seed(21)
ng <- 40L
per <- 6L
g <- factor(rep(seq_len(ng), each = per))
x <- rnorm(ng * per)
u <- rnorm(ng, 0, 0.8)
y <- rpois(ng * per, exp(1.1 + 0.4 * x + u[g]))
dd <- data.frame(y = y, x = x, g = g)[y >= 1L, ]
bfm <- bf(y | trunc(lb = 1) ~ x + (1 | g)) + poisson()
t0 <- proc.time()[["elapsed"]]
fi <- frm(bfm, dd, importance = 500L)
cat(sprintf("importance 500: objective %.6f seconds %.1f\n", fi$opt$objective,
            proc.time()[["elapsed"]] - t0))
t0 <- proc.time()[["elapsed"]]
fl <- frm(bfm, dd)
cat(sprintf("laplace: objective %.6f seconds %.1f\n", fl$opt$objective,
            proc.time()[["elapsed"]] - t0))
mem <- system2("powershell", c("-NoProfile", "-Command",
  sprintf("(Get-Process -Id %d).PeakWorkingSet64 / 1MB", Sys.getpid())),
  stdout = TRUE)
cat("peak memory (MB, memory.size):", format(mem), "\n")
cat("END\n")
