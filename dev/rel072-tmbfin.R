# Consolidation item 2: TMB's ADFun finalizer prints "promise already
# under evaluation" once per session when a garbage collection forces
# TMB:::isNullPointer's lazy-load promise with other finalizers pending
# (dev/illse2-findings.md, section 4; TMB-4 in dev/upstream-bugs.md).
# The fit is lane illse2's: a degree-7 raw polynomial beside K = 150
# standard normal covariates, gaussian, n = 2000, so 159 outer
# parameters (dev/illse2-p4wide.R 150 2000 1:2). Run it in a fresh
# session per library and count the lines on stderr:
#
#   Rscript dev/rel072-tmbfin.R <lib> 2> err.txt
#
# <lib> is put ahead of the user library: wt-release-lib (frmtmb 0.72.0
# with the .onLoad() guard) or wt-illse2-lib (the lane, without it).
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
K <- 150
n <- 2000
for (seed in 1:2) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(sample(1:50, n, TRUE)))
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  Z <- matrix(stats::rnorm(n * K), n, K)
  colnames(Z) <- paste0("z", seq_len(K))
  d <- cbind(d, Z)
  d$y <- sin(3 * d$x) + 0.05 * rowSums(Z[, 1:10]) + stats::rnorm(n, 0, 0.5)
  rhs <- paste(c(paste0("x", 1:7), colnames(Z)), collapse = " + ")
  fo <- stats::as.formula(paste("y ~", rhs))
  f <- suppressWarnings(frm(fo, data = d))
  cat(sprintf("seed %d | outer parameters %d | logLik %.6f\n", seed,
              length(f$opt$par), as.numeric(logLik(f))))
  flush(stdout())
}
invisible(gc())
cat("DONE\n")
