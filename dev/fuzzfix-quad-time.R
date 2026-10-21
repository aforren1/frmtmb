# Stage timings of a nested quadrature fit (verbose = TRUE).
#   Rscript dev/fuzzfix-quad-time.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
d <- readRDS("dev/fuzzfix-out/data-cum.rds")
t0 <- proc.time()[[3]]
f <- frm(bf(y | weights(w) ~ 1 + x + (1 | ga/gb)) + cumulative(), data = d,
         quadrature = TRUE, verbose = TRUE)
cat("fit total", proc.time()[[3]] - t0, "s\n")
p <- f$opt$par
t0 <- proc.time()[[3]]; for (i in 1:5) f$obj$fn(p); cat("fn x5", proc.time()[[3]] - t0, "\n")
t0 <- proc.time()[[3]]; for (i in 1:5) f$obj$gr(p); cat("gr x5", proc.time()[[3]] - t0, "\n")
t0 <- proc.time()[[3]]; s <- summary(f); cat("summary", proc.time()[[3]] - t0, "\n")
print(f$opt[c("iterations", "evaluations")])
