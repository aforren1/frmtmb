# Taping transient of single RTMB operations on long vectors: the peak
# working set during MakeADFun against the tape size.
#   Rscript dev/nestquad-transient-probe.R <op> <n>
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-nestquad-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
op <- a[1]; n <- as.integer(a[2])
set.seed(1)
y <- rnorm(n)
f <- switch(op,
  dnorm = function(p) -sum(dnorm(y, p$m + y * 0.1, exp(p$s), log = TRUE)),
  dt = function(p) -sum(dt((y - p$m) / exp(p$s), exp(p$s) + y * 0, log = TRUE)),
  lgamma = function(p) -sum(lgamma(exp(p$s) + y^2 * 0.1 + p$m)),
  plogis = function(p) -sum(log(plogis(p$m + y))),
  exp = function(p) -sum(exp(p$m + y * p$s)))
gc()
m0 <- ps::ps_memory_info()
obj <- MakeADFun(f, list(m = 0, s = 0), silent = TRUE)
m1 <- ps::ps_memory_info()
info <- utils::getFromNamespace("info", "TMB")(obj$env$ADFun)
cat(sprintf("TRANS op %s n %d values %d peak_delta %.0f MB rss_delta %.0f MB bytes_per_value %.0f\n",
            op, n, info$values_size,
            (m1[["peak_wset"]] - m0[["peak_wset"]]) / 2^20,
            (m1[["rss"]] - m0[["rss"]]) / 2^20,
            (m1[["peak_wset"]] - m0[["rss"]]) / info$values_size))
