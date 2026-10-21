# Cost of the iterated Gauss-Kronrod tape by nesting depth: tape time,
# fit time and the process's peak working set, one depth per process.
#   Rscript dev/fuzzfix-quad-depth.R <lib|base> <depth 1|2|3> <ng>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
depth <- as.integer(a[2]); ng <- as.integer(a[3])
set.seed(7)
d <- expand.grid(r = 1:2, gc = c("p", "q"), gb = c("i", "ii"),
                 ga = factor(seq_len(ng)))
d$x <- rnorm(nrow(d))
d$y <- rpois(nrow(d), exp(0.5 + 0.3 * d$x + rnorm(ng, 0, 0.6)[d$ga]))
re <- c("(1 | ga)", "(1 | ga/gb)", "(1 | ga/gb/gc)")[depth]
fo <- stats::as.formula(paste("y ~ x +", re))
t0 <- proc.time()[[3]]
f <- tryCatch(suppressMessages(suppressWarnings(
  frm(bf(fo) + poisson(), data = d, quadrature = TRUE))),
  error = function(e) e)
el <- proc.time()[[3]] - t0
peak <- system2("powershell", c("-NoProfile", "-Command",
  shQuote(sprintf("(Get-Process -Id %d).PeakWorkingSet64 / 1MB",
                  Sys.getpid()))), stdout = TRUE)
cat(sprintf("DEPTH %d rows %d nb %s time %.1fs peakMB %s %s\n", depth,
            nrow(d), if (inherits(f, "error")) NA else length(f$estimates$b),
            el, round(as.numeric(peak)),
            if (inherits(f, "error")) paste("ERROR", conditionMessage(f))
            else sprintf("logLik %.6f", as.numeric(logLik(f)))))
