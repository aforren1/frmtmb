# Probe the cumulative quadrature summary() error on a given library.
#   Rscript dev/fuzzfix-cum-probe.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
d <- readRDS("dev/fuzzfix-out/data-cum.rds")
fit <- suppressWarnings(frm(bf(y | weights(w) ~ 1 + x + (1 | ga/gb)) +
                              cumulative(), data = d, quadrature = TRUE))
o <- fit$obj
cat("random:", deparse(o$env$random), "\n")
cat(".random:", deparse(o$env$.random), "\n")
cat("integrate:", names(o$env$integrate), "\n")
cat("par names:", unique(names(o$env$par)), "\n")
cat("opt par:", names(fit$opt$par), "\n")
sdr <- frmtmb:::sdr_of(fit)
cat("par.random:", length(sdr$par.random), "\n")
str(sdr$env$random)
print(names(sdr))
