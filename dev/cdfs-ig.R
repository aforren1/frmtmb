# inverse.gaussian's lcdf where exp(2 shape / mu) overflows: the base
# build against this lane's, one process each (argument: library), and
# statmod as the reference
a <- commandArgs(TRUE)
.libPaths(c(if (length(a)) a[1], "C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
f <- frmtmb:::as_frmtmb_family(inverse.gaussian())
q <- c(0.5, 0.9, 1, 1.1, 2)
F <- f$lcdf(q, list(mu = rep(1, 5), shape = rep(1000, 5)), list())
cat(dirname(find.package("frmtmb")), "\n")
cat("lcdf at mu = 1, shape = 1000:", format(F, digits = 10), "\n")
cat("statmod:                     ", format(statmod::pinvgauss(q, 1, 1000), digits = 10), "\n")
