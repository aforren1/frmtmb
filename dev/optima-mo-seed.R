# Lane optima, item 1: one seed of dev/optima-mo-study.R in detail: the
# optimizer's message, the gradient and the simplexes at the end, with
# verbose = 2's stages.
#   Rscript dev/optima-mo-seed.R base|lane seed
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
seed <- as.integer(args[2])
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
src <- readLines("dev/optima-mo-study.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^simplex_of", src) - 1L)]))
d <- mk(seed)
ex <- exact_max(d)
f <- frm(ls ~ mo(income) * age, data = d, verbose = TRUE)
ms <- utils::getFromNamespace("mo_simplex", "frmtmb")
cat("logLik", format(as.numeric(logLik(f)), digits = 12), "exact",
    format(ex$ll, digits = 12), "\n")
print(f$opt[c("convergence", "message", "iterations", "evals")])
print(f$opt$par)
print(f$obj$gr(f$opt$par))
z <- f$estimates[grepl("^zeta", names(f$estimates))]
print(lapply(z, ms))
print(ex[c("b_mo", "w1", "w2")])
