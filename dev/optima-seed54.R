# Lane optima: mo_search()'s record on seed 54 of brms_monotonic's data
# code, on the R that runs this script (reference BLAS or the OpenBLAS
# build of dev/optima-openblas.sh).
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
src <- readLines("dev/optima-mo-study.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^# exact maximum", src) - 1L)]))
f <- suppressWarnings(frm(ls ~ mo(income) * age, data = mk(54)))
print(f$opt$mo_search, digits = 17)
cat("logLik", format(as.numeric(logLik(f)), digits = 15), "\n")
