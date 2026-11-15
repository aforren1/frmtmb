.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
src <- deparse(TMB::MakeADFun)
h <- grep("integrate", src)
cat(src[min(h):(max(h) + 30)], sep = "\n")
cat("---- GK\n"); print(TMB:::GK)
