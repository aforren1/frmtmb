.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
print(extSoftVersion()["BLAS"]); print(La_library()); print(La_version())
set.seed(1); m <- matrix(rnorm(4e6), 2000); print(system.time(m %*% m))
