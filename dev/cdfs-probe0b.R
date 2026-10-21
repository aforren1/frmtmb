.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
print(RTMBdist::dskewnorm2); print(RTMBdist::dexgauss); print(RTMBdist::pinvgauss); print(RTMBdist::dbetabinom); print(RTMBdist::dvm)
print(sapply(c("Rmpfr","sn","extraDistr","circular","CircStats","statmod","gsl","pracma","numDeriv"), requireNamespace, quietly=TRUE))
