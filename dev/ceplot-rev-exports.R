.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib","C:/Users/adf44/source/r/rellib-r4","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
e <- getNamespaceExports("frmtmb"); print(intersect(c("gaussian","cumulative","zero_inflated_poisson","lognormal","bf","poisson"), e))
e2 <- getNamespaceExports("brms"); print(intersect(c("gaussian","cumulative","zero_inflated_poisson","lognormal","bf"), e2))