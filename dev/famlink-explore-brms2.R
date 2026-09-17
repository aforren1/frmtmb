.libPaths(c("C:/Users/adf44/source/r/pinlib",
             "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
print(brms:::.brmsfamily)
print(brms::mixture)
print(brms:::links_dpars)
print(brms:::family_names)
print(brms:::family.brmsfit)
