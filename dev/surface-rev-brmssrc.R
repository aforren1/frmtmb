# Reviewer: print brms 2.23.0 sources the review reads.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("brms")
for (f in commandArgs(TRUE)) { cat("=====", f, "\n"); print(get(f, ns)) }
