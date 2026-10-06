# Reviewer: lme4's allFit() as installed: version, the start_from_mle
# default, and how each setting refits.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
cat("lme4", as.character(packageVersion("lme4")), "\n")
f <- lme4::allFit
cat("start_from_mle default:", deparse(formals(f)$start_from_mle), "\n")
b <- deparse(f)
i <- grep("start_from_mle|update", b)
cat(b[sort(unique(c(i - 1, i, i + 1)))], sep = "\n")
