.libPaths(c("C:/Users/adf44/source/r/rellib-r9", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
# Reviewer of lane sampleprior: what brms 2.23.0 does with uniform() on
# class b without lb/ub
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
r <- withCallingHandlers(brms::stancode(y ~ x, data = D$wide,
       prior = brms::set_prior("uniform(2.2, 2.6)", class = "b", coef = "x")),
       warning = function(w) { cat("WARNING:", conditionMessage(w), "\n");
         invokeRestart("muffleWarning") })
s <- strsplit(as.character(r), "\n")[[1]]
cat(grep("uniform|vector\[Kc\] b", s, value = TRUE), sep = "\n")
