# What the six fixture fits in tests.brmsfit-methods.R actually are.
# tests.brmsfit-methods.R, tests.posterior_epred.R, tests.emmeans.R,
# tests.data-helpers.R and tests.priorsense.R all run on these, so the
# port cost of those files is the cost of rebuilding these designs.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n\n")
for (i in 1:6) {
  f <- get(paste0("brmsfit_example", i), envir = asNamespace("brms"))
  cat("== fit", i, "==\n")
  print(f$formula)
  cat("ndraws:", brms::ndraws(f), " nobs:", nobs(f), "\n")
  cat("vars:", paste(head(brms::variables(f), 40), collapse = " "), "\n\n")
}
