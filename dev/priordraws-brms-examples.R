# The priors and prior draw columns of brms 2.23.0's stored example
# fits, which the ported brms suite's fixtures stand in for.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
for (k in 1:6) {
  ex <- get(paste0("brmsfit_example", k), asNamespace("brms"))
  cat("\n===== example", k, "=====\n")
  cat("sample_prior attr:", format(attr(ex$prior, "sample_prior")), "\n")
  print(ex$prior[nzchar(ex$prior$prior) & ex$prior$source == "user", ])
  v <- variables(ex)
  cat("prior vars:", paste(grep("^prior_", v, value = TRUE), collapse = " "),
      "\n")
  cat("sim: iter", ex$fit@sim$iter, "warmup", ex$fit@sim$warmup,
      "warmup2", ex$fit@sim$warmup2, "chains", ex$fit@sim$chains, "\n")
}
