.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
ns <- asNamespace("brms")
show <- function(nm) {
  cat("\n######## ", nm, "\n")
  f <- get(nm, envir = ns)
  print(f)
}
for (nm in c("resp_rate", "resp_subset", "resp_index", "resp_cat",
             "resp_thres", "mi", "data_response", "stan_log_lik_rate",
             "get_ad_values")) {
  if (exists(nm, envir = ns)) show(nm) else cat("\n## missing:", nm, "\n")
}
cat("\n######## grep for rate in namespace function bodies\n")
fns <- ls(ns)
hit <- function(pat) {
  out <- character()
  for (f in fns) {
    obj <- get(f, envir = ns)
    if (!is.function(obj)) next
    txt <- paste(deparse(obj), collapse = "\n")
    if (grepl(pat, txt)) out <- c(out, f)
  }
  out
}
cat("rate:", hit("\\$rate|\"rate\"|denom"), "\n")
cat("subset:", hit("\\$subset|\"subset\"|resp_subset"), "\n")
cat("index:", hit("\\$index|\"index\"|idxl|Jmi"), "\n")
cat("cat:", hit("\\$cat\\b|\"cat\""), "\n")
