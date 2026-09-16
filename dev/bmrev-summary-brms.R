## Reviewer: does brms's OWN summary() report the same R-hat and ESS as
## brms's rhat()/neff_ratio()? If it does, "summary(ds) still reports
## rstan's numbers" is a choice this lane made, not a trade brms makes.
LIB <- "C:/Users/adf44/source/r/bmrev-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
b <- asNamespace("brms")
for (nm in c("summary.brmsfit")) {
  cat("=========== brms:::", nm, " ===========\n", sep = "")
  print(get(nm, envir = b))
}
cat("\n=========== the summary helper it calls ===========\n")
src <- utils::capture.output(print(get("summary.brmsfit", envir = b)))
hits <- grep("posterior_summary|summarise_draws|rhat|ess_|Rhat|Bulk|Tail",
             src, value = TRUE)
cat(paste(hits, collapse = "\n"), "\n")
for (nm in c("brmssummary", "summarise_draws", "posterior_summary.default",
             "get_criterion")) {
  f <- tryCatch(get(nm, envir = b), error = function(e) NULL)
  if (is.null(f)) next
  cat("\n=========== brms:::", nm, " ===========\n", sep = "")
  print(f)
}
cat("\n=========== grep brms sources for Bulk_ESS ===========\n")
for (nm in ls(b)) {
  f <- tryCatch(get(nm, envir = b), error = function(e) NULL)
  if (!is.function(f)) next
  s <- tryCatch(paste(deparse(f), collapse = " "),
                error = function(e) "")
  if (grepl("Bulk_ESS", s, fixed = TRUE)) {
    cat("---- ", nm, " ----\n", sep = "")
    print(f)
  }
}
cat("DONE\n")
