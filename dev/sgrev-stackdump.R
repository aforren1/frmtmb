.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
suppressMessages(library(frmtmb.sample))
obj <- structure(list(), class = "frmtmb_draws")
for (nm in c("log_posterior", "neff_ratio", "nuts_params", "rhat")) {
  g <- get(nm, envir = globalenv())
  cat("\n=====", nm, " generic from ",
      environmentName(environment(g)), "\n", sep = "")
  tryCatch(withCallingHandlers(g(obj), error = function(e) {
    cat("  msg: ", conditionMessage(e), "\n")
    n <- sys.nframe()
    for (i in seq_len(n)) {
      f <- sys.function(i)
      cat(sprintf("   [%2d] %-28s env=%s\n", i,
                  paste(deparse(sys.call(i)[[1]]), collapse = "")[1],
                  environmentName(environment(f))))
    }
  }), error = function(e) invisible(NULL))
}
