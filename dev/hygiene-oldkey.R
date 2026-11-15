# F1: the new factory block of test-register-once.R against the lane's
# first key (namespace scan, else code alone), put back into the lane
# build with assignInNamespace(). It must fail there.
#   Rscript dev/hygiene-oldkey.R
.libPaths(c("C:/Users/adf44/source/r/wt-hygiene-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
old <- function(fn) {
  env <- environment(fn) %||% baseenv()
  ns <- topenv(env)
  if (isNamespace(ns)) {
    for (nm in ls(ns, all.names = TRUE)) {
      if (!bindingIsActive(nm, ns) &&
          identical(get(nm, envir = ns, inherits = FALSE), fn)) {
        return(paste0(environmentName(ns), "::", nm))
      }
    }
  }
  paste0(environmentName(ns), "|", paste(deparse(fn), collapse = "\n"))
}
environment(old) <- asNamespace("frmtmb")
assignInNamespace("registration_key", old, "frmtmb")
r <- as.data.frame(test_file("tests/testthat/test-register-once.R",
                             package = "frmtmb",
                             env = test_env("frmtmb"), reporter = "silent"))
print(r[, c("test", "passed", "failed", "error")], row.names = FALSE)
