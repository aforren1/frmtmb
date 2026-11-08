# Prints marginaleffects 1.0.0's delta-method internals, to read how it forms
# the Jacobian and the standard error.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("marginaleffects")
a <- commandArgs(TRUE)
for (f in a) {
  if (exists(f, envir = ns)) {
    cat("\n=====", f, "\n")
    print(get(f, envir = ns))
  }
}
