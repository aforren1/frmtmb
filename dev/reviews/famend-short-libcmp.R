# Short check of lane famend: does the installed core match the tree?
#   Rscript dev/reviews/famend-short-libcmp.R <lib>
# Sources each changed R/ file into a scratch environment and compares
# every function's body and formals with the installed namespace's copy.
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famend"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
files <- c("R/families-extreme.R", "R/compat.R", "R/families-cdf.R",
           "R/families.R", "R/fit.R", "R/predict.R", "R/priors.R")
strip <- function(f) {
  if (is.function(f)) {
    f <- utils::removeSource(f)
    paste(deparse(f, control = c("keepNA", "keepInteger", "niceNames",
                                 "showAttributes")), collapse = "\n")
  } else {
    paste(deparse(f), collapse = "\n")
  }
}
for (fl in files) {
  e <- new.env(parent = ns)
  sys.source(file.path(wt, fl), envir = e, keep.source = FALSE)
  nm <- ls(e, all.names = TRUE)
  diff <- character(0)
  miss <- character(0)
  for (n in nm) {
    if (!exists(n, envir = ns, inherits = FALSE)) {
      miss <- c(miss, n)
      next
    }
    if (!identical(strip(get(n, e)), strip(get(n, ns)))) diff <- c(diff, n)
  }
  cat(sprintf("%-22s objects %4d differ %d missing %d %s %s\n", fl,
              length(nm), length(diff), length(miss),
              paste(diff, collapse = ","), paste(miss, collapse = ",")))
}
