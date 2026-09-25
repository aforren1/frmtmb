# Does the lane library match the worktree source? Compare every
# function defined in R/*.R of frmtmb and frmtmb.sample with the
# installed namespace, by deparse.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-prelude.R")
suppressPackageStartupMessages(library(frmtmb.sample))
check_pkg <- function(pkg, dir) {
  ns <- asNamespace(pkg)
  env <- new.env(parent = ns)
  n_same <- 0; diff <- character(0); missing <- character(0)
  for (f in list.files(dir, "[.][Rr]$", full.names = TRUE)) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
          is.name(e[[2]]) && is.call(e[[3]]) &&
          identical(e[[3]][[1]], as.name("function"))) {
        nm <- as.character(e[[2]])
        src <- eval(e[[3]], env)
        if (!exists(nm, envir = ns, inherits = FALSE)) {
          missing <- c(missing, nm); next
        }
        inst <- get(nm, envir = ns)
        if (!is.function(inst)) next
        a <- deparse(removeSource(src)); b <- deparse(removeSource(inst))
        if (identical(a, b)) n_same <- n_same + 1 else diff <- c(diff, nm)
      }
    }
  }
  cat(pkg, ": identical", n_same, " differ", length(diff), " missing",
      length(missing), "\n")
  if (length(diff)) cat("  differ:", head(diff, 30), "\n")
}
check_pkg("frmtmb", file.path(R2_ROOT, "R"))
check_pkg("frmtmb.sample", file.path(R2_ROOT, "extensions/frmtmb.sample/R"))
cat("lp_prior_dpar exists:", exists("lp_prior_dpar", asNamespace("frmtmb")), "\n")
cat("brms", format(packageVersion("brms")), "\n")
