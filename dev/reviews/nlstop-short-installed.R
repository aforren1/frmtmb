# Short check: the installed prior_conflict_escape() and fit_assembled()
# equal the tree's.
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
e <- new.env(parent = ns)
for (f in c("R/fit-restart.R", "R/fit.R")) sys.source(f, envir = e)
for (fn in c("prior_conflict_escape", "fit_assembled", "nl_limit_message"))
  cat(fn, identical(deparse(get(fn, e)), deparse(get(fn, ns))), "\n")
