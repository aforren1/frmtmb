# Is the installed lane build the worktree source? Compare every function
# defined in the lane-touched R files against the installed namespace.
source("dev/simnewdata-review/rv-prelude.R")
ns <- asNamespace("frmtmb")
for (f in c("R/simulate-newdata.R", "R/predict.R", "R/re-formula.R",
            "R/conditional-effects.R")) {
  e <- new.env()
  sys.source(f, envir = e, keep.source = FALSE)
  nm <- ls(e, all.names = TRUE)
  diff <- nm[!vapply(nm, function(n) {
    if (!exists(n, ns, inherits = FALSE)) return(FALSE)
    a <- get(n, e); b <- get(n, ns)
    if (!is.function(a)) return(identical(a, b))
    identical(deparse(body(a)), deparse(body(b))) &&
      identical(formals(a), formals(b))
  }, NA)]
  cat(f, ": ", length(nm), " objects, differ: ",
      paste(diff, collapse = ", "), "\n", sep = "")
}
