# Instrument check: the lane library's namespace against the worktree's R/.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
ns <- asNamespace("frmtmb")
files <- list.files("C:/Users/adf44/source/r/frmtmb-wt-bigx/R", "\\.R$",
                    full.names = TRUE)
n <- 0; diff <- character(0); miss <- character(0)
for (f in files) {
  ex <- parse(f, keep.source = FALSE)
  for (e in ex) {
    if (is.call(e) && (identical(e[[1]], as.name("<-")) ||
                       identical(e[[1]], as.name("="))) &&
        is.name(e[[2]]) && is.call(e[[3]]) &&
        identical(e[[3]][[1]], as.name("function"))) {
      nm <- as.character(e[[2]])
      fn <- eval(e[[3]])
      n <- n + 1
      if (!exists(nm, envir = ns, inherits = FALSE)) {
        miss <- c(miss, nm); next
      }
      g <- get(nm, envir = ns)
      if (!identical(deparse(body(fn)), deparse(body(g))) ||
          !identical(deparse(formals(fn)), deparse(formals(g)))) {
        diff <- c(diff, nm)
      }
    }
  }
}
cat("functions", n, "differing", length(diff), "missing", length(miss), "\n")
print(head(diff)); print(head(miss))
cat("design_alias_tol", get("design_alias_tol", ns), "\n")
