# Reviewer (check of punch 2c): does the installed lane build match the tree's R/ sources?
.libPaths(c("C:/Users/adf44/source/r/wt-seflags-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-seflags/R"
nd <- 0
for (f in c("se-check.R", "predict.R", "confint.R", "brms-shapes.R", "conditional-smooths.R", "extra-deriv.R", "interop.R", "methods-fit.R", "predict-brms.R")) {
  e <- new.env(parent = ns)
  sys.source(file.path(wt, f), envir = e, keep.source = FALSE)
  for (nm in ls(e, all.names = TRUE)) {
    a <- get(nm, e)
    if (!exists(nm, ns, inherits = FALSE)) { cat("missing", nm, "\n"); next }
    b <- get(nm, ns)
    if (is.function(a)) {
      same <- identical(deparse(a), deparse(b))
    } else same <- identical(a, b)
    if (!same) { nd <- nd + 1; cat("DIFF", f, nm, "\n") }
  }
}
cat("differences:", nd, "\n")
