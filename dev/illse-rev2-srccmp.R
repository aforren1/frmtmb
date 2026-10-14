# Reviewer of lane illse, round 2: every function of the changed R files, as the
# tree has it, against the lane library's installed namespace.
.libPaths(c("C:/Users/adf44/source/r/wt-illse-lib", "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.spline)})
cmp <- function(files, nsname) {
  ns <- asNamespace(nsname)
  nd <- 0; nm <- 0; tot <- 0
  for (f in files) {
    e <- new.env(parent = ns)
    sys.source(f, envir = e, keep.source = FALSE)
    for (n in ls(e, all.names = TRUE)) {
      o <- get(n, e)
      if (!is.function(o)) next
      tot <- tot + 1
      if (!exists(n, ns, inherits = FALSE)) { nm <- nm + 1; cat("missing", n, "\n"); next }
      a <- deparse(o); b <- deparse(get(n, ns))
      if (!identical(a, b)) { nd <- nd + 1; cat("differs", n, "\n") }
    }
  }
  cat(nsname, ": functions", tot, "differ", nd, "missing", nm, "\n")
}
W <- "C:/Users/adf44/source/r/frmtmb-wt-illse"
cmp(file.path(W, "R", c("se-reparam.R", "se-check.R", "fit.R", "predict.R", "confint.R",
                        "interop.R", "methods-fit.R", "multiple.R", "predict-brms.R",
                        "brms-shapes.R", "conditional-smooths.R", "sandwich.R", "insight.R")), "frmtmb")
cmp(file.path(W, "extensions/frmtmb.spline/R", c("curve-cov.R", "curve-deriv.R", "curve-feature.R")),
    "frmtmb.spline")
