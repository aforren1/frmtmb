# Does the lane library hold the tree's R code? Source each package's
# R/ into an environment and compare every function with the installed
# namespace's (deparsed, srcref dropped).
.libPaths(c("C:/Users/adf44/source/r/wt-sigmainn-lib", "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
for (p in c("frmtmb", "frmtmb.sample")) {
  dir <- if (p == "frmtmb") file.path(wt, "R") else file.path(wt, "extensions", p, "R")
  ns <- asNamespace(p)
  env <- new.env(parent = ns)
  for (f in sort(list.files(dir, "[.]R$", full.names = TRUE)))
    tryCatch(sys.source(f, envir = env, keep.source = FALSE), error = function(e) cat("source error", f, conditionMessage(e), "\n"))
  nms <- ls(env, all.names = TRUE)
  fn <- nms[vapply(nms, function(n) is.function(get(n, env)), NA)]
  dd <- function(f) paste(deparse(f, control = NULL), collapse = "\n")
  diff <- fn[vapply(fn, function(n) {
    if (!exists(n, ns, inherits = FALSE)) return(TRUE)
    !identical(dd(get(n, env)), dd(get(n, ns)))
  }, NA)]
  cat(p, "from", dirname(find.package(p)), ": functions", length(fn), "differing", length(diff), "\n")
  if (length(diff)) print(head(diff, 20))
}
