# Reviewer: does the worker's library hold the worktree's code? Every
# function defined in R/ of core and of frmtmb.sample is sourced into a
# fresh environment and compared, body and formals, with the installed
# namespace object. Log: dev/aterms2-rev-log-00-integrity.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-aterms2"
cmp <- function(pkg, rdir) {
  ns <- asNamespace(pkg)
  env <- new.env(parent = ns)
  for (f in list.files(rdir, "[.]R$", full.names = TRUE)) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
          is.name(e[[2]]) && is.call(e[[3]]) &&
          identical(e[[3]][[1]], as.name("function"))) {
        nm <- as.character(e[[2]])
        assign(nm, eval(e[[3]], env), envir = env)
      }
    }
  }
  nms <- ls(env, all.names = TRUE)
  bad <- character(0)
  miss <- character(0)
  for (nm in nms) {
    if (!exists(nm, envir = ns, inherits = FALSE)) {
      miss <- c(miss, nm)
      next
    }
    a <- get(nm, envir = env)
    b <- get(nm, envir = ns)
    if (!is.function(b)) next
    if (!identical(deparse(body(a)), deparse(body(b))) ||
        !identical(deparse(formals(a)), deparse(formals(b)))) {
      bad <- c(bad, nm)
    }
  }
  cat(pkg, ": functions compared", length(nms), " differing", length(bad),
      " absent from namespace", length(miss), "\n")
  if (length(bad)) cat("  differ:", bad, "\n")
  if (length(miss)) cat("  absent:", head(miss, 20), "\n")
}
cmp("frmtmb", file.path(wt, "R"))
cmp("frmtmb.sample", file.path(wt, "extensions/frmtmb.sample/R"))
cat("frmtmb path", find.package("frmtmb"), "\n")
cat("sample path", find.package("frmtmb.sample"), "\n")
