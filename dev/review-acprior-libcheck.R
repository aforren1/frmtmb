# Reviewer: is the lane library the worktree source? Every function of
# R/*.R of frmtmb, frmtmb.sample and frmtmb.spatial, deparsed, against
# the installed namespace. Log: dev/review-acprior-log/libcheck.txt
.libPaths(c("C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
for (p in c("frmtmb", "frmtmb.sample", "frmtmb.spatial")) {
  dir <- if (p == "frmtmb") file.path(wt, "R") else
    file.path(wt, "extensions", p, "R")
  ns <- asNamespace(p)
  e <- new.env(parent = ns)
  for (f in list.files(dir, "[.]R$", full.names = TRUE)) {
    tryCatch(sys.source(f, envir = e, keep.source = FALSE),
             error = function(err) NULL)
  }
  fns <- Filter(function(n) is.function(e[[n]]), ls(e, all.names = TRUE))
  diffs <- Filter(function(n) {
    !exists(n, envir = ns, inherits = FALSE) ||
      !identical(deparse(body(e[[n]])), deparse(body(get(n, envir = ns))))
  }, fns)
  cat(p, ":", length(fns), "source functions;", length(diffs), "differ:",
      paste(head(diffs, 10), collapse = ", "), "\n")
}
