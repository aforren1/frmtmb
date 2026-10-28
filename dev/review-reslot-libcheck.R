# Reviewer: does the lane library hold the worktree source? Deparse every
# function of R/ and compare with the installed namespace. Run from the
# worktree root.
.libPaths(c("C:/Users/adf44/source/r/wt-reslot-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (p in c("frmtmb", "frmtmb.sample")) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  root <- if (p == "frmtmb") "." else file.path("extensions", p)
  ns <- asNamespace(p)
  env <- new.env(parent = ns)
  for (f in list.files(file.path(root, "R"), "[.]R$", full.names = TRUE)) {
    try(sys.source(f, envir = env, keep.source = FALSE), silent = TRUE)
  }
  nms <- Filter(function(n) is.function(env[[n]]), ls(env, all.names = TRUE))
  diff <- 0; missing <- 0
  for (n in nms) {
    if (!exists(n, envir = ns, inherits = FALSE)) {
      missing <- missing + 1
      next
    }
    a <- deparse(env[[n]]); b <- deparse(get(n, envir = ns))
    if (!identical(a, b)) { diff <- diff + 1; cat("  differs:", n, "\n") }
  }
  cat(p, ": functions", length(nms), " differing", diff, " missing", missing,
      "\n")
}
