# Reviewer: does the lane library hold the tree's code? Sources every
# R file of the tree into an environment and compares each function
# with the installed namespace's, sources removed.
#   Rscript dev/review-priordraws-libmatch.R [lib]
lib <- commandArgs(TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/wt-priordraws-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
W <- "C:/Users/adf44/source/r/frmtmb-wt-priordraws"
for (p in list(c("frmtmb", W), c("frmtmb.sample", file.path(W, "extensions/frmtmb.sample")))) {
  suppressPackageStartupMessages(loadNamespace(p[1]))
  ns <- asNamespace(p[1])
  env <- new.env(parent = ns)
  for (f in list.files(file.path(p[2], "R"), "[.][Rr]$", full.names = TRUE)) {
    tryCatch(sys.source(f, env, keep.source = FALSE),
             error = function(e) cat("  source error", basename(f), conditionMessage(e), "\n"))
  }
  nms <- ls(env, all.names = TRUE)
  fn <- nms[vapply(nms, function(n) is.function(get(n, env)), NA)]
  dep <- function(x) {
    x <- utils::removeSource(x)
    environment(x) <- globalenv()
    deparse(x)
  }
  diff <- character(0); miss <- character(0)
  for (n in fn) {
    if (!exists(n, envir = ns, inherits = FALSE)) { miss <- c(miss, n); next }
    if (!identical(dep(get(n, env)), dep(get(n, ns)))) diff <- c(diff, n)
  }
  cat(sprintf("%s (%s): %d functions in the tree, %d differ, %d missing from the namespace\n",
              p[1], find.package(p[1]), length(fn), length(diff), length(miss)))
  if (length(diff)) cat("  differ:", head(diff, 20), "\n")
  if (length(miss)) cat("  missing:", head(miss, 20), "\n")
}
