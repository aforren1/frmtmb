## Reviewer: does the installed frmtmb match the tree's R/ sources?
## Deparse-compares every function the tree defines with the namespace's.
.libPaths(c("C:/Users/adf44/source/r/wt-mointerop-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-mointerop"
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
cat("built:", packageDescription("frmtmb")$Built, "\n")
env <- new.env(parent = ns)
files <- list.files(file.path(root, "R"), "[.][Rr]$", full.names = TRUE)
nerr <- 0
for (f in files) {
  ex <- parse(f, keep.source = FALSE)
  for (e in ex) {
    # only top-level assignments of functions
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
        is.call(e[[3]]) && identical(as.character(e[[3]][[1]]), "function")) {
      r <- tryCatch(eval(e, env), error = function(err) NULL)
      if (is.null(r)) nerr <- nerr + 1
    }
  }
}
nm <- ls(env, all.names = TRUE)
same <- 0; diff <- character(0); missing <- character(0)
for (n in nm) {
  if (!exists(n, envir = ns, inherits = FALSE)) {
    missing <- c(missing, n); next
  }
  if (!is.function(get(n, envir = ns))) { diff <- c(diff, paste0(n, "(not fn)")); next }
  a <- deparse(get(n, envir = env), control = NULL)
  b <- deparse(removeSource(get(n, envir = ns)), control = NULL)
  if (identical(a, b)) same <- same + 1 else diff <- c(diff, n)
}
cat("functions in tree:", length(nm), "identical:", same,
    "different:", length(diff), "not in namespace:", length(missing),
    "eval errors:", nerr, "\n")
if (length(diff)) cat("DIFFERENT:", diff, "\n")
if (length(missing)) cat("MISSING:", head(missing, 20), "\n")
