# Reviewer of lane defects: does the installed namespace match the
# worktree's R/ sources, function by function? Control: the same
# comparison against rellib-r3 must report exactly the changed ones.
#   Rscript dev/defects-rev-installed.R <lib> <pkg> <srcdir>
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]; pkg <- a[2]; src <- a[3]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(loadNamespace(pkg))
ns <- asNamespace(pkg)
cat("loaded", pkg, format(packageVersion(pkg)), "from",
    dirname(getNamespaceInfo(ns, "path")), "\n")
env <- new.env(parent = ns)
files <- list.files(file.path(src, "R"), pattern = "[.][Rr]$",
                    full.names = TRUE)
for (f in files) {
  tryCatch(sys.source(f, envir = env, keep.source = FALSE),
           error = function(e) cat("SOURCE ERROR", f, conditionMessage(e),
                                   "\n"))
}
nms <- ls(env, all.names = TRUE)
fns <- nms[vapply(nms, function(n) is.function(get(n, env)), NA)]
norm <- function(f) {
  c(deparse(args(f)), deparse(body(f)))
}
diffs <- character(0); missing <- character(0)
for (n in fns) {
  if (!exists(n, envir = ns, inherits = FALSE)) {
    missing <- c(missing, n); next
  }
  g <- get(n, envir = ns, inherits = FALSE)
  if (!is.function(g)) { diffs <- c(diffs, n); next }
  if (!identical(norm(get(n, env)), norm(g))) diffs <- c(diffs, n)
}
# non-function objects too (tables such as residuals_unsupported)
objs <- setdiff(nms, fns)
odiff <- character(0)
for (n in objs) {
  if (!exists(n, envir = ns, inherits = FALSE)) {
    missing <- c(missing, n); next
  }
  if (!identical(get(n, env), get(n, envir = ns, inherits = FALSE))) {
    odiff <- c(odiff, n)
  }
}
cat("functions compared:", length(fns), " objects:", length(objs), "\n")
cat("DIFFERENT functions:", length(diffs), "\n ", paste(diffs, collapse = " "),
    "\n")
cat("DIFFERENT objects:", length(odiff), "\n ", paste(odiff, collapse = " "),
    "\n")
cat("MISSING from namespace:", length(missing), "\n ",
    paste(missing, collapse = " "), "\n")
