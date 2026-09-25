# Does an installed build match the worktree's R/ sources? Every
# function and constant defined in the sources is compared with the
# installed namespace's, by deparse.
# Usage: Rscript dev/phase3b-libmatch.R <lib> <pkg>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(a[2], character.only = TRUE))
ns <- asNamespace(a[2])
cat("installed from", find.package(a[2]), "\n")
env <- new.env(parent = ns)
for (f in list.files(file.path("extensions", a[2], "R"), full.names = TRUE)) {
  sys.source(f, envir = env, keep.source = FALSE)
}
nm <- ls(env, all.names = TRUE)
diff <- character(0)
for (n in nm) {
  if (!exists(n, envir = ns, inherits = FALSE)) {
    diff <- c(diff, paste("missing in build:", n)); next
  }
  x <- get(n, env); y <- get(n, ns)
  if (is.function(x)) {
    if (!identical(deparse(x), deparse(y))) diff <- c(diff, paste("differs:", n))
  } else if (is.list(x) && all(vapply(x, is.function, TRUE))) {
    if (!identical(lapply(x, deparse), lapply(y, deparse))) {
      diff <- c(diff, paste("differs:", n))
    }
  } else if (!is.environment(x) && !identical(x, y)) {
    diff <- c(diff, paste("differs:", n))
  }
}
cat(length(nm), "source objects compared;", length(diff), "differ\n")
if (length(diff)) cat(diff, sep = "\n")
