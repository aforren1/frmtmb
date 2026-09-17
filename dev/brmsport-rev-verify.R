# Reviewer (brmsport): does the lane build match the worktree source?
# Compares every top-level function assignment in R/ against the
# installed namespace, by deparse without srcref.
lib <- "C:/Users/adf44/source/r/brmsport-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-brmsport"
chk <- function(pkg, dir) {
  ns <- asNamespace(pkg)
  cat(pkg, as.character(packageVersion(pkg)), find.package(pkg), "\n")
  n <- 0; bad <- character()
  for (f in list.files(file.path(dir, "R"), full.names = TRUE)) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
          is.name(e[[2]]) && is.call(e[[3]]) &&
          identical(e[[3]][[1]], as.name("function"))) {
        nm <- as.character(e[[2]])
        if (!exists(nm, envir = ns, inherits = FALSE)) next
        obj <- get(nm, envir = ns)
        if (!is.function(obj)) next
        src <- eval(e[[3]])
        n <- n + 1
        a <- deparse(removeSource(src)); b <- deparse(removeSource(obj))
        if (!identical(a, b)) bad <- c(bad, nm)
      }
    }
  }
  cat("compared", n, "functions; differ:", length(bad), "\n")
  if (length(bad)) print(head(bad, 20))
}
chk("frmtmb", root)
chk("frmtmb.sample", file.path(root, "extensions", "frmtmb.sample"))
cat("StanHeaders", as.character(packageVersion("StanHeaders")),
    "brms", as.character(packageVersion("brms")), "\n")
