# Reviewer, lane sampfix: does the worker's library hold the worktree's
# code? Every top-level function in the worktree's R/ files is compared
# with the installed namespace's (formals and body, srcref dropped).
LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
.libPaths(c(LANE, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix"
chk <- function(pkg, dir) {
  ns <- asNamespace(pkg)
  cat(pkg, "from", dirname(find.package(pkg)), "\n")
  env <- new.env()
  for (f in list.files(dir, "\\.R$", full.names = TRUE)) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (is.call(e) && (identical(e[[1]], as.name("<-")) ||
                         identical(e[[1]], as.name("="))) &&
          is.name(e[[2]]) && is.call(e[[3]]) &&
          identical(e[[3]][[1]], as.name("function"))) {
        assign(as.character(e[[2]]), eval(e[[3]], baseenv()), env)
      }
    }
  }
  nm <- ls(env, all.names = TRUE)
  bad <- character()
  for (n in nm) {
    if (!exists(n, ns, inherits = FALSE)) { bad <- c(bad, paste("MISSING", n)); next }
    a <- get(n, ns); b <- get(n, env)
    if (!is.function(a)) next
    a <- utils::removeSource(a); b <- utils::removeSource(b)
    if (!identical(deparse(formals(a)), deparse(formals(b))) ||
        !identical(deparse(body(a)), deparse(body(b)))) bad <- c(bad, n)
  }
  cat("  functions compared:", length(nm), " differing:", length(bad), "\n")
  if (length(bad)) print(bad)
  cat("  installed version:", as.character(packageVersion(pkg)), "\n")
}
chk("frmtmb", file.path(WT, "R"))
chk("frmtmb.sample", file.path(WT, "extensions/frmtmb.sample/R"))
