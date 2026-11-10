# Compare the worktree's definitions of the m1-m4 functions with the
# copies in the lane library's installed namespace. Only the named
# top-level assignments are evaluated, so the rest of each file need not
# run outside the package.
#   Rscript dev/compat-m-nscheck.R
.libPaths(c("C:/Users/adf44/source/r/wt-compat-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(getNamespaceInfo(ns, "path")),
    as.character(utils::packageVersion("frmtmb")), "\n")
fns <- c("compat_ordinal_only_rules", "registered_aterm_of",
         "aterm_spelling")
src <- new.env()
for (f in c("R/compat.R", "R/parse.R")) {
  for (e in parse(f, keep.source = FALSE)) {
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
        is.name(e[[2]]) && as.character(e[[2]]) %in% fns)
      assign(as.character(e[[2]]), eval(e[[3]]), envir = src)
  }
}
for (fn in fns) {
  a <- get0(fn, envir = src, inherits = FALSE)
  b <- get0(fn, envir = ns, inherits = FALSE)
  same <- !is.null(a) && !is.null(b) &&
    identical(deparse(a), deparse(b))
  cat(sprintf("%-28s source:%s installed:%s identical:%s\n", fn,
              !is.null(a), !is.null(b), same))
}
