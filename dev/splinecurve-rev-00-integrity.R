# Reviewer, lane splinecurve: is each library what it is said to be?
#   Rscript splinecurve-rev-00-integrity.R <lib>
# Prints an md5 over every function of the frmtmb and frmtmb.spline
# namespaces (deparsed, no srcref), and, for frmtmb.spline, compares each
# installed function with the tree's current source.
a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
options(keep.source = FALSE)
ns_hash <- function(pkg) {
  ns <- asNamespace(pkg)
  nm <- sort(ls(ns, all.names = TRUE))
  txt <- unlist(lapply(nm, function(n) {
    f <- get(n, envir = ns)
    if (is.function(f)) c(n, deparse(f)) else NULL
  }))
  tf <- tempfile()
  writeLines(txt, tf)
  unname(tools::md5sum(tf))
}
for (p in c("frmtmb", "frmtmb.spline")) {
  cat(p, format(packageVersion(p, lib.loc = LIB)), "from",
      dirname(find.package(p)), "md5", ns_hash(p), "\n")
}
cat("frm_curve formals:", paste(names(formals(frmtmb.spline::frm_curve)),
                                collapse = ", "), "\n")

# installed frmtmb.spline against the tree's source
src <- "C:/Users/adf44/source/r/frmtmb-wt-release/extensions/frmtmb.spline/R"
e <- new.env(parent = asNamespace("frmtmb.spline"))
for (f in list.files(src, full.names = TRUE, pattern = "[.]R$")) {
  sys.source(f, envir = e, keep.source = FALSE)
}
ns <- asNamespace("frmtmb.spline")
diffs <- character(0)
n_cmp <- 0L
for (n in ls(e, all.names = TRUE)) {
  s <- get(n, envir = e)
  if (!is.function(s)) next
  if (!exists(n, envir = ns, inherits = FALSE)) {
    diffs <- c(diffs, paste0(n, " (not installed)"))
    next
  }
  i <- get(n, envir = ns)
  n_cmp <- n_cmp + 1L
  if (!identical(deparse(s), deparse(i))) diffs <- c(diffs, n)
}
cat("functions compared with the tree's source:", n_cmp,
    " differing:", length(diffs), "\n")
if (length(diffs)) cat("  ", paste(diffs, collapse = ", "), "\n")
