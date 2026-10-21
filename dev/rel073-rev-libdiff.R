# Release reviewer, 0.73.0: are wt-release-lib (the build every scan and
# session ran on) and rellib-r11 (built from the final tree) the same
# code? Every object of each package's namespace, deparsed (whitespace
# and comments do not survive), compared by package.
#   Rscript dev/rel073-rev-libdiff.R <pkg> <libA> <libB>
a <- commandArgs(trailingOnly = TRUE)
dump_ns <- function(lib, pkg) {
  out <- tempfile(fileext = ".rds")
  code <- sprintf(paste0(
    ".libPaths(c('%s', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'));",
    "ns <- asNamespace('%s'); nm <- sort(ls(ns, all.names = TRUE));",
    "v <- vapply(nm, function(n) { x <- get(n, envir = ns);",
    "if (is.function(x)) x <- removeSource(x);",
    "paste(deparse(x, control = c('keepNA', 'keepInteger')), collapse = '\n') },",
    "''); saveRDS(list(v = v, path = find.package('%s')), '%s')"),
    lib, pkg, pkg, gsub("\\\\", "/", out))
  system2(file.path(R.home("bin"), "Rscript.exe"), c("-e", shQuote(code)))
  readRDS(out)
}
A <- dump_ns(a[2], a[1]); B <- dump_ns(a[3], a[1])
cat(a[1], ":", A$path, "vs", B$path, "\n")
nm <- union(names(A$v), names(B$v))
only <- setdiff(nm, intersect(names(A$v), names(B$v)))
diffs <- Filter(function(n) !identical(A$v[[n]], B$v[[n]]),
                intersect(names(A$v), names(B$v)))
cat("objects", length(nm), "| only in one:", length(only),
    paste(only, collapse = " "), "| differ:", length(diffs),
    paste(diffs, collapse = " "), "\n")
