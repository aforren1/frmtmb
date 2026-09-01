## Third probe: confirm the RTMBp tape is actually parallelized, and diff
## the R sources of the two packages.
out <- file("dev/rtmbp-probe3.txt", open = "wt")
p <- function(...) cat(..., "\n", file = out, sep = "")

## --- R-source diff -----------------------------------------------------
dep <- function(pkg, nm) {
  ns <- asNamespace(pkg)
  ob <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
  if (is.null(ob)) return("<absent>")
  paste(deparse(ob), collapse = "\n")
}
ns1 <- ls(asNamespace("RTMB"), all.names = TRUE)
ns2 <- ls(asNamespace("RTMBp"), all.names = TRUE)
p("objects only in RTMBp: ", paste(setdiff(ns2, ns1), collapse = ", "))
p("objects only in RTMB : ", paste(setdiff(ns1, ns2), collapse = ", "))
diffs <- character()
for (nm in intersect(ns1, ns2)) {
  a <- dep("RTMB", nm); b <- dep("RTMBp", nm)
  if (!identical(a, b)) diffs <- c(diffs, nm)
}
p("\nR objects whose source DIFFERS between RTMB and RTMBp:")
p(paste(diffs, collapse = ", "))
for (nm in diffs) {
  p("\n### ", nm, " -- RTMB ###"); p(dep("RTMB", nm))
  p("### ", nm, " -- RTMBp ###"); p(dep("RTMBp", nm))
}
close(out)
cat("wrote dev/rtmbp-probe3.txt\n")
