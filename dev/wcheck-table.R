# Lane wcheck: the before/after verdict table, generated from the logs
# of dev/wcheck-cases.R (dev/wcheck-out/cases-<part>-<base|lane>.txt).
#   Rscript dev/wcheck-table.R
arm_of <- function(f) {
  x <- readLines(f, warn = FALSE)
  x <- x[startsWith(x, "ARM ")]
  nm <- trimws(substr(x, 5, 50))
  rest <- trimws(substring(x, 51))
  out <- vapply(rest, function(r) {
    if (startsWith(r, "REFUSED")) return("refused")
    if (startsWith(r, "ERROR")) {
      return(paste0("error: ", substr(sub("^ERROR \\| ", "", r), 1, 40)))
    }
    ll <- sub("^ok logLik ([^ ]+) .*", "\\1", r)
    bd <- sub(".*betad ([^|]*)\\|.*", "\\1", r)
    cv <- sub(".*conv ([^ ]+) .*", "\\1", r)
    sprintf("fit %s, betad %s, conv %s", ll, trimws(bd), cv)
  }, "")
  stats::setNames(unname(out), nm)
}
for (p in c("mv", "cens", "prior")) {
  b <- arm_of(sprintf("dev/wcheck-out/cases-%s-base.txt", p))
  l <- arm_of(sprintf("dev/wcheck-out/cases-%s-lane.txt", p))
  cat(sprintf("\nPart %s (%d arms)\n\n| arm | rellib-r20 | lane |\n|---|---|---|\n",
              p, length(b)))
  for (a in names(b)) {
    cat(sprintf("| %s | %s | %s |\n", gsub("|", "\\|", a, fixed = TRUE),
                b[[a]], if (is.null(l[[a]])) "absent" else l[[a]]))
  }
}
