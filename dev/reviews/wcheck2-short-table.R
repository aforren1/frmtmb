# Short check of lane wcheck2's fix round: the 58-row case table rerun on
# the lane library (logs in <dir>) against the lane's recorded rellib-r21
# logs (dev/wcheck2-out/accept/cases-<part>-base.txt), by verdict and value
# as dev/wcheck2-table.R reads them; warning text is ignored.
#   Rscript dev/reviews/wcheck2-short-table.R <dir>
dir <- commandArgs(TRUE)[1]
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
tot <- 0
same <- 0
for (p in c("mv", "cens", "prior")) {
  b <- arm_of(sprintf("dev/wcheck2-out/accept/cases-%s-base.txt", p))
  l <- arm_of(file.path(dir, sprintf("cases-%s-lane.txt", p)))
  stopifnot(identical(names(b), names(l)))
  for (a in names(b)) {
    tot <- tot + 1
    if (identical(b[[a]], l[[a]])) same <- same + 1 else
      cat(sprintf("DIFF %s %s\n  base %s\n  lane %s\n", p, a, b[[a]], l[[a]]))
  }
  cat(sprintf("part %s: %d rows\n", p, length(b)))
}
cat(sprintf("case table: %d of %d rows identical (rellib-r21 record vs lane rerun)\n",
            same, tot))
