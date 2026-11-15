# Lane wcheck2: wcheck's case table on rellib-r21 against the lane, from
# the logs of dev/wcheck2-accept.sh (dev/wcheck2-out/accept/), and the
# line-by-line comparison of the guard-absent models and the route
# table. Prints the rows that differ and the count of rows that agree.
#   Rscript dev/wcheck2-table.R
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
  b <- arm_of(file.path(dir, sprintf("cases-%s-base.txt", p)))
  l <- arm_of(file.path(dir, sprintf("cases-%s-lane.txt", p)))
  cat(sprintf("\nPart %s: %d arms\n", p, length(b)))
  for (a in names(b)) {
    tot <- tot + 1
    lv <- if (is.null(l[[a]])) "absent" else l[[a]]
    if (identical(b[[a]], lv)) {
      same <- same + 1
    } else {
      cat(sprintf("| %s | %s | %s |\n", gsub("|", "\\|", a, fixed = TRUE),
                  b[[a]], lv))
    }
  }
}
cat(sprintf("\ncase table: %d of %d rows identical on rellib-r21 and the lane\n",
            same, tot))
for (nm in c("common", "routes")) {
  b <- readLines(file.path(dir, sprintf("%s-base.txt", nm)), warn = FALSE)
  l <- readLines(file.path(dir, sprintf("%s-lane.txt", nm)), warn = FALSE)
  b <- b[!startsWith(b, "lib:")]
  l <- l[!startsWith(l, "lib:")]
  d <- which(b != l)
  cat(sprintf("%s: %d lines, %d differ (apart from lib:)\n", nm, length(b),
              length(d)))
  for (i in d) cat("  base:", b[i], "\n  lane:", l[i], "\n")
}
