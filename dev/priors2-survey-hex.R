# Lane priors2: the 33-fit survey of two builds compared in hex. Each
# value is written with %.17g, which reads back as the same double, so
# the hex of the value read is the hex of the value fitted.
# Rscript dev/priors2-survey-hex.R <before.tsv> <after.tsv>
a <- commandArgs(trailingOnly = TRUE)
rd <- function(f) {
  s <- strsplit(readLines(f)[-1L], "\t", fixed = TRUE)
  stats::setNames(lapply(s, function(x) {
    if (x[2L] != "OK") return(NA_character_)
    v <- as.numeric(c(x[4L], strsplit(x[5L], ",", fixed = TRUE)[[1L]]))
    paste(sprintf("%a", v), collapse = " ")
  }), vapply(s, `[`, "", 1L))
}
b <- rd(a[1])
f <- rd(a[2])
same <- vapply(names(b), function(k) identical(b[[k]], f[[k]]), NA)
for (k in names(b)) {
  cat(sprintf("%-18s %s\n", k, if (same[[k]]) "bitwise identical" else
    if (is.na(b[[k]])) "before: error" else "differs"))
}
cat(sum(same), "of", length(same), "bitwise identical\n")
