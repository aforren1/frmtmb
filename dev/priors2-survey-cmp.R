# Compare two runs of dev/priors2-survey.R: per case, the codes, the
# objectives, and the largest move of an outer parameter, relative to
# max(1, |before|).
# Rscript dev/priors2-survey-cmp.R <before.tsv> <after.tsv>
a <- commandArgs(trailingOnly = TRUE)
rd <- function(f) {
  x <- readLines(f)[-1L]
  s <- strsplit(x, "\t", fixed = TRUE)
  data.frame(case = vapply(s, `[`, "", 1L), status = vapply(s, `[`, "", 2L),
             code = vapply(s, `[`, "", 3L), obj = vapply(s, `[`, "", 4L),
             par = vapply(s, `[`, "", 5L), stringsAsFactors = FALSE)
}
b <- rd(a[1])
f <- rd(a[2])
stopifnot(identical(b$case, f$case))
num <- function(s) as.numeric(strsplit(s, ",", fixed = TRUE)[[1L]])
cat(sprintf("%-18s %-6s %-6s %5s %5s %22s %22s %10s\n", "case", "before",
            "after", "c0", "c1", "obj before", "obj after", "max move"))
for (i in seq_len(nrow(b))) {
  mv <- NA_real_
  if (b$status[i] == "OK" && f$status[i] == "OK") {
    p0 <- num(b$par[i])
    p1 <- num(f$par[i])
    if (length(p0) == length(p1)) {
      mv <- max(abs(p1 - p0) / pmax(1, abs(p0)))
    }
  }
  cat(sprintf("%-18s %-6s %-6s %5s %5s %22s %22s %10.3g\n", b$case[i],
              b$status[i], f$status[i], b$code[i], f$code[i],
              substr(b$obj[i], 1, 22), substr(f$obj[i], 1, 22), mv))
}
