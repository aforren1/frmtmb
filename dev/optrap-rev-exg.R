# Reviewer of lane optrap: exgaussian fits with a convergence warning per
# scale, lane as is against the trip-wire in full par_units.
cw <- function(w) grepl("maximum absolute gradient|did not report", w)
for (b in c("r7", "lane", "lanegu")) {
  x <- read.delim(sprintf("dev/optrap-rev-log/scale-exgauss-%s.tsv", b),
                  quote = "")
  cat(b, "convergence-warned by scale:",
      paste(tapply(cw(x$warn), x$scale, sum), collapse = " "),
      "| code != 0:", paste(tapply(x$code != 0, x$scale, sum, na.rm = TRUE),
                           collapse = " "),
      "| errors:", paste(tapply(is.na(x$code), x$scale, sum), collapse = " "),
      "\n")
}
x <- read.delim("dev/optrap-rev-log/scale-exgauss-r7.tsv", quote = "")
print(unique(substr(x$warn[is.na(x$code)], 1, 150)))
