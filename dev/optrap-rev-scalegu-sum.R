# Reviewer of lane optrap: fits with a convergence warning per case and
# scale, lane build as is against the lane with the trip-wire in full
# par_units (dev/optrap-rev-scale.R with OPTRAP_REV_GU=full), and the
# objectives of the two.
cw <- function(w) grepl("maximum absolute gradient|did not report", w)
for (cs in c("gauss_log", "lognormal", "gamma", "student_sigx",
             "gauss_sigx", "skewn", "weights", "offset", "mixture",
             "bigmean", "mv")) {
  a <- read.delim(sprintf("dev/optrap-rev-log/scale-%s-lane.tsv", cs),
                  quote = "")
  b <- read.delim(sprintf("dev/optrap-rev-log/scale-%s-lanegu.tsv", cs),
                  quote = "")
  ta <- tapply(cw(a$warn), a$scale, sum)
  tb <- tapply(cw(b$warn), b$scale, sum)
  cat(sprintf("%-12s lane %s | full units %s | max |obj diff| %.2g\n", cs,
              paste(ta, collapse = " "), paste(tb, collapse = " "),
              max(abs(a$ll - b$ll), na.rm = TRUE)))
}
