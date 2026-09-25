# Lane wt-mvprior: the verdict of each specification in brms against
# each frmtmb arm, read off the three probe logs. Prints one row per
# specification and the counts of agreement.
#
# Usage: Rscript dev/mvprior-compare.R > dev/mvprior-log/compare.txt
root <- "C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-log"
verdicts <- function(f) {
  x <- readLines(file.path(root, f))
  case <- ""
  out <- list()
  for (l in x) {
    if (startsWith(l, "== ")) {
      case <- sub("^== ([^:]+):.*$", "\\1", l)
      next
    }
    m <- regmatches(l, regexec("^(.*\\]) +(ACCEPT|REFUSE)", l))[[1]]
    if (length(m) == 3L) {
      out[[length(out) + 1L]] <- data.frame(case = case,
                                            spec = trimws(m[2]),
                                            v = m[3])
    }
  }
  do.call(rbind, out)
}
options(width = 200)
b <- verdicts("brms-probe.txt")
base <- verdicts("probe-base.txt")
lane <- verdicts("probe-lane.txt")
stopifnot(identical(b$spec, base$spec), identical(b$spec, lane$spec))
tab <- data.frame(case = b$case, spec = b$spec, brms = b$v,
                  base = base$v, lane = lane$v)
tab$base_agrees <- tab$brms == tab$base
tab$lane_agrees <- tab$brms == tab$lane
print(tab, right = FALSE, row.names = FALSE)
cat("\nspecifications:", nrow(tab), "\n")
cat("base agrees with brms:", sum(tab$base_agrees), "of", nrow(tab), "\n")
cat("lane agrees with brms:", sum(tab$lane_agrees), "of", nrow(tab), "\n")
cat("\nlane disagreements:\n")
print(tab[!tab$lane_agrees, c("case", "spec", "brms", "lane")],
      right = FALSE, row.names = FALSE)
