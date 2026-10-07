# Reviewer of lane setier, re-check of item 8: across the singular and
# small-sd studies, fits at optimizer code != 0 that got no convergence
# warning, and what the standard-error analysis lost on them.
f <- function(path, warncol) {
  x <- read.delim(path, stringsAsFactors = FALSE)
  x$lost[is.na(x$lost)] <- ""
  q <- x[x$code != 0, ]
  silent <- q[q[[warncol]] == 0, ]
  allb <- vapply(strsplit(silent$lost, "[;,]"), function(v) {
    length(v) > 0 && all(grepl(":boundary$", v))
  }, NA)
  cat(sprintf("%-40s code!=0 %3d | no convergence warning %3d | of those all-boundary %3d, other %d\n",
              basename(path), nrow(q), nrow(silent), sum(allb), sum(!allb)))
}
for (p in c("dev/setier-rev2-log/sing-lane-ref.tsv",
            "dev/setier-rev2-log/sing-lane-ob.tsv")) f(p, "other_warn")
f("dev/setier-rev2-log/smallsd-lane-ref.tsv", "other_warn")
a <- read.delim("dev/setier-rev2-log/sing-lane-ref.tsv")
b <- read.delim("dev/setier-rev2-log/sing-lane-ob.tsv")
cat("singular study, ref vs OpenBLAS rows differing in lost or message:",
    sum(a$lost != b$lost | a$boundary_msg != b$boundary_msg), "\n")
x <- a[a$lme4_singular & a$boundary_msg == 0, ]
print(x[, c("design", "seed", "code", "lost", "other_warn")])
