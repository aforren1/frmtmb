# Counts from dev/nanse-mo-escape.R's TSV for the findings.
#   Rscript dev/nanse-escape-sum.R [tsv]
f <- commandArgs(TRUE)[1]
if (is.na(f)) f <- "dev/nanse-log/mo-escape.tsv"
X <- utils::read.delim(f)
cat("seeds", nrow(X), "\n")
for (v in c("gap_fit", "gap_escape")) {
  cat(v, ": > 1e-3", sum(X[[v]] > 1e-3), " > 1e-2", sum(X[[v]] > 1e-2),
      " > 0.1", sum(X[[v]] > 0.1), " max", format(max(X[[v]]), digits = 3),
      "\n")
}
r <- X$rounds > 0
cat("seeds where a vertex direction improved:", sum(r), "; log-likelihood",
    "gained there:", format(X$gap_fit[r] - X$gap_escape[r], digits = 3), "\n")
