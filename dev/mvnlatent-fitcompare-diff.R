# Compare two saves of dev/mvnlatent-fitcompare.R, fit by fit and
# quantity by quantity: identical(), else the largest absolute and
# relative difference.
#
#   Rscript dev/mvnlatent-fitcompare-diff.R <a.rds> <b.rds>
a <- commandArgs(trailingOnly = TRUE)
A <- readRDS(a[1])
B <- readRDS(a[2])
stopifnot(identical(names(A), names(B)))
num <- function(x) {
  if (is.character(x)) return(NULL)
  if (is.list(x)) x <- unlist(x)
  suppressWarnings(as.numeric(unclass(x)))
}
q <- c("logLik", "df", "est", "opt", "se", "fixef", "predict", "fitted",
       "probs", "simulate", "conv")
tab <- NULL
for (id in names(A)) {
  for (k in q) {
    x <- A[[id]][[k]]
    y <- B[[id]][[k]]
    same <- identical(x, y)
    xa <- num(x)
    ya <- num(y)
    d <- if (!same && length(xa) && length(xa) == length(ya)) {
      ok <- is.finite(xa) & is.finite(ya)
      c(max(abs(xa - ya)[ok]), max((abs(xa - ya) / pmax(abs(xa), 1e-300))[ok]))
    } else c(NA, NA)
    tab <- rbind(tab, data.frame(fit = id, what = k, identical = same,
                                 n = length(xa), max_abs = d[1],
                                 max_rel = d[2]))
  }
}
cat("fits", length(A), "; quantities compared", nrow(tab),
    "; identical", sum(tab$identical), "\n")
cat("non-identical by quantity:\n")
print(table(tab$what[!tab$identical]))
if (any(!tab$identical)) {
  cat("largest differences:\n")
  print(utils::head(tab[!tab$identical, ][order(-tab$max_rel[!tab$identical]),
                                          ], 30), row.names = FALSE)
}
# a message (say, a refused vcov) must be the same message
chr <- vapply(names(A), function(id) {
  all(vapply(q, function(k) {
    !is.character(A[[id]][[k]]) || identical(A[[id]][[k]], B[[id]][[k]])
  }, TRUE))
}, TRUE)
cat("character results equal in every fit:", all(chr), "\n")
cat("logLik range:", format(range(vapply(A, `[[`, 0, "logLik")),
                            digits = 10), "\n")
