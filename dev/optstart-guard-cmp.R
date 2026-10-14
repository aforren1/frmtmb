# Lane optstart: compare two runs of dev/optstart-guard.R, case by case:
# how many starting values are identical in hex, and the values that
# moved.
#   Rscript dev/optstart-guard-cmp.R <base.tsv> <lane.tsv>
a <- commandArgs(TRUE)
B <- read.delim(a[1], stringsAsFactors = FALSE)
N <- read.delim(a[2], stringsAsFactors = FALSE)
M <- merge(B, N, by = c("case", "par"), all = TRUE, suffixes = c(".b", ".n"))
M$same <- !is.na(M$hex.b) & !is.na(M$hex.n) & M$hex.b == M$hex.n
for (cs in unique(M$case)) {
  x <- M[M$case == cs, ]
  cat(sprintf("%-28s %2d of %2d identical in hex\n", cs, sum(x$same),
              nrow(x)))
  for (i in which(!x$same)) {
    cat(sprintf("    %-22s base %-14s lane %-14s\n", x$par[i],
                format(x$value.b[i], digits = 7),
                format(x$value.n[i], digits = 7)))
  }
}
