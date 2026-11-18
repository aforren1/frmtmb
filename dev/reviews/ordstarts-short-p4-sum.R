# Reviewer of lane ordstarts, build round: summary of probe 4. The build's
# starts = 10 (the opt-in route the warning names) against 0.83.0's
# starts = 10, no random effects.
#   Rscript ordstarts-short-p4-sum.R > ordstarts-short-p4-sum.txt
O <- "dev/reviews/ordstarts-build-review-p4-"
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE)
b <- rd(paste0(O, "base.tsv"))
l <- rd("dev/reviews/ordstarts-short-p4-lane.tsv")
key <- function(X) paste(X$fam, X$n, X$seed, X$arm)
l <- l[match(key(b), key(l)), ]
stopifnot(identical(key(b), key(l)))
s0b <- b$arm == "s0"
cat("starts = 0: build vs 0.83.0, max |dll|:",
    format(max(abs(l$ll[s0b] - b$ll[s0b]))), "\n\n")
B <- b[b$arm == "s10", ]
L <- l[l$arm == "s10", ]
L0 <- l[l$arm == "s0", ]
silent <- function(X) X$nw == 0 & X$snw == 0
for (des in unique(paste(B$fam, B$n))) {
  i <- paste(B$fam, B$n) == des
  d <- L$ll[i] - B$ll[i]
  cat(sprintf(paste0("%-12s fits %2d | build ran no start %2d (0.83.0: %2d)",
                     " | build below 0.83.0 by > 1e-3 %2d, of them 0.83.0",
                     " silent %2d and the build warned %2d | build above %2d |",
                     " degenerate warning build %2d, 0.83.0 %2d\n"),
              des, sum(i), sum(L$ran[i] == 0), sum(B$ran[i] == 0),
              sum(d < -1e-3), sum(d < -1e-3 & silent(B[i, ])),
              sum(d < -1e-3 & silent(B[i, ]) & L$degen[i]),
              sum(d > 1e-3), sum(L$degen[i]), sum(B$degen[i])))
  j <- which(i & L$ll < B$ll - 1e-3 & L$ran == 0)
  for (k in j) {
    cat(sprintf("   seed %2d: build %.4f (ran 0, warns %s) ; 0.83.0 %.4f (ran %d, silent %s)\n",
                B$seed[k], L$ll[k], L$degen[k], B$ll[k], B$ran[k],
                silent(B[k, ])))
  }
}
