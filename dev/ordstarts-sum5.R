# Lane ordstarts, build round: the K >= 3 rule. On the converging K = 3
# design (dev/ordstarts-swap.R k3 none 3000, sep 2; three families,
# dev/ordstarts-out/swapk3s.tsv): per seed, the default fit A (starts =
# 0), the 10 jittered starts B, and each of the 5 assignments in both
# variants; which assignments reach a good point above A, whether the
# three transpositions (132, 213, 321) alone would have, and the cost.
# Good: code 0, proper, no degenerate and no lost-standard-error warning.
#   Rscript dev/ordstarts-sum5.R > dev/ordstarts-out/sum5.txt
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE, quote = "",
                                    comment.char = "")
X <- rd("dev/ordstarts-out/swapk3s.tsv")
good <- function(D) D$code %in% 0 & D$cls %in% "proper" &
  !(D$degW %in% TRUE) & !(D$lostW %in% TRUE)
tr <- c("perm132", "perm213", "perm321")
A <- X[X$what == "A", ]
cat("seeds", nrow(A), "; A code 0", sum(A$code %in% 0), "; A class",
    paste(names(table(A$cls)), table(A$cls), collapse = ", "), "\n\n")
res <- NULL
for (s in A$seed) {
  a <- A[A$seed == s, ]
  b <- X[X$seed == s & X$what == "B", ]
  p <- X[X$seed == s & grepl("^perm", X$what), ]
  g <- good(p) & p$ll - a$ll > 1e-3
  best_all <- if (any(g)) max(p$ll[g]) else a$ll
  gt <- g & p$what %in% tr
  best_tr <- if (any(gt)) max(p$ll[gt]) else a$ll
  win <- if (any(g)) paste(unique(p$what[g & p$ll >= best_all - 1e-3]),
                           collapse = ",") else ""
  res <- rbind(res, data.frame(seed = s, codeA = a$code, clsA = a$cls,
    up_B = b$ll - a$ll, up_all = best_all - a$ll, up_tr = best_tr - a$ll,
    n_good_up = sum(g), winners = win,
    cost_all = sum(p$elapsed[p$variant == "raw"]) / a$elapsed,
    cost_tr = sum(p$elapsed[p$variant == "raw" & p$what %in% tr]) /
      a$elapsed,
    cost_B = b$elapsed / a$elapsed))
}
print(res, digits = 3, row.names = FALSE)
cat("\nseeds where a good swap end is above A by > 1e-3: all assignments",
    sum(res$up_all > 1e-3), "; transpositions only", sum(res$up_tr > 1e-3),
    "; the 10 jittered starts", sum(res$up_B > 1e-3), "\n")
cat("seeds where the transpositions fall short of all assignments by",
    "> 1e-3:", sum(res$up_all - res$up_tr > 1e-3), "\n")
cat("seeds where all assignments reach the starts' gain or more:",
    sum(res$up_all >= res$up_B - 1e-3), "of", nrow(res), "\n")
cat("cost over the fit (raw variant), median: all 5",
    signif(median(res$cost_all), 3), "; 3 transpositions",
    signif(median(res$cost_tr), 3), "; 10 starts",
    signif(median(res$cost_B), 3), "\n")
