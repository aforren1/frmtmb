# Review of lane slopetrig: summarize dev/reviews/slopetrig-review-oos-*.tsv.
# Short = more than 1e-3 below the best of every arm of both libraries.
# Usage: Rscript slopetrig-review-oos-sum.R <base.tsv> <lane.tsv>
a <- commandArgs(TRUE)
B <- utils::read.delim(a[1], quote = "", comment.char = "")
L <- utils::read.delim(a[2], quote = "", comment.char = "")
stopifnot(identical(paste(B$design, B$seed), paste(L$design, L$seed)))
best <- pmax(B$ll_ref, B$ll_F, B$ll_T, B$ll_D, L$ll_ref, L$ll_F, L$ll_T,
             L$ll_D, na.rm = TRUE)
k <- 1e-3
sh <- function(ll, i = seq_along(best)) ifelse(is.na(ll), NA, best[i] - ll > k)
out <- do.call(rbind, lapply(split(seq_len(nrow(B)), factor(B$design, unique(B$design))),
  function(i) {
    data.frame(design = B$design[i[1]], n = length(i),
      engB = sum(B$as_D[i], na.rm = TRUE), engL = sum(L$as_D[i], na.rm = TRUE),
      shF = sum(sh(B$ll_F[i], i), na.rm = TRUE),
      shT = sum(sh(L$ll_T[i], i), na.rm = TRUE),
      shDb = sum(sh(B$ll_D[i], i), na.rm = TRUE),
      silDb = sum(sh(B$ll_D[i], i) & B$w_D[i] == 0, na.rm = TRUE),
      shDl = sum(sh(L$ll_D[i], i), na.rm = TRUE),
      silDl = sum(sh(L$ll_D[i], i) & L$w_D[i] == 0, na.rm = TRUE),
      maxDl = signif(max(best[i] - L$ll_D[i], na.rm = TRUE), 3),
      maxDb = signif(max(best[i] - B$ll_D[i], na.rm = TRUE), 3),
      shRef = sum(sh(L$ll_ref[i], i), na.rm = TRUE),
      wOK_Dl = sum(!sh(L$ll_D[i], i) & L$w_D[i] > 0, na.rm = TRUE),
      wOK_Db = sum(!sh(B$ll_D[i], i) & B$w_D[i] > 0, na.rm = TRUE),
      naD = sum(is.na(L$ll_D[i])),
      tDl_tF = signif(stats::median(L$t_D[i] / pmax(L$t_F[i], 0.01)), 3))
  }))
options(width = 220)
print(out, row.names = FALSE)
cat("\nTOTAL n", nrow(B), "; short D base", sum(sh(B$ll_D), na.rm = TRUE),
    "(silent", sum(sh(B$ll_D) & B$w_D == 0, na.rm = TRUE), ") ; short D lane",
    sum(sh(L$ll_D), na.rm = TRUE), "(silent",
    sum(sh(L$ll_D) & L$w_D == 0, na.rm = TRUE), ")\n")
cat("F arm identical base/lane:", sum(B$ll_F == L$ll_F, na.rm = TRUE), "of",
    nrow(B), "\n")
cat("\nlane D short or newly warned fits:\n")
i <- which(sh(L$ll_D) | (L$w_D > 0 & B$w_D == 0))
print(data.frame(design = L$design[i], seed = L$seed[i],
                 short = signif(best[i] - L$ll_D[i], 3), wD = L$w_D[i],
                 wDb = B$w_D[i], wt = substr(L$wt_D[i], 1, 120)), row.names = FALSE)
cat("\nerrors:\n")
e <- grepl("^ERROR", L$wt_D) | grepl("^ERROR", L$wt_F)
print(unique(data.frame(design = L$design[e], D = substr(L$wt_D[e], 1, 120),
                        F = substr(L$wt_F[e], 1, 60))), row.names = FALSE)
