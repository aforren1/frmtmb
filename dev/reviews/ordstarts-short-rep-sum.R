# Reviewer of lane ordstarts, build round: summary of the reproduction
# (ordstarts-build-review-rep.R), questions 5 and 8.
# silent: no warning at the fit and none at the first summary().
# clean B: code 0, no degenerate and no lost-standard-error warning at the
#   fit or summary(), largest threshold Est.Error at most 100 (fixef()).
# shortfall: the better point above by more than 1e-3.
#   Rscript ordstarts-short-rep-sum.R > ordstarts-short-rep-sum.txt
O <- "dev/reviews/ordstarts-build-review-rep-"; OL <- "dev/reviews/ordstarts-short-rep-"
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE)
for (des in c("cs1000", "cc300")) {
  b <- rd(paste0(O, "base-", des, ".tsv"))
  l <- rd(paste0(OL, "lane-", des, ".tsv"))
  A <- b[b$arm == "A", ]
  B <- b[b$arm == "B", ]
  A0 <- l[l$arm == "A0", ]
  D <- l[l$arm == "D", ]
  stopifnot(identical(A$seed, B$seed), identical(A$seed, A0$seed),
            identical(A$seed, D$seed))
  silent <- function(X) X$nw == 0 & X$snw == 0
  clean <- function(X) X$code %in% 0 & !X$degen & !X$lost & X$maxtse <= 100
  before <- B$ll > A$ll + 1e-3 & silent(A) & clean(B)
  after <- B$ll > D$ll + 1e-3 & silent(D) & clean(B)
  cat("==", des, "== seeds", nrow(A), "\n")
  cat("A0 (lane starts = 0) vs 0.83.0 default A, max |dll|:",
      format(max(abs(A0$ll - A$ll))), "\n")
  cat("silent shortfalls before (B over A):", sum(before), "seeds",
      paste(A$seed[before], collapse = ","), "\n")
  cat("silent shortfalls after (B over D):", sum(after), "seeds",
      paste(A$seed[after], collapse = ","), "\n")
  cat("D up over A0:", sum(D$ll > A0$ll + 1e-3), "; down:",
      sum(D$ll < A0$ll - 1e-3), "; D above B:", sum(D$ll > B$ll + 1e-3),
      "; D below B:", sum(D$ll < B$ll - 1e-3), "; swap kept:",
      sum(D$swapped == 1), "\n")
  cat("degenerate warning (fit or summary): A0", sum(A0$degen), "D",
      sum(D$degen), "B", sum(B$degen), "\n")
  r <- D$el / A0$el
  k <- D$swapped == 1
  ran <- D$swaps > 0
  cat("cost D/A0: median", format(median(r), digits = 3), "; where the swap",
      "was kept", if (any(k)) format(median(r[k]), digits = 3) else "-",
      "(n", sum(k), ") ; ran and not kept",
      if (any(ran & !k)) format(median(r[ran & !k]), digits = 3) else "-",
      "(n", sum(ran & !k), ") ; B/A median",
      format(median(B$el / A$el), digits = 3), "\n")
  print(data.frame(seed = A$seed, A = round(A$ll, 4), B = round(B$ll, 4),
                   A0 = round(A0$ll, 4), D = round(D$ll, 4),
                   silentA = silent(A), cleanB = clean(B),
                   silentD = silent(D), swapped = D$swapped,
                   DoverA0 = round(r, 2)), row.names = FALSE)
}
