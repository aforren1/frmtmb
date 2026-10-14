# Lane optstart: the battery's two kinked families (asym_laplace, huber),
# before against after: how far short of the best the fits stop, in
# bands, over all five designs (dev/optstart-battery-sum.R's best,
# degenerate stops excluded).
#   Rscript dev/optstart-battery-kink.R <before label> <after label>
a <- commandArgs(TRUE)
rd <- function(lab) {
  fs <- Sys.glob(sprintf("dev/optstart-log/battery/%s-*.tsv", lab))
  do.call(rbind, lapply(fs, read.delim, stringsAsFactors = FALSE, quote = ""))
}
M <- merge(rd(a[1]), rd(a[2]), by = c("family", "design", "seed"),
           suffixes = c(".b", ".n"))
tr <- pmax(M$ll_truth.b, M$ll_truth.n, na.rm = TRUE)
deg <- function(l) !is.na(l) & !is.na(tr) & l > tr + 50
best <- pmax(M$ll_ref.b, M$ll_truth.b, M$ll_truth.n,
             ifelse(deg(M$ll.b), NA, M$ll.b), ifelse(deg(M$ll.n), NA, M$ll.n),
             na.rm = TRUE)
for (fam in c("asym_laplace", "huber")) {
  k <- M$family == fam
  for (arm in c("b", "n")) {
    l <- M[[paste0("ll.", arm)]][k]
    dg <- deg(M[[paste0("ll.", arm)]])[k]
    g <- (l - best[k])[!dg & !is.na(l)]
    cat(sprintf("%-13s %s: n %d | degenerate %d | gap > -1e-4 %d | (-1e-4, -0.1] %d | (-0.1, -1] %d | (-1, -10] %d | < -10 %d | sum of gaps %.1f\n",
                fam, if (arm == "b") "before" else "after ", sum(k),
                sum(dg), sum(g > -1e-4), sum(g <= -1e-4 & g > -0.1),
                sum(g <= -0.1 & g > -1), sum(g <= -1 & g > -10),
                sum(g <= -10), sum(g)))
  }
}
