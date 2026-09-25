# Merge the four 400-bit reference chunks written by
# dev/phase3b-cdf-reference-v2.R (run as: Rscript
# dev/phase3b-cdf-reference-v2.R 400 <1..4> 4) into one table, and
# say how far the two independent 400-bit series agree.
# Output: dev/phase3b-log/cdf-reference-400.csv and .txt, and the test
# fixture extensions/frmtmb.eam/tests/testthat/fixtures/wiener-cdf-ref.csv
x <- do.call(rbind, lapply(1:4, function(k) {
  utils::read.csv(sprintf("dev/phase3b-log/cdf-ref2-400-chunk%d.csv", k),
                  colClasses = c(F_ref = "character", S_ref = "character"))
}))
x <- x[order(x$grid != "density", x$v, x$w, x$a, x$t), ]
stopifnot(nrow(x) == 1125 + 192, !anyDuplicated(x[, c("t", "a", "w", "v",
                                                      "grid")]))
utils::write.csv(x, "dev/phase3b-log/cdf-reference-400.csv", row.names = FALSE)
fx <- x[, c("t", "a", "w", "v", "grid", "F_ref", "S_ref")]
utils::write.csv(fx, "extensions/frmtmb.eam/tests/testthat/fixtures/wiener-cdf-ref.csv",
                 row.names = FALSE)
out <- c(sprintf("400-bit reference: %d rows (%d density grid, %d tail)",
                 nrow(x), sum(x$grid == "density"), sum(x$grid == "tail")),
         sprintf("rows where the small-time image series was also summed: %d",
                 sum(!is.na(x$agree))),
         sprintf("worst relative disagreement of the two 400-bit series: %.3g",
                 max(x$agree, na.rm = TRUE)),
         sprintf("smallest F: %.3g; smallest S: %.3g",
                 min(as.numeric(x$F_ref)), min(as.numeric(x$S_ref))))
writeLines(out, "dev/phase3b-log/cdf-reference-400.txt")
cat(out, sep = "\n")
