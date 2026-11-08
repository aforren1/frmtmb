# Summarize dev/mfxse-rate.R's output, one arm: the off-rate by quantity, by
# family and quantity, the ratio of each level's SE to fitted()'s Est.Error,
# the warnings marginaleffects calls raised, and the seconds per fit.
#   Rscript dev/mfxse-ratesum.R <arm>    (before | after)
arm <- commandArgs(TRUE)[1]
O <- "C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-rate-out/"
fs <- list.files(O, pattern = paste0("^", arm, "-.*tsv$"), full.names = TRUE)
r0 <- do.call(rbind, lapply(fs, read.delim))
cat("arm", arm, " files", length(fs), "\n")
key <- c("family", "J", "m", "sdg", "slope", "seed")
cat("families", paste(sort(unique(r0$family)), collapse = " "), " fits",
    nrow(unique(r0[, key])), " rows", nrow(r0), " non-finite ratio",
    sum(!is.finite(r0$ratio)), "\n")
r <- r0[is.finite(r0$ratio), ]
off <- function(x, lo = 0.8, hi = 1.25) mean(x < lo | x > hi)
tab <- function(by) {
  a <- aggregate(ratio ~ ., r[, c(by, "ratio")], function(x) {
    c(n = length(x), off = round(off(x), 2), low = round(mean(x < 0.8), 2),
      high = round(mean(x > 1.25), 2), med = round(median(x), 2),
      min = round(min(x), 2), max = round(max(x), 2))
  })
  print(do.call(data.frame, a), row.names = FALSE)
}
cat("\nBy quantity:\n")
tab("qty")
cat("\nBy family and quantity:\n")
tab(c("family", "qty"))
cat("\nA fit with any material quantity, share:\n")
k <- split(r, interaction(r[, key], drop = TRUE))
av <- vapply(k, function(d) {
  x <- d[d$qty != "pred", ]
  c(any_avg = any(x$ratio < 0.8 | x$ratio > 1.25),
    any_pred = any(d$qty == "pred" & (d$ratio < 0.8 | d$ratio > 1.25)))
}, c(NA, NA))
print(round(rowMeans(av), 3))
cat("\nLevels against fitted()'s Est.Error (me_se / fitted_se):\n")
p <- r0[r0$qty == "pred" & is.finite(r0$fitted_se) & is.finite(r0$me_se), ]
q <- p$me_se / p$fitted_se
cat(sprintf("n %d  min %.6f  median %.6f  max %.6f  outside 1e-3: %d\n",
            length(q), min(q), median(q), max(q), sum(abs(q - 1) > 1e-3)))
cat("\nWarnings raised per fit (class of the first of each kind):\n")
w <- unique(r0[, c(key, "warn")])
print(table(ifelse(is.na(w$warn) | w$warn == "", "(none)", w$warn)))
cat("\nSeconds for the five marginaleffects calls, per fit, by J:\n")
s <- unique(r0[, c(key, "secs")])
print(aggregate(secs ~ J, s, function(x) round(c(med = median(x),
                                                  max = max(x)), 2)))
cat("\nConvergence codes:\n")
print(table(unique(r0[, c(key, "conv")])$conv))
