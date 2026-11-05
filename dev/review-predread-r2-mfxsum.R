# Short check of fix round 2: the rate of materially off SEs from the review's designs.
O <- "dev/review-predread-out/r2/"
fs <- list.files(O, pattern = "^mfxrate-.*tsv$", full.names = TRUE)
r <- do.call(rbind, lapply(fs, read.delim))
r <- r[is.finite(r$ratio), ]
cat("families", paste(unique(r$family), collapse = " "), " fits", nrow(unique(r[, c("family", "J", "m", "sdg", "slope", "seed")])), "\n")
off <- function(x, lo = 0.8, hi = 1.25) mean(x < lo | x > hi)
tab <- function(by) {
  a <- aggregate(ratio ~ ., r[, c(by, "ratio")], function(x) c(n = length(x), off = round(off(x), 2), low = round(mean(x < 0.8), 2),
                                                               high = round(mean(x > 1.25), 2), med = round(median(x), 2),
                                                               min = round(min(x), 2), max = round(max(x), 2)))
  print(do.call(data.frame, a), row.names = FALSE)
}
cat("\nBy quantity:\n")
tab("qty")
cat("\nBy family and quantity:\n")
tab(c("family", "qty"))
cat("\nBy design, sdg and slope:\n")
tab(c("qty", "sdg", "slope"))
cat("\nLinear link, slope with no random slope, the wording's 'reads those coefficients alone':\n")
s <- r[r$family == "gaussian" & !r$slope & r$qty %in% c("slope_x", "cmp_f", "cmp_w"), ]
print(aggregate(ratio ~ qty, s, function(x) round(c(min = min(x), med = median(x), max = max(x)), 3)))
cat("\nA fit with any material quantity, share by quantity set:\n")
k <- split(r, interaction(r$family, r$J, r$m, r$sdg, r$slope, r$seed, drop = TRUE))
av <- vapply(k, function(d) {
  x <- d[d$qty != "pred", ]
  c(any_avg = any(x$ratio < 0.8 | x$ratio > 1.25),
    any_slope_cmp = any(x$qty %in% c("slope_x", "cmp_f", "cmp_w") & (x$ratio < 0.8 | x$ratio > 1.25)),
    any_pred = any(d$qty == "pred" & (d$ratio < 0.8 | d$ratio > 1.25)))
}, c(NA, NA, NA))
print(round(rowMeans(av), 2))
