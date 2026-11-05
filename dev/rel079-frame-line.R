# The corrected-weights line of the survey under each frame: 0.77.0's, which the
# sample was drawn from and whose strata.tsv the headline's shares are, and this
# build's. A frame moves the line only through its strata shares.
#   Rscript dev/rel079-frame-line.R
rd <- function(f) utils::read.delim(f, colClasses = "character", na.strings = "NA")
res <- rd("dev/equiv-data/results.tsv")
res$draws <- as.numeric(res$draws)
res$agree <- as.logical(res$agree)
wn <- rd("dev/equiv-data/weights.tsv")
wd <- rd("dev/equiv-data/weights-drawn.tsv")
wn <- stats::setNames(as.numeric(wn$weight), wn$feature)
wd <- stats::setNames(as.numeric(wd$weight), wd$feature)
strata <- c("untested", "conditional", "works", "refused")
ratio <- (wn[res$feature_a] * wn[res$feature_b]) / (wd[res$feature_a] * wd[res$feature_b])
rr <- vapply(strata, function(s) {
  i <- res$status == s
  sum((res$draws * ratio * res$agree)[i]) / sum((res$draws * ratio)[i])
}, 0)
for (nm in c("dev/equiv-frame-out/compat.tsv", "dev/rel079-log/frame-r17/compat.tsv")) {
  cf <- rd(nm)
  cf <- cf[cf$in_core_only == "TRUE", ]
  pw <- wn[cf$feature_a] * wn[cf$feature_b]
  W <- vapply(strata, function(s) sum(pw[cf$status == s], na.rm = TRUE), 0)
  W <- W / sum(W)
  cat(nm, ": shares", paste(sprintf("%s %.3f", strata, W), collapse = ", "),
      "; line", sprintf("%.1f", 100 * sum(W * rr)), "\n")
}
