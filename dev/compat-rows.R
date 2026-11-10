# Print frm_compat() rows for the ten pairs of the compat lane, on the
# library named by the first argument (default: the base build).
#   Rscript dev/compat-rows.R [lib;lib...]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[1] else
  "C:/Users/adf44/source/r/rellib-r19"
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
pairs <- list(
  p167 = c("poisson", "arma()"), p158 = c("fcor()", "nl"),
  p171 = c("zero_inflated_beta", "nl"), p169 = c("von_mises", "mi()"),
  p166 = c("nbinom1", "nl"), p014 = c("cs_pred()", "nl"),
  p146 = c("trials()", "cs_pred()"), p010 = c("cens()", "disc"),
  p011 = c("cens()", "equidistant"), p059 = c("vint()", "sum_to_zero"))
cp <- as.data.frame(frm_compat())
for (id in names(pairs)) {
  a <- pairs[[id]][1]; b <- pairs[[id]][2]
  hit <- cp[(cp$feature_a == a & cp$feature_b == b) |
              (cp$feature_a == b & cp$feature_b == a), , drop = FALSE]
  cat("\n==", id, a, "x", b, "\n")
  if (nrow(hit) == 0L) cat("  (no row)\n")
  for (i in seq_len(nrow(hit))) {
    cat("  status:", hit$status[i], "\n")
    for (nm in setdiff(names(hit), c("feature_a", "feature_b", "status"))) {
      cat("  ", nm, ": ", paste(hit[[nm]][i], collapse = " "), "\n", sep = "")
    }
  }
}
