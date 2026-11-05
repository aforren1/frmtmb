for (f in list.files("dev", "^review-nlfix-latres-.*tsv$", full.names = TRUE)) {
  x <- read.delim(f)
  cat(sprintf("%-46s fits %d err %d nonzero code %d ll_max %.6g ll_median %.5g ldisp_max %.1f max|ll - ref| %.3g\n", basename(f), sum(is.finite(x$ll)), sum(nzchar(as.character(x$err)) & !is.na(x$err)), sum(x$code != 0, na.rm = TRUE), max(x$ll, na.rm = TRUE), median(x$ll, na.rm = TRUE), max(x$ldisp, na.rm = TRUE), suppressWarnings(max(abs(x$ll - x$ref), na.rm = TRUE))))
}
a <- read.delim("dev/review-nlfix-latres-merge-none-bb.tsv"); b <- read.delim("dev/review-nlfix-latres-merge-box-bb.tsv")
cat("merge: none against box at 15, largest |ll difference|:", max(abs(a$ll - b$ll), na.rm = TRUE), " seeds ", paste(a$seed[abs(a$ll - b$ll) > 1e-3], collapse = ","), "\n")
