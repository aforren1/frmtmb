# Lane slopetrig: the calibration grids (dev/slopetrig-calib.R and
# dev/slopetrig-calib2.R, base arm) as one table of short fits under
# autoscale = FALSE, among the columns the spread rule does not engage
# (spread >= 0.05): rows are the column's root mean square (binned),
# columns the ratio of root mean square to sd (binned); each entry is
# "short / fits". Then the same for autoscale = TRUE.
# Usage: Rscript slopetrig-calib-tab.R <tsv> [<tsv> ...]
fs <- commandArgs(TRUE)
d <- do.call(rbind, lapply(fs, utils::read.delim))
best <- pmax(d$ll_ref, d$ll_F, d$ll_T, d$ll_D, na.rm = TRUE)
d$short_F <- best - d$ll_F > 1e-3
d$short_T <- best - d$ll_T > 1e-3
d <- d[d$spread >= 0.05, ]
d$rms_bin <- cut(d$rms, c(0, 2, 5, 20, 200, Inf),
                 labels = c("rms<2", "2-5", "5-20", "20-200", ">200"))
d$ratio_bin <- cut(d$ratio, c(1, 1.8, 2.5, 3.5, 4.5, 6, 9, 15, 40, Inf),
                   right = FALSE)
tab <- function(v) {
  n <- table(d$rms_bin, d$ratio_bin)
  s <- tapply(v, list(d$rms_bin, d$ratio_bin), sum)
  s[is.na(s)] <- 0
  out <- matrix(paste0(s, "/", n), nrow(n), dimnames = dimnames(n))
  out[n == 0] <- "."
  print(noquote(out))
}
options(width = 200)
cat("fits with spread >= 0.05:", nrow(d), "\nshort under autoscale = FALSE:\n")
tab(d$short_F)
cat("\nshort under autoscale = TRUE:\n")
tab(d$short_T)
cat("\nby family, short F:", tapply(d$short_F, d$fam, sum), "of",
    tapply(d$short_F, d$fam, length), "\n")
