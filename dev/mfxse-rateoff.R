# The rows dev/mfxse-rate.R counts as material on one arm, with the design,
# the estimate and the reference, to tell a delta-method limit from a defect.
#   Rscript dev/mfxse-rateoff.R <arm>
arm <- commandArgs(TRUE)[1]
O <- "C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-rate-out/"
fs <- list.files(O, pattern = paste0("^", arm, "-.*tsv$"), full.names = TRUE)
r <- do.call(rbind, lapply(fs, read.delim))
x <- r[is.finite(r$ratio) & (r$ratio < 0.8 | r$ratio > 1.25), ]
x$fit_ratio <- x$fitted_se / x$ref
x$mfx_fit <- x$me_se / x$fitted_se
print(x[, c("family", "J", "m", "sdg", "slope", "seed", "qty", "est", "me_se",
            "ref", "ratio", "fit_ratio", "mfx_fit", "conv")], row.names = FALSE,
      digits = 4)
cat("\nmaterial rows", nrow(x), "on", nrow(unique(x[, c("family", "J", "m",
                                                       "sdg", "slope",
                                                       "seed")])),
    "fits; with conv != 0:", sum(x$conv != 0), "\n")
