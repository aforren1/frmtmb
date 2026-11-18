# Lane nlband: the band of poles of 5115, seeds 1 to 3 (the review's) and
# 20 fresh seeds 101 to 120. For each seed: the default fit, and the fit
# from the review's start (1, -1.9) with the pole in the gap of the data.
# NLBAND_ARM=before|after Rscript dev/nlband-repro.R [seeds]
args <- commandArgs(trailingOnly = TRUE)
seeds <- if (length(args)) eval(parse(text = args[1])) else c(1:3, 101:120)
source("dev/nlband-common.R")
for (s in seeds) {
  d <- band_data(s)
  t0 <- proc.time()[[3]]
  r <- tt(frm(band_form, data = d))
  el <- proc.time()[[3]] - t0
  g <- tt(frm(band_form, data = d, start = list(beta = c(1, -1.9))))
  b <- if (iserr(r)) c(NA, NA) else r$estimates$beta
  bg <- if (iserr(g)) c(NA, NA) else g$estimates$beta
  cat(sprintf(paste0("ROW seed %d n %d default %.6f code %s b %.4f a %.4f ",
                     "warn %s msg %s t %.2f | gap-start %.6f code %s b %.4f ",
                     "warn %s msg %s | short %.3f\n"),
              s, nrow(d), llv(r), if (iserr(r)) "ERR" else r$opt$convergence,
              b[2], b[1], wsum(r), msum(r), el, llv(g),
              if (iserr(g)) "ERR" else g$opt$convergence, bg[2], wsum(g),
              msum(g), llv(g) - llv(r)))
}
cat("DONE\n")
