# Lane `learnhier`: print one replicate record.
#   Rscript dev/learnhier-peek.R <rds>
source("dev/learnhier-env.R")
a <- commandArgs(trailingOnly = TRUE)
r <- readRDS(a[[1L]])
cat("design", r$design, "seed", r$seed, r$ns, "x", r$nt,
    "rows", r$rows, "\n")
cat("ok", r$ok, "fit_s", round(r$fit_s, 1), "conv", r$conv,
    "maxgrad", format(r$max_grad, digits = 3), "pdHess", r$pdHess,
    "nbadse", r$n_bad_se, "npar", r$n_par, "\n")
cat("logLik", round(r$logLik, 3), " oracle", round(r$oracle, 3),
    " diff", round(r$logLik - r$oracle, 4), "\n")
cat("peak R heap MB", round(r$peak_r_mb), "\n\n")
cat("fixed:\n"); print(round(r$fixed, 4))
cat("\nvarcorr:\n"); print(round(r$vc, 4))
cat("\ntruth, drawn block:\n"); print(round(r$truth_drawn, 4))
cat("\ntruth, fitted block:\n"); print(round(r$truth_fitted, 4))
if (!is.null(r$ndt)) {
  cat("\nndt:\n"); print(round(r$ndt, 4))
}
