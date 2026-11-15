# Lane predraw2, backlog 4911: 3269's zib, zib zi ~ x (A1) and
# hurdle_poisson hu ~ x (A4) at 4000 draws, seed 1: mean, P(y = 0), 90%
# and 98% widths, the side law's message, the draw check's verdict, and
# the parameters the standard-error check holds (variance 0 in the draw
# covariance). brms 2.23.0 (dev/ppdraw-fix-brms.R, lane ppdraw's log):
# A1 5.7036 / 0.0098 / 4.9707 / 7.1353; A4 12.3835 / 0.0033 / 11.6305 /
# 16.7241; 3269 5.7022 / 0.0097 / 4.9752 / 7.0005.
# Usage: Rscript dev/predraw2-a14.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
source("dev/predraw2-cases.R")
for (id in c("zib3269", "A1", "A4")) {
  fit <- case_fit(id)$fit
  said <- character(0)
  set.seed(1)
  p <- withCallingHandlers(predict(fit, ndraws = 4000, summary = FALSE),
    message = function(m) {
      said <<- c(said, class(m)[1])
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      said <<- c(said, class(w)[1])
      invokeRestart("muffleWarning")
    })
  w <- function(a, b) mean(apply(p, 2, stats::quantile, b) -
                             apply(p, 2, stats::quantile, a))
  ds <- ns$fit_draw_space(fit)
  held <- rownames(ds$V)[diag(ds$V) == 0]
  chk <- fit$cache[["draw_check"]]
  cat(sprintf(paste("%-8s mean %.4f P(y=0) %.4f w90 %.4f w98 %.4f fitted",
                    "%.4f | check excess %s fire %s | held in V: %s |",
                    "said: %s\n"), id, mean(p), mean(p == 0), w(0.05, 0.95),
              w(0.01, 0.99), mean(suppressWarnings(fitted(fit))[,
                                                               "Estimate"]),
              format(chk$excess, digits = 3), format(chk$fire),
              paste(held, collapse = " "), paste(unique(said),
                                                 collapse = " ")))
}
