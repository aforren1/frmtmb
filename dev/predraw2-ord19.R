# Lane predraw2, fix round (review FIX 2): the 19 fired ordinal-mixture
# fits of the review (dev/reviews/predraw2-review-e.R, seeds below), on
# this build, whose importance weights carry the side law's log Jacobian
# itself. The mean absolute error of the predicted category shares
# against the true probabilities, 1000 draws, seed 1; the review's
# figures for the anchored weights are in predraw2-review-e.txt ("mae
# lane").
# Usage: Rscript dev/predraw2-ord19.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
n <- 300L
source("dev/ppdraw-sweep-gen.R")
cuts <- c(-Inf, -1, 0.5, 2, Inf)
truthP <- function(x) {
  0.6 * t(vapply(x, function(xi) diff(stats::plogis(cuts - xi)),
                 numeric(4))) +
    0.4 * t(vapply(x, function(xi) diff(stats::plogis(cuts - xi - 3)),
                   numeric(4)))
}
seeds <- c(3, 4, 8, 9, 10, 11, 13, 15, 18, 25, 26, 28, 29, 31, 32, 35, 36,
           39, 40)
for (seed in seeds) {
  d <- gen("mixture_cumulative", "interior", seed)
  s <- spec("mixture_cumulative")
  fit <- tryCatch(suppressMessages(suppressWarnings(
    frm(s[[1]], family = s[[2]], data = d))), error = function(e) NULL)
  if (is.null(fit)) {
    cat("seed", seed, "fit failed\n")
    next
  }
  set.seed(1)
  p <- suppressWarnings(suppressMessages(
    predict(fit, ndraws = 1000L, summary = FALSE)))
  sh <- t(apply(p, 2L, function(col) tabulate(col, 4L) / 1000))
  chk <- fit$cache[["draw_check"]]
  law <- fit$cache[["draw_edge_law"]]
  cat(sprintf("seed %2d fire %-5s excess %8.3f rel %6.3f narrowed %d | mae %.4f\n",
              seed, chk$fire, chk$excess, chk$rel,
              if (is.list(law)) length(law$narrowed) else 0L,
              mean(abs(sh - truthP(d$x)))))
}
