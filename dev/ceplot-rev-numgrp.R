# Reviewer check (lane ceplot): crossed (1 | g) + (1 | h) + (1 | g:h)
# with INTEGER grouping columns against the same model with factor
# columns. The fits are the same model; every display should agree.
#   Rscript dev/ceplot-rev-numgrp.R lane|base
# Data seed 49 (the worker's crossed case), boot seed 5.
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib",
          "C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
f4 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
set.seed(49)
dc <- expand.grid(g = factor(1:6), h = factor(1:5), r = 1:4)
dc <- dc[!(dc$g == "1" & dc$h == "1"), ]
dc$x <- rnorm(nrow(dc))
dc$y <- rnorm(nrow(dc), 1 + 0.5 * dc$x + rnorm(6)[dc$g] + rnorm(5)[dc$h] +
                rnorm(30, 0, 0.7)[as.integer(interaction(dc$g, dc$h))], 0.5)
di <- dc
di$g <- as.integer(as.character(di$g))
di$h <- as.integer(as.character(di$h))
fml <- bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h))
ff <- frm(fml, family = gaussian(), data = dc)
fi <- frm(fml, family = gaussian(), data = di)
cat("logLik factor", logLik(ff), "integer", logLik(fi), "\n")
run <- function(fit, cond, ...) {
  r <- tryCatch(withCallingHandlers(
    conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                        conditions = cond, ...)$x,
    warning = function(w) {
      cat("   WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(r, "error")) return(paste("ERROR", conditionMessage(r)))
  paste("est", f4(r$estimate__), "| width", f4(r$upper__ - r$lower__))
}
for (cc in list(c(2, 5), c(1, 1), c(3, 2))) {
  for (band in c("wald", "boot")) {
    args <- if (band == "boot") list(band = "boot", boot = 40, seed = 5) else
      list()
    cf <- list(g = as.character(cc[1]), h = as.character(cc[2]))
    ci <- list(g = as.integer(cc[1]), h = as.integer(cc[2]))
    cat(sprintf("g=%d h=%d %s\n  factor : %s\n  integer: %s\n", cc[1], cc[2],
                band, do.call(run, c(list(ff, cf), args)),
                do.call(run, c(list(fi, ci), args))))
  }
}
nd <- data.frame(x = c(0, 1, 2), g = c(2, 3, 4), h = c(5, 2, 1))
cat("fitted newdata (2,5),(3,2),(4,1), factor:\n")
ndf <- nd
ndf$g <- factor(ndf$g, levels = 1:6)
ndf$h <- factor(ndf$h, levels = 1:5)
print(fitted(ff, newdata = ndf)[, 1:2])
cat("integer:\n")
print(tryCatch(fitted(fi, newdata = nd)[, 1:2],
               error = function(e) conditionMessage(e)))
