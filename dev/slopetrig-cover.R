# Lane slopetrig, fix round: does reading the ratio over the rows the
# block covers change any calibration design? Every grid design is one
# (1 + x | g) block, whose intercept is nonzero on every row, so it
# should not. For each cell of dev/slopetrig-calib.R and
# dev/slopetrig-calib2.R (gaussian, seed 1, both forms), the plan's
# ratio from frm(dry_run = "frame") against the ratio over every row,
# and whether the default's decision changes. Also the dose on half the
# rows at r = 0.2 under (1 + x | g) and under (0 + x | g).
# Usage: Rscript slopetrig-cover.R
.libPaths(c("C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
cells <- list()
for (c0 in c(0.1, 1, 10, 100, 1000))
  for (r in c(1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1))
    cells[[length(cells) + 1L]] <- list(design = "unc", c0 = c0, r = r)
for (s in c(1, 3, 10, 30))
  cells[[length(cells) + 1L]] <- list(design = "unc", c0 = 2000, r = s / 2000)
for (rho in c(2, 3, 4, 5, 7))
  for (c0 in c(1, 10, 1000))
    cells[[length(cells) + 1L]] <- list(design = "unc", c0 = c0,
                                        r = 1 / sqrt(rho^2 - 1))
for (r in c(0.01, 0.03, 0.1, 0.3, 0.2))
  cells[[length(cells) + 1L]] <- list(design = "dose", c0 = 1, r = r)
cells[[length(cells) + 1L]] <- list(design = "twonz", c0 = 1, r = NA)
cells[[length(cells) + 1L]] <- list(design = "ulpind", c0 = 0.3, r = NA)
for (c0 in c(0.5, 1, 2, 3, 10, 100))
  cells[[length(cells) + 1L]] <- list(design = "binlev", c0 = c0, r = NA)
mk <- function(cell) {
  set.seed(1)
  n <- 300
  g <- factor(rep(1:20, each = 15))
  trt <- rep(c(TRUE, FALSE), length.out = n)
  z <- rnorm(n)
  x <- switch(cell$design,
    unc = cell$c0 + cell$r * cell$c0 * z,
    dose = ifelse(trt, 1 + cell$r * z, 0),
    twonz = ifelse(trt, sample(c(1, 1.01), n, TRUE), 0),
    ulpind = ifelse(trt, sample(c(0.3, 0.1 * 3), n, TRUE), 0),
    binlev = cell$c0 + trt + 0.01 * z)
  data.frame(g = g, x = x, y = rnorm(n))
}
n <- 0L
nchg <- 0L
dmax <- 0
for (cell in cells) for (fo in c("y ~ 1 + (1 + x | g)", "y ~ x + (1 + x | g)")) {
  d <- mk(cell)
  fr <- frm(stats::as.formula(fo), data = d, dry_run = "frame")
  zz <- ns$autoscale_plan(fr)$y.mu$z
  if (!length(zz)) next
  all_rows <- sqrt(mean(d$x^2)) / sd(d$x)
  n <- n + 1L
  dmax <- max(dmax, abs(zz[[1]]$ratio / all_rows - 1))
  nchg <- nchg + ((zz[[1]]$ratio > 2) != (all_rows > 2))
}
cat("designs", n, "; max |plan ratio / every-row ratio - 1|", signif(dmax, 3),
    "; decisions changed", nchg, "\n")
for (fo in c("y ~ 1 + (1 + x | g)", "y ~ 1 + (0 + x | g)")) {
  d <- mk(list(design = "dose", c0 = 1, r = 0.2))
  fr <- frm(stats::as.formula(fo), data = d, dry_run = "frame")
  zz <- ns$autoscale_plan(fr)$y.mu$z[[1]]
  cat(sprintf("dose 0 or 1 + 0.2 z, %s: ratio %.3g (every row %.3g, nonzero rows %.3g)\n",
              fo, zz$ratio, sqrt(mean(d$x^2)) / sd(d$x),
              sqrt(mean(d$x[d$x != 0]^2)) / sd(d$x[d$x != 0])))
}
