# Lane slopetrig, item 5033: calibration grid for the ratio of root mean
# square to sd of a random-slope column (over its nonzero rows, as
# autoscale_block_cols() reads it).
#
# Designs (z ~ N(0, 1), trt alternating TRUE/FALSE):
#   unc     x = c0 + r c0 z, c0 in 0.1 .. 1000, r in 1e-3 .. 1
#   year    x = 2000 + s z, s in 1, 3, 10, 30
#   dose    x = trt (1 + r z), r in 0.01 .. 0.3 (indicator times a
#           near-constant level), and the short check's twonz and ulpind
#   binlev  x = c0 + trt + 0.01 z (a 0/1 column on a level)
# eta = 0.5 + b0[g] + b1[g] z + 0.3 z with b0, b1 ~ N(0, 0.5) (unc, year,
# binlev), or 0.5 + b0[g] + (0.3 + b1[g]) x (dose, as
# dev/reviews/slopescale-short-trigger.R); gaussian (noise sd 1) or
# poisson; 20 groups of 15; y ~ 1 + (1 + x | g) and y ~ x + (1 + x | g);
# seeds 1 to 5.
#
# Arms per case, interleaved in one process: ref (the same model on the
# centered, standardized column, default control), F (autoscale =
# FALSE), T (autoscale = TRUE), D (the default of the library loaded).
# Each arm's log-likelihood, warning count, convergence code and elapsed
# time. "short" is read later against the best arm.
# Usage: Rscript slopetrig-calib.R <gau|pois> <base|lane> <out.tsv>
a <- commandArgs(TRUE)
fam <- a[1]
arm <- a[2]
out <- a[3]
stopifnot(fam %in% c("gau", "pois"), arm %in% c("base", "lane"))
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "fam", fam, "frmtmb", as.character(packageVersion("frmtmb")),
    "lib:", dirname(find.package("frmtmb")), "\n")
cells <- list()
for (c0 in c(0.1, 1, 10, 100, 1000))
  for (r in c(1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1))
    cells[[length(cells) + 1L]] <- list(design = "unc", c0 = c0, r = r)
for (s in c(1, 3, 10, 30))
  cells[[length(cells) + 1L]] <- list(design = "year", c0 = 2000, r = s / 2000)
for (r in c(0.01, 0.03, 0.1, 0.3))
  cells[[length(cells) + 1L]] <- list(design = "dose", c0 = 1, r = r)
cells[[length(cells) + 1L]] <- list(design = "twonz", c0 = 1, r = NA)
cells[[length(cells) + 1L]] <- list(design = "ulpind", c0 = 0.3, r = NA)
for (c0 in c(1, 10, 100))
  cells[[length(cells) + 1L]] <- list(design = "binlev", c0 = c0, r = NA)
mk <- function(cell, seed) {
  set.seed(seed)
  n <- 300
  g <- factor(rep(1:20, each = 15))
  trt <- rep(c(TRUE, FALSE), length.out = n)
  z <- rnorm(n)
  b0 <- rnorm(20, 0, 0.5)[g]
  b1 <- rnorm(20, 0, 0.5)[g]
  c0 <- cell$c0
  r <- cell$r
  x <- switch(cell$design,
    unc = c0 + r * c0 * z,
    year = c0 + r * c0 * z,
    dose = ifelse(trt, 1 + r * z, 0),
    twonz = ifelse(trt, sample(c(1, 1.01), n, TRUE), 0),
    ulpind = ifelse(trt, sample(c(0.3, 0.1 * 3), n, TRUE), 0),
    binlev = c0 + trt + 0.01 * z)
  eta <- if (cell$design %in% c("dose", "twonz", "ulpind")) {
    0.5 + b0 + (0.3 + b1) * x
  } else if (cell$design == "binlev") {
    0.5 + b0 + (0.3 + b1) * (x - c0)
  } else {
    0.5 + b0 + b1 * z + 0.3 * z
  }
  y <- if (fam == "gau") eta + rnorm(n) else rpois(n, exp(eta))
  data.frame(g = g, x = x, y = y)
}
run1 <- function(expr) {
  nw <- 0L
  t0 <- proc.time()[["elapsed"]]
  f <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(w) {
      nw <<- nw + 1L
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(f, "error")) {
    return(list(ll = NA_real_, nw = nw, conv = NA_integer_, t = el,
                as = NA))
  }
  list(ll = as.numeric(logLik(f)), nw = nw,
       conv = as.integer(f$opt$convergence), t = el,
       as = isTRUE(f$autoscaled))
}
fm <- if (fam == "gau") gaussian() else poisson()
hdr <- TRUE
for (cell in cells) for (fo in c("y ~ 1 + (1 + x | g)", "y ~ x + (1 + x | g)"))
  for (seed in 1:5) {
    d <- mk(cell, seed)
    nz <- d$x[d$x != 0]
    fml <- stats::as.formula(fo)
    dc <- d
    dc$x <- (d$x - mean(d$x)) / sd(d$x)
    ref <- run1(frm(fml, family = fm, data = dc))
    rF <- run1(frm(fml, family = fm, data = d,
                   control = frmtmb_control(autoscale = FALSE)))
    rT <- run1(frm(fml, family = fm, data = d,
                   control = frmtmb_control(autoscale = TRUE)))
    rD <- run1(frm(fml, family = fm, data = d))
    row <- data.frame(
      design = cell$design, c0 = cell$c0, r = cell$r, fam = fam,
      form = if (grepl("^y ~ x", fo)) "x" else "1", seed = seed,
      ndist = length(unique(d$x)), spread = sd(nz),
      rms = sqrt(mean(nz^2)), ratio = sqrt(mean(nz^2)) / sd(nz),
      ll_ref = ref$ll, ll_F = rF$ll, ll_T = rT$ll, ll_D = rD$ll,
      w_ref = ref$nw, w_F = rF$nw, w_T = rT$nw, w_D = rD$nw,
      c_F = rF$conv, c_T = rT$conv, c_D = rD$conv, as_D = rD$as,
      t_ref = ref$t, t_F = rF$t, t_T = rT$t, t_D = rD$t)
    utils::write.table(row, out, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = hdr, append = !hdr)
    hdr <- FALSE
  }
cat("DONE", length(cells), "cells\n")
