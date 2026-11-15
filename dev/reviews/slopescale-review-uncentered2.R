# Review of lane slopescale, follow-up of slopescale-review-uncentered.R:
# why does autoscale = TRUE reach the maximum on an uncentered slope
# column (x = 5 + 0.01 z) when the plan's factor is near 5? Compare, on
# cases the default leaves short: the plain fit, the plain fit restarted
# from its own stop, the pre-fit alone (the scaled frame fitted by
# itself), and autoscale = TRUE. Lane library only.
# Usage: Rscript slopescale-review-uncentered2.R
.libPaths(c("C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
mk <- function(seed, s, c0) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  z <- rnorm(300)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * z + 0.3 * z
  data.frame(g = g, x = c0 + s * z, y = eta + rnorm(300))
}
q <- function(expr) suppressMessages(suppressWarnings(expr))
cases <- list(list(1, 0.01, 5, y ~ 1 + (1 + x | g)),
              list(2, 0.03, 5, y ~ x + (1 + x | g)),
              list(2, 0.03, 5, y ~ 1 + (1 + x | g)),
              list(9, 0.01, 5, y ~ 1 + (1 + x | g)))
for (cs in cases) {
  d <- mk(cs[[1]], cs[[2]], cs[[3]])
  fo <- cs[[4]]
  fF <- q(frm(fo, data = d, control = frmtmb_control(autoscale = FALSE)))
  fR <- q(frm(fo, data = d, start = fF$estimates,
              control = frmtmb_control(autoscale = FALSE)))
  fT <- q(frm(fo, data = d, control = frmtmb_control(autoscale = TRUE)))
  fr <- frm(fo, data = d, dry_run = "frame")
  pl <- frmtmb:::autoscale_plan(fr)
  zs <- unlist(lapply(pl, function(p) vapply(p$z, `[[`, 0, "scale")))
  # the pre-fit alone: the data with the slope's column divided by the
  # plan's factor (and the fixed column centered and scaled, as the plan
  # does, when there is one)
  dsc <- d
  dsc$x <- d$x / zs[1]
  fP <- q(frm(fo, data = dsc, control = frmtmb_control(autoscale = FALSE)))
  # the column centered (what no plan does for a Z column)
  dc <- d
  dc$x <- (d$x - mean(d$x)) / sd(d$x)
  fC <- q(frm(fo, data = dc, control = frmtmb_control(autoscale = FALSE)))
  cat(sprintf("seed %d s %g %s: z factor %.4g | plain %.4f restart %.4f prefit-alone %.4f TRUE %.4f centered %.4f\n",
              cs[[1]], cs[[2]], deparse(fo), zs[1], logLik(fF), logLik(fR),
              logLik(fP), logLik(fT), logLik(fC)))
  cat("   plain theta", signif(fF$opt$par[names(fF$opt$par) == "theta"], 4),
      "; TRUE theta", signif(fT$opt$par[names(fT$opt$par) == "theta"], 4),
      "\n")
}
