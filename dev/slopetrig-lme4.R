# Lane slopetrig, item 5033: frmtmb's default against lme4 and glmmTMB
# on the item's poorly centered designs (x = 2000 + 3 z, x = 5 + 0.1 z;
# dev/slopetrig-repro.R's data), gaussian, both forms, seeds 1 to 5, on
# the raw column (_raw) and on the centered, standardized column (_c).
# ML throughout (REML = FALSE). Each entry is the shortfall against the
# best of the six fits of the same data.
# Usage: Rscript slopetrig-lme4.R base|lane
arm <- commandArgs(TRUE)[1]
stopifnot(arm %in% c("base", "lane"))
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(lme4)
  library(glmmTMB)
})
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")),
    "lib:", dirname(find.package("frmtmb")),
    "lme4", as.character(packageVersion("lme4")),
    "glmmTMB", as.character(packageVersion("glmmTMB")), "\n")
mk <- function(seed, s, c0) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  z <- rnorm(300)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * z + 0.3 * z
  data.frame(g = g, x = c0 + s * z, y = eta + rnorm(300))
}
q <- function(expr) {
  r <- try(suppressWarnings(suppressMessages(expr)), silent = TRUE)
  if (inherits(r, "try-error")) NA_real_ else as.numeric(logLik(r))
}
rows <- list()
for (cc in list(c(2000, 3), c(5, 0.1)))
  for (fo in c("y ~ 1 + (1 + x | g)", "y ~ x + (1 + x | g)"))
    for (seed in 1:5) {
      d <- mk(seed, cc[2], cc[1])
      dc <- d
      dc$x <- (d$x - mean(d$x)) / sd(d$x)
      fml <- stats::as.formula(fo)
      one <- function(dd) {
        c(frmtmb = q(frm(fml, data = dd)),
          lme4 = q(lmer(fml, data = dd, REML = FALSE)),
          glmmTMB = q(glmmTMB(fml, data = dd, REML = FALSE)))
      }
      a <- one(d)
      b <- one(dc)
      best <- max(c(a, b), na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        c0 = cc[1], s = cc[2], form = fo, seed = seed,
        frm_raw = best - a[["frmtmb"]], frm_c = best - b[["frmtmb"]],
        lme4_raw = best - a[["lme4"]], lme4_c = best - b[["lme4"]],
        tmb_raw = best - a[["glmmTMB"]], tmb_c = best - b[["glmmTMB"]])
    }
r <- do.call(rbind, rows)
options(width = 200)
print(r, digits = 3, row.names = FALSE)
sh <- function(v) sprintf("%d of %d short > 1e-6 (max %.3g, NA %d)",
                          sum(v > 1e-6, na.rm = TRUE), length(v),
                          max(v, na.rm = TRUE), sum(is.na(v)))
cat("\nshort of the best fit, by column:\n")
for (nm in names(r)[5:10]) cat(" ", nm, ":", sh(r[[nm]]), "\n")
cat("max |frm_raw - lme4_c|:", signif(max(abs(r$frm_raw - r$lme4_c)), 3),
    "; max |frm_raw - tmb_c|:",
    signif(max(abs(r$frm_raw - r$tmb_c), na.rm = TRUE), 3), "\n")
