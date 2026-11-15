# Lane slopescale, item 1018, brief step 4: frmtmb against lme4 and
# glmmTMB on the triage's design, with the slope's column at spread s
# (1e-4, 1e-5) and at spread 1. Each optimum is the log-likelihood the
# package reports (ML: REML = FALSE for lmer and glmmTMB). The reference
# is the best of all fits of the same data on either scale.
# Usage: Rscript slopescale-lme4.R before|after
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(lme4)
  library(glmmTMB)
})
cat("frmtmb", as.character(packageVersion("frmtmb")), find.package("frmtmb"),
    "lme4", as.character(packageVersion("lme4")),
    "glmmTMB", as.character(packageVersion("glmmTMB")), "arm", arm, "\n")
q <- function(expr) {
  r <- try(suppressWarnings(suppressMessages(expr)), silent = TRUE)
  if (inherits(r, "try-error")) NA_real_ else as.numeric(logLik(r))
}
rows <- list()
for (fam in c("gau", "pois")) for (form in c("0", "1"))
  for (s in c(1e-4, 1e-5)) for (seed in 1:6) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x0 <- rnorm(300)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * x0
  y <- if (fam == "gau") eta + rnorm(300) else rpois(300, exp(eta))
  d_s <- data.frame(g = g, x = x0 * s, y = y)
  d_1 <- data.frame(g = g, x = x0, y = y)
  fo <- if (form == "0") y ~ 1 + (0 + x | g) else y ~ 1 + (1 + x | g)
  one <- function(d) {
    c(frmtmb = q(frm(fo, family = if (fam == "gau") gaussian() else
                       poisson(), data = d)),
      lme4 = if (fam == "gau") q(lmer(fo, data = d, REML = FALSE)) else
        q(glmer(fo, data = d, family = poisson())),
      glmmTMB = q(glmmTMB(fo, data = d, REML = FALSE,
                          family = if (fam == "gau") gaussian() else
                            poisson())))
  }
  a <- one(d_s)
  b <- one(d_1)
  best <- max(c(a, b), na.rm = TRUE)
  rows[[length(rows) + 1L]] <- data.frame(
    fam = fam, form = form, s = s, seed = seed,
    frm_s = best - a[["frmtmb"]], frm_1 = best - b[["frmtmb"]],
    lme4_s = best - a[["lme4"]], lme4_1 = best - b[["lme4"]],
    tmb_s = best - a[["glmmTMB"]], tmb_1 = best - b[["glmmTMB"]])
}
r <- do.call(rbind, rows)
options(width = 200)
print(r, digits = 3)
sh <- function(v) sprintf("%d of %d short > 1e-6 (max %.3g, NA %d)",
                          sum(v > 1e-6, na.rm = TRUE), length(v),
                          max(v, na.rm = TRUE), sum(is.na(v)))
cat("\nshort of the best fit, by column:\n")
for (nm in names(r)[5:10]) cat(" ", nm, ":", sh(r[[nm]]), "\n")
