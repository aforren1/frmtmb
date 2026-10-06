# Lane fixes, item 4e: what variables() lists for disc, frmtmb's fit,
# frmtmb.sample's draws, and brms 2.23.0 (a fit made with empty = TRUE
# carries no draws, so brms's names come from its Stan program).
#   Rscript dev/fixes-disc-vars.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(62)
n <- 150
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8)
for (fo in list(y ~ x, bf(y ~ x, disc ~ 0 + z))) {
  cat("== frmtmb", deparse(if (inherits(fo, "formula")) fo else fo$formula),
      if (!inherits(fo, "formula")) "+ disc ~ 0 + z", "\n")
  fit <- frm(fo, family = cumulative(), data = d)
  cat("variables(fit):", variables(fit), "\n")
  ds <- suppressMessages(suppressWarnings(
    frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 1)))
  cat("variables(draws):", variables(ds), "\n")
  code <- brms::stancode(if (inherits(fo, "formula")) brms::bf(fo) else
    brms::bf(y ~ x, disc ~ 0 + z), data = d, family = brms::cumulative())
  cat("brms Stan lines naming disc:\n")
  cat(grep("disc", strsplit(code, "\n")[[1]], value = TRUE), sep = "\n")
}
