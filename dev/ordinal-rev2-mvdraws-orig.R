# Reviewer: frmtmb.sample on multivariate ordinal fits. Which structures
# fail, where, and whether the flexible case fails on the base build too.
# Usage: Rscript ... <base|lane>. Data seed 20261009, sampler seed 5.
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "\n")
set.seed(20261009)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- 1L + (stats::rlogis(n) + d$x > -1) + (stats::rlogis(n) + d$x > 0) +
  (stats::rlogis(n) + d$x > 1)
d$y2 <- 1L + (stats::rlogis(n) + d$x > 0) + (stats::rlogis(n) + d$x > 1)
d$yg <- rnorm(n) + d$x
cases <- list(
  flex_flex = quote(list(cumulative(), sratio())),
  equi_flex = quote(list(cumulative(threshold = "equidistant"), sratio())),
  flex_stz = quote(list(cumulative(), sratio(threshold = "sum_to_zero"))),
  equi_stz = quote(list(cumulative(threshold = "equidistant"),
                        sratio(threshold = "sum_to_zero"))))
for (nm in names(cases)) {
  cat("\n==", nm, "\n")
  fam <- tryCatch(eval(cases[[nm]]), error = function(e) NULL)
  if (is.null(fam)) { cat("  family refused on this arm\n"); next }
  r <- tryCatch({
    fit <- frm(bf(y ~ x) + bf(y2 ~ x), data = d, family = fam)
    ds <- suppressWarnings(suppressMessages(
      frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 5)))
    v <- variables(ds)
    cat("  variables:", v[!grepl("^(lp__|lprior)", v)], "\n")
    ep <- posterior_epred(ds, ndraws = 2, resp = "y")
    cat("  epred(resp = y) dim:", dim(ep), "\n")
    "ok"
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n  call:", deparse(conditionCall(e))[1], "\n")
    tb <- sys.calls()
  })
}
