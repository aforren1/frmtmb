# Puts require_fit_packages(<fit>) at the head of the frmtmb_fit
# methods that reached one part of a model on their own, or answered
# with NA, on a fit with a missing package
# (dev/smallfix-methods-probe.R; review finding 1). Run once.
#   Rscript dev/smallfix-head-checks.R
m <- list(
  c("R/confint.R", "hypothesis.frmtmb_fit", "x"),
  c("R/predict.R", "residuals.frmtmb_fit", "object"),
  c("R/confint.R", "drop1.frmtmb_fit", "object"),
  c("R/methods-fit.R", "fixef.frmtmb_fit", "object"),
  c("R/sugar.R", "ngrps.frmtmb_fit", "object"),
  c("R/confint.R", "parnames.frmtmb_fit", "x"),
  c("R/confint.R", "variables.frmtmb_fit", "x"),
  c("R/influence.R", "influence.frmtmb_fit", "model"),
  c("R/confint.R", "profile.frmtmb_fit", "fitted"),
  c("R/conditional-smooths.R", "conditional_smooths.frmtmb_fit", "x"),
  c("R/brms-postfit.R", "pp_mixture.frmtmb_fit", "x"),
  c("R/conditional-effects.R", "plot.frmtmb_fit", "x"),
  c("R/interop.R", "emm_basis.frmtmb_fit", "object"),
  c("R/interop.R", "recover_data.frmtmb_fit", "object"))
for (r in m) {
  x <- readLines(r[1])
  i <- grep(paste0("^", gsub(".", "[.]", r[2], fixed = TRUE),
                   " <- function"), x)
  stopifnot(length(i) == 1L)
  j <- i
  while (!grepl("[{]\\s*$", x[j])) j <- j + 1L
  x <- append(x, paste0("  require_fit_packages(", r[3], ")"), j)
  writeLines(x, r[1])
  cat(r[1], r[2], "line", j + 1L, "\n")
}
