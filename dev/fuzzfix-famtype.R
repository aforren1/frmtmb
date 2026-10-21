# dpar names and frame linpreds of the families the refusal names
.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
d <- data.frame(y = c(0.2, 0.5, 0.7, 0.4), x = 1:4)
for (f in list(gaussian(), student(), lognormal(), Gamma(), weibull(), Beta(),
               inverse.gaussian(), skew_normal(), shifted_lognormal(), huber(),
               von_mises(), frechet(), gen_extreme_value())) {
  sp <- frm(bf(y ~ x), family = f, data = d, dry_run = "spec")
  fr <- frm(bf(y ~ x), family = f, data = d, dry_run = "frame")
  fm <- sp$responses[[1]]$family
  cat(sprintf("%-18s dpars=%s linpreds=%s\n", fm$family,
              paste(fm$dpars, collapse = ","), paste(names(fr$linpreds), collapse = ",")))
}
fr <- frm(bf(y ~ x, sigma = 2), family = gaussian(), data = d, dry_run = "frame")
str(fr$linpreds$y.sigma$constant)
