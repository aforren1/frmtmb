# Reviewer probe: does frmtmb.sample distinguish `sigma ~ 1` from a
# scalar sigma (it would need to for `xi ~ 1`)?
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1)
d <- data.frame(y = rnorm(100), x = rnorm(100))
for (f in list(bf(y ~ x), bf(y ~ x, sigma ~ 1))) {
  cat("\n", deparse(f$formula), if (length(f$pforms)) deparse(f$pforms[[1]]), "\n")
  print(as.data.frame(default_prior(f, data = d, family = gaussian(),
                                    route = "sample"))[, c("prior", "class", "dpar")])
  print(as.data.frame(suppressWarnings(brms::default_prior(
    if (length(f$pforms)) brms::bf(y ~ x, sigma ~ 1) else brms::bf(y ~ x),
    data = d, family = gaussian())))[, c("prior", "class", "dpar")])
}
fit <- frm(bf(y ~ x, xi ~ 1), family = gen_extreme_value(),
           data = transform(d, y = y + 3))
lp <- Filter(function(l) identical(l$dpar, "xi"), fit$frame$linpreds)[[1]]
cat("xi ~ 1 linpred: ncol(X) =", ncol(lp$X), "; names:", names(lp), "\n")
