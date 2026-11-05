# Item 7, control for the sampling failure: same call on mi() alone and cens() alone.
source("dev/review-refuse-lib.R")
set.seed(10004); n <- 180; x <- rnorm(n)
base <- data.frame(x = x, y = rnorm(n, 0.5 * x), ce = 0L)
base$y[seq(5, n, 9)] <- NA
for (lab in c("mi", "cens", "mi+cens")) {
  form <- switch(lab, mi = brms::bf(y | mi() ~ x), cens = brms::bf(y | cens(ce) ~ x), `mi+cens` = brms::bf(y | mi() + cens(ce) ~ x))
  dd <- if (lab == "cens") base[!is.na(base$y), ] else base
  pri <- E$brms_flat_prior(form, data = dd, family = stats::gaussian())
  code <- brms::make_stancode(form, data = dd, prior = pri)
  sdat <- E$brms_standata(form, data = dd, prior = pri)
  mod <- E$brms_stan_model(code)
  s <- tr(suppressWarnings(rstan::sampling(mod, data = sdat, chains = 1, iter = 200, refresh = 0, seed = 1)))
  pr("%-8s draws kept %s", lab, if (is_err(s)) msg(s, 100) else sum(s@sim$n_save))
  s2 <- tr(suppressWarnings(rstan::sampling(rstan::stan_model(model_code = code, auto_write = FALSE), data = sdat, chains = 1, iter = 200, refresh = 0, seed = 1)))
  pr("%-8s fresh compile, draws kept %s", lab, if (is_err(s2)) msg(s2, 100) else sum(s2@sim$n_save))
}
