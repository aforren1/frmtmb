# Area 4: D7, the cov = FALSE divergence. Seed 1582 as dev/nlfix-arnl2.R, plus
# Stan's log_prob for each spelling, a by-hand density, the poisson latent form
# (brms-45) and what simulate() draws.
source("dev/review-nlfix-common.R")
suppressPackageStartupMessages(library(testthat))
E <- new.env(parent = asNamespace("frmtmb"))
E$`%||%` <- function(a, b) if (is.null(a)) b else a
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) sys.source(f, envir = E)
d <- E$eq_dat(seed = 1582)
set.seed(1582)
d$y <- 2 + 0.8 * d$x + rnorm(nrow(d), 0, 0.5)
# real AR: y gets an AR(1) residual of 0.6 within each group
e <- ave(rnorm(nrow(d), 0, 0.5), d$g, FUN = function(z) as.numeric(stats::filter(z, 0.6, "recursive")))
d$y2 <- 2 + 0.8 * d$x + e
run <- function(yv) {
  cat("\n=== response", yv, "\n")
  d$y <- d[[yv]]
  fl <- frm(bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = d)
  fn <- suppressWarnings(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE,
                                autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = d,
                             start = list(beta = c("b0_(Intercept)" = 2, "b1_(Intercept)" = 0.5))))
  cat("frm logLik lin, nl:", format(c(as.numeric(logLik(fl)), as.numeric(logLik(fn))), digits = 12), "\n")
  b <- fn$estimates$beta
  bodyv <- b[[1]] + b[[2]] * d$x
  ar <- as.numeric(fn$estimates$thetaac %||% NA)
  print(fn$estimates[setdiff(names(fn$estimates), "b")])
  # by hand: brms's Err[n, 1] is the previous row's residual inside the group
  res <- d$y - bodyv
  lag <- ave(res, d$g, FUN = function(z) c(0, z[-length(z)]))
  sig <- as.numeric(sigma(fn))[1]
  phi <- as.numeric(unlist(fn$estimates[grep("ac|ar", names(fn$estimates))]))
  fx <- function(p) sum(dnorm(d$y, bodyv + p * lag, sig, log = TRUE))
  cv <- tt(confint_varcorr(fn))
  phi_hat <- optimize(fx, c(-0.99, 0.99), maximum = TRUE)$maximum
  cat("by hand phi at sigma:", phi_hat, "ll", fx(phi_hat), " ll with phi=0:", fx(0), "\n")
  bn <- brms::bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = FALSE))
  bl <- brms::bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = FALSE))
  fam <- brms::brmsfamily("gaussian")
  fxn <- E$brms_fixed_fit(bn, fam, d, fn)
  fxl <- E$brms_fixed_fit(bl, fam, d, fl)
  en <- brms::posterior_epred(fxn)[1, ]
  el <- brms::posterior_epred(fxl)[1, ]
  ln <- brms::log_lik(fxn)[1, ]
  ll <- brms::log_lik(fxl)[1, ]
  cat("brms nl epred - body, max:", max(abs(en - bodyv)), " brms lin epred - frm nl fitted, max:",
      max(abs(el - fitted(fn)[, "Estimate"])), "\n")
  cat("brms log_lik sum nl, lin:", format(c(sum(ln), sum(ll)), digits = 12), "\n")
  for (arm in c("nl", "lin")) {
    f <- if (arm == "nl") bn else bl
    fit <- if (arm == "nl") fn else fl
    prior <- E$brms_flat_prior(f, data = d, family = fam)
    sdat <- E$brms_standata(f, data = d, family = fam, prior = prior)
    code <- brms::make_stancode(f, data = d, family = fam, prior = prior)
    sf <- suppressMessages(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 0))
    pars <- E$stan_pars_from_fit(fit, sdat, code, NULL)
    lp <- rstan::log_prob(sf, rstan::unconstrain_pars(sf, pars), adjust_transform = FALSE)
    cat("Stan log_prob", arm, format(lp, digits = 12), "\n")
    if (arm == "nl") cat(grep("Err\\[|mu\\[n\\]|mu\\[", strsplit(as.character(code), "\n")[[1]], value = TRUE), sep = "\n")
  }
  # predictive: brms's and frm's, mean and variance per row at the same values
  pl <- brms::posterior_predict(fxl, ndraws = 1)
  set.seed(9)
  ppl <- do.call(rbind, lapply(1:2000, function(i) brms::posterior_predict(fxl, ndraws = 1)))
  ppn <- do.call(rbind, lapply(1:2000, function(i) brms::posterior_predict(fxn, ndraws = 1)))
  shift <- fitted(fn)[, "Estimate"] - bodyv
  cat("brms pp lin slope on shift:", coef(lm(I(colMeans(ppl) - (fixef(fl)[1, 1] + fixef(fl)[2, 1] * d$x)) ~ shift))[2],
      " row sd mean", mean(apply(ppl, 2, sd)), "\n")
  cat("brms pp nl slope on shift:", coef(lm(I(colMeans(ppn) - bodyv) ~ shift))[2], " row sd mean", mean(apply(ppn, 2, sd)), "\n")
  sl <- as.matrix(simulate(fl, nsim = 2000, seed = 3))
  sn <- as.matrix(simulate(fn, nsim = 2000, seed = 3))
  cat("frm sim lin slope on shift:", coef(lm(I(rowMeans(sl) - (fixef(fl)[1, 1] + fixef(fl)[2, 1] * d$x)) ~ shift))[2],
      " row sd mean", mean(apply(sl, 1, sd)), "\n")
  cat("frm sim nl slope on shift:", coef(lm(I(rowMeans(sn) - bodyv) ~ shift))[2], " row sd mean", mean(apply(sn, 1, sd)), "\n")
  cat("sigma:", sig, " marginal if AR:", sig / sqrt(1 - phi_hat^2), "\n")
  # within-group lag-1 correlation of a simulated series about the body
  rs <- sn - bodyv
  lagc <- mean(sapply(1:200, function(j) {
    r <- rs[, j]
    l <- ave(r, d$g, FUN = function(z) c(NA, z[-length(z)]))
    cor(r, l, use = "complete.obs")
  }))
  cat("frm simulate within-group lag-1 cor about the body:", lagc, " fitted phi:", phi_hat, "\n")
  ob <- ave(res, d$g, FUN = function(z) c(NA, z[-length(z)]))
  cat("data lag-1 cor about the body:", cor(res, ob, use = "complete.obs"), "\n")
  pfrm <- predict(fn, summary = FALSE, ndraws = 2000, propagate_error = FALSE)
  cat("frm predict nl slope on shift:", coef(lm(I(colMeans(pfrm) - bodyv) ~ shift))[2], " row sd mean", mean(apply(pfrm, 2, sd)), "\n")
  # pp_check-style: var(y sim about the mean) vs data
  cat("resid sd about one-step mean, data:", sd(d$y - fitted(fn)[, "Estimate"]), "\n")
  # refit from a simulated data set
  ds <- d; ds$y <- sn[, 1]
  fs <- suppressWarnings(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = FALSE)),
                             data = ds, start = list(beta = c("b0_(Intercept)" = 2, "b1_(Intercept)" = 0.5))))
  cat("refit on one simulate() replicate: ar estimate:", paste(format(unlist(fs$estimates["thetaac"]), digits = 5), collapse = ","),
      " original:", paste(format(unlist(fn$estimates["thetaac"]), digits = 5), collapse = ","), "\n")
  phis <- sapply(1:40, function(j) {
    ds$y <- sn[, j]
    r <- suppressWarnings(frm(bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = ds))
    v <- r$estimates[grep("ac", names(r$estimates))]
    as.numeric(unlist(v))[1]
  })
  cat("mean AR estimate over 40 simulate() replicates:", mean(phis), " sd", sd(phis), "\n")
}
run("y")
run("y2")
cat("\n=== brms-45 side: linear gaussian against poisson latent, same machinery\n")
set.seed(167)
dp <- expand.grid(t = 1:5, g = factor(1:30)); dp <- dp[order(dp$g, dp$t), ]
n <- nrow(dp); dp$x <- rnorm(n)
ee <- as.vector(replicate(30, arima.sim(list(ar = 0.5), 5, sd = 0.5)))
dp$count <- rpois(n, exp(0.8 + 0.4 * dp$x + ee))
dp$yg <- 0.8 + 0.4 * dp$x + ee
for (spec in list(list(f = count ~ x + ar(t, g), fam = "poisson", fm = poisson()),
                  list(f = yg ~ x + ar(t, g), fam = "gaussian", fm = gaussian()))) {
  f <- spec$f
  fit <- suppressWarnings(frm(frmtmb::bf(f), family = spec$fm, data = dp))
  fb <- brms::bf(f); fam <- brms::brmsfamily(spec$fam)
  bfit <- suppressWarnings(E$brms_fixed_fit(fb, fam, dp, fit))
  prior <- E$brms_flat_prior(fb, data = dp, family = fam)
  sdat <- E$brms_standata(fb, data = dp, family = fam, prior = prior)
  code <- brms::make_stancode(fb, data = dp, family = fam, prior = prior)
  sf <- suppressMessages(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 0))
  pars <- E$stan_pars_from_fit(fit, sdat, code, NULL)
  lp <- rstan::log_prob(sf, rstan::unconstrain_pars(sf, pars), adjust_transform = FALSE)
  zp <- if (!is.null(pars$zerr)) sum(dnorm(pars$zerr, log = TRUE)) else 0
  ll <- sum(brms::log_lik(bfit)[1L, ])
  cat(sprintf("%s: log_prob %.6f; sum(brms log_lik) + log N(zerr) %.6f; gap %.6g\n", deparse1(f), lp, ll + zp, lp - ll - zp))
}
