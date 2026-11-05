# Item 3: on the family with the largest shift, is brms's posterior for the
# other parameters the complete-case one, and where do frm()'s two routes sit?
# exgaussian, sigma 0.5, beta 2, 200 rows, 60 missing. Seed 10600.
source("dev/review-refuse-lib.R")
set.seed(10600)
n <- 200; x <- rnorm(n)
d <- data.frame(x = x, y = 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2))
rows <- sample(n, 60); d$y[rows] <- NA; dc <- d[-rows, ]
fb <- brms::exgaussian()
one <- function(form, dat) {
  pri <- brms::default_prior(form, data = dat, family = fb)
  code <- brms::make_stancode(form, data = dat, family = fb, prior = pri)
  sdat <- E$brms_standata(form, data = dat, family = fb, prior = pri)
  suppressWarnings(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 4, iter = 2000, seed = 11, refresh = 0, cores = 1))
}
sj <- one(brms::bf(y | mi() ~ x), d); sc <- one(brms::bf(y ~ x), dc)
ej <- rstan::extract(sj); ec <- rstan::extract(sc)
for (nm in c("b", "Intercept", "sigma", "beta")) {
  pr("brms %-10s joint %.4f (sd %.4f); observed rows %.4f (sd %.4f); z %.2f", nm, mean(ej[[nm]]), sd(ej[[nm]]), mean(ec[[nm]]), sd(ec[[nm]]),
     (mean(ej[[nm]]) - mean(ec[[nm]])) / sd(ec[[nm]]))
}
a <- frm(bf(y | mi() ~ x), family = exgaussian(), data = d)
b <- frm(bf(y ~ x), family = exgaussian(), data = dc)
pr("frm: joint b %.4f sigma %.4f beta %.4f; observed rows b %.4f sigma %.4f beta %.4f", a$estimates$beta[[2]], exp(a$estimates$betad[[1]]), exp(a$estimates$betad[[2]]),
   b$estimates$beta[[2]], exp(b$estimates$betad[[1]]), exp(b$estimates$betad[[2]]))
pr("frm: logLik joint %.3f, observed rows %.3f", as.numeric(logLik(a)), as.numeric(logLik(b)))
pr("brms value at the posterior-mean Ymi vs fitted: Ymi mean %.3f, sd %.3f at row 1; frm mode %.3f, fitted mean %.3f", mean(ej$Ymi[, 1]), sd(ej$Ymi[, 1]),
   frmtmb:::mi_values(a, "y")[sort(rows)[1]], fitted(a)[sort(rows)[1], 1])

# mi() + cens() with no missing value: a guard on a correct model?
set.seed(10004); n <- 180; x <- rnorm(n)
d0 <- data.frame(x = x, y = rnorm(n, 0.5 * x), ce = rbinom(n, 1, 0.1))
form <- brms::bf(y | mi() + cens(ce) ~ x)
pri <- E$brms_flat_prior(form, data = d0, family = stats::gaussian())
code <- brms::make_stancode(form, data = d0, prior = pri)
sdat <- E$brms_standata(form, data = d0, prior = pri)
s <- tr(suppressWarnings(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 1, iter = 200, refresh = 0, seed = 1)))
pr("brms mi() + cens() with no NA: draws kept %s", if (is_err(s)) msg(s, 100) else sum(s@sim$n_save))
form2 <- brms::bf(y | mi() + trunc(lb = -50) ~ x)
a <- tr(frm(bf(y | mi() + cens(ce) ~ x), data = d0)); pr("frm mi() + cens() no NA: %s", msg(a, 120))
a <- tr(frm(bf(y | mi() + trunc(lb = -50) ~ x), data = d0)); pr("frm mi() + trunc() no NA: %s", msg(a, 120))
a <- tr(frm(bf(y | mi() ~ x), family = asym_laplace(), data = d0)); pr("frm mi() asym_laplace no NA: %s", msg(a, 120))
