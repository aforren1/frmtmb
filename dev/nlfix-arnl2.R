# nl + cov = FALSE: which paths read the ARMA shift? Seed 1582. brms's nl
# path has no predictor_ac() (brms/R/predictor.R, predictor.bprepnl), so its
# posterior_epred is the body alone, while its Stan program's mu carries the
# shift. Measured here at one parameter point, frm() against brms on the nl
# and the linear spelling of one model.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
for (f in list.files("tests/testthat", "helper-brms", full.names = TRUE)) {
  sys.source(f, envir = E)
}
d <- E$eq_dat(seed = 1582)
d$y <- 2 + 0.8 * d$x + rnorm(nrow(d), 0, 0.5)
fl <- frm(bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = d)
fn <- suppressWarnings(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE,
                              autocor = ~ ar(time = t, gr = g, cov = FALSE)),
                           data = d,
                           start = list(beta = c("b0_(Intercept)" = 2,
                                                 "b1_(Intercept)" = 0.5))))
b <- fn$estimates$beta
body <- b[[1]] + b[[2]] * d$x
fv <- fitted(fn)[, "Estimate"]
sh <- fv - body
cat("nl fitted - body, sd:", sd(sh), "\n")
sn <- simulate(fn, nsim = 2000, seed = 3)
sl <- simulate(fl, nsim = 2000, seed = 3)
mn <- rowMeans(as.matrix(sn))
ml <- rowMeans(as.matrix(sl))
se <- fn$estimates$betad
cat("simulate nl: slope of (mean - body) on the shift:",
    coef(lm(I(mn - body) ~ sh))[2], "\n")
cat("simulate linear: slope of (mean - linear body) on the shift:",
    coef(lm(I(ml - (fixef(fl)[1, 1] + fixef(fl)[2, 1] * d$x)) ~ sh))[2], "\n")
pn <- predict(fn, summary = FALSE, ndraws = 2000, propagate_error = FALSE)
cat("predict nl: slope:", coef(lm(I(colMeans(pn) - body) ~ sh))[2], "\n")
pl <- predict(fl, summary = FALSE, ndraws = 2000, propagate_error = FALSE)
cat("predict linear: slope:", coef(lm(I(colMeans(pl) -
  (fixef(fl)[1, 1] + fixef(fl)[2, 1] * d$x)) ~ sh))[2], "\n")
# brms at those values
bn <- brms::bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE,
               autocor = ~ ar(time = t, gr = g, cov = FALSE))
bl <- brms::bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = FALSE))
fxn <- E$brms_fixed_fit(bn, gaussian(), d, fn)
fxl <- E$brms_fixed_fit(bl, gaussian(), d, fl)
en <- brms::posterior_epred(fxn)[1, ]
el <- brms::posterior_epred(fxl)[1, ]
cat("brms nl epred - body, max:", max(abs(en - body)), "\n")
cat("brms linear epred - frm nl fitted, max:", max(abs(el - fv)), "\n")
cat("brms linear epred - frm linear fitted, max:",
    max(abs(el - fitted(fl)[, "Estimate"])), "\n")
fl2 <- E$eq_try
ln <- brms::log_lik(fxn)[1, ]
ll <- brms::log_lik(fxl)[1, ]
cat("brms log_lik sum, nl and linear:", sum(ln), sum(ll), "\n")
cat("frm logLik:", as.numeric(logLik(fn)), "\n")
