# Reviewer, lane ordinal: outputs of 26 models under one arm, for an
# identical() comparison of base (rellib-r4) and lane.
# Usage: Rscript dev/ordinal-rev-bitwise.R <base|lane>
# Data seed 20261002. Output: dev/ordinal-rev-bitwise-<arm>.rds
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
  library(emmeans)
})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
cat("arm", arm, "frmtmb", find.package("frmtmb"), "\n")

set.seed(20261002)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(letters[1:6], n, TRUE)),
                h = factor(sample(c("p", "q"), n, TRUE)))
u <- stats::rlogis(n) + 0.8 * d$x + rnorm(6, 0, 0.4)[as.integer(d$g)]
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yo <- factor(d$y, ordered = TRUE)
d$yh <- ifelse(runif(n) < stats::plogis(-1 + 0.3 * d$x), 0L, d$y)
d$yg <- 2 + 0.5 * d$x + rnorm(n, 0, exp(0.2 * d$z))
d$yc <- rpois(n, exp(0.3 + 0.4 * d$x + rnorm(6, 0, 0.3)[as.integer(d$g)]))
d$yb <- rbinom(n, 1, stats::plogis(0.2 + d$x))
d$ycat <- factor(sample(c("A", "B", "C"), n, TRUE))
d$yz <- ifelse(runif(n) < 0.3, 0L, rpois(n, exp(0.5 + 0.3 * d$x)))
d$yb01 <- stats::plogis(0.3 * d$x + rnorm(n, 0, 0.5))
d$yga <- rgamma(n, 2, 2 / exp(0.3 * d$x))
d$ymix <- ifelse(runif(n) < 0.4, rnorm(n, -2, 1), rnorm(n, 2, 1)) + 0.5 * d$x
M <- t(rmultinom(n, 10, c(0.2, 0.3, 0.5)))
d$ym <- M
d$tot <- 10L

models <- list(
  cum_logit = quote(frm(y ~ x, family = cumulative(), data = d)),
  cum_probit_re = quote(frm(y ~ x + (1 | g), family = cumulative("probit"),
                            data = d)),
  cum_cloglog_gr = quote(frm(y | thres(gr = h) ~ x,
                             family = cumulative("cloglog"), data = d)),
  cum_ordered_factor = quote(frm(yo ~ x + z, family = cumulative(), data = d)),
  sratio_cs = quote(frm(y ~ x + cs(z), family = sratio(), data = d)),
  sratio_cauchit = quote(frm(y ~ x, family = sratio("cauchit"), data = d)),
  cratio_cloglog_re = quote(frm(y ~ x + (1 | g), family = cratio("cloglog"),
                                data = d)),
  cratio_pa = quote(frm(y ~ x, family = cratio("probit_approx"), data = d)),
  acat_cs = quote(frm(y ~ x + cs(z), family = acat(), data = d)),
  acat_gr = quote(frm(y | thres(gr = h) ~ x, family = acat(), data = d)),
  hurdle = quote(frm(bf(yh ~ x, hu ~ x), family = hurdle_cumulative(),
                     data = d)),
  hurdle_disc = quote(frm(bf(yh ~ x, disc ~ 0 + z),
                          family = hurdle_cumulative(), data = d)),
  mv_ord_gauss = quote(frm(bf(y ~ x) + bf(yg ~ x), data = d,
                           family = list(cumulative(), gaussian()))),
  categorical = quote(frm(ycat ~ x, family = categorical(), data = d)),
  multinomial = quote(frm(ym | trials(tot) ~ x, family = multinomial(K = 3),
                          data = d)),
  mixture = quote(frm(ymix ~ x, family = mixture(gaussian(), gaussian()),
                      data = d)),
  poisson_re = quote(frm(yc ~ x + (1 | g), family = poisson(), data = d)),
  gauss_sigma = quote(frm(bf(yg ~ x, sigma ~ z), family = gaussian(),
                          data = d)),
  negbin = quote(frm(yc ~ x, family = negbinomial(), data = d)),
  bernoulli = quote(frm(yb ~ x, family = bernoulli(), data = d)),
  zip = quote(frm(yz ~ x, family = zero_inflated_poisson(), data = d)),
  student = quote(frm(yg ~ x, family = student(), data = d)),
  beta = quote(frm(yb01 ~ x, family = Beta(), data = d)),
  gamma = quote(frm(yga ~ x, family = Gamma("log"), data = d)),
  gauss_nu_fixed = quote(frm(bf(yg ~ x, nu = 4), family = student(),
                             data = d)),
  poisson_mix = quote(frm(yc ~ x, family = mixture(poisson(), poisson()),
                          data = d))
)

grab <- function(expr) tryCatch(suppressWarnings(expr),
                                error = function(e) paste("ERROR:",
                                                          conditionMessage(e)))
out <- list()
for (nm in names(models)) {
  cat("model", nm, "\n")
  fit <- tryCatch(suppressWarnings(eval(models[[nm]])), error = function(e) e)
  if (inherits(fit, "error")) { out[[nm]] <- list(fit_error = conditionMessage(fit)); next }
  r <- list()
  r$par <- fit$opt$par
  r$objective <- fit$opt$objective
  r$logLik <- grab(logLik(fit))
  r$fitted <- grab(fitted(fit))
  r$predict <- grab({set.seed(1); predict(fit, ndraws = 20)})
  r$fixef <- grab(fixef(fit))
  r$fixef_flat <- grab(fixef(fit, flatten = TRUE))
  r$vcov <- grab(vcov(fit))
  r$confint <- grab(confint(fit))
  r$variables <- grab(variables(fit))
  r$summary <- grab(utils::capture.output(print(summary(fit))))
  r$print <- grab(utils::capture.output(print(fit)))
  r$sim <- grab(simulate(fit, nsim = 2, seed = 1))
  r$resid <- grab(residuals(fit))
  r$default_prior <- grab(as.data.frame(default_prior(fit)))
  r$coef <- grab(coef(fit))
  r$hyp <- grab(hypothesis(fit, "x = 0")$hypothesis)
  r$loglik_rows <- grab(log_lik(fit))
  r$emm <- grab(as.data.frame(emmeans(fit, ~ x)))
  r$boot <- grab(frm_bootstrap(fit, nsim = 3, seed = 5)$t)
  r$ce <- grab(conditional_effects(fit, "x")[[1]])
  if (nm %in% c("cum_logit", "hurdle", "acat_cs", "sratio_cs",
                "cum_probit_re", "gauss_sigma", "poisson_re")) {
    ds <- grab(suppressMessages(frm_sample(fit, chains = 1, iter = 200,
                                           refresh = 0, seed = 7)))
    r$draws_vars <- grab(variables(ds))
    r$draws <- grab(as.matrix(ds))
    r$draws_fixef <- grab(fixef(ds))
    r$draws_summary <- grab(utils::capture.output(print(summary(ds))))
    r$draws_epred <- grab(posterior_epred(ds, ndraws = 3))
  }
  out[[nm]] <- r
}
saveRDS(out, file.path(wt, paste0("dev/ordinal-rev2-bitwise-", arm, ".rds")))
cat("done\n")
