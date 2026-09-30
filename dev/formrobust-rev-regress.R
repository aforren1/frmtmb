# Reviewer: base (rellib-r4) against lane on a battery of models.
# Usage: REVLIB="" Rscript formrobust-rev-regress.R base
#        Rscript formrobust-rev-regress.R lane
# Writes dev/formrobust-rev-log/regress-<arm>.rds. Data seed 71; every
# stochastic call carries its own seed.
LIB <- Sys.getenv("REVLIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
arm <- commandArgs(trailingOnly = TRUE)[1]
suppressMessages({library(frmtmb); library(emmeans)})
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(71)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n), u = runif(n),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(rep(1:24, each = 10)),
                t = rep(1:10, 24),
                time = runif(n, 1, 3), w = runif(n, 0.5, 2),
                s = runif(n, 0.2, 0.5), nt = rpois(n, 6) + 2L)
re <- rnorm(24, 0, 0.5)[d$g]
d$y <- 1 + 0.5 * d$x - 0.3 * d$z + re + rnorm(n)
d$y2 <- 0.5 - 0.4 * d$x + rnorm(n)
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x + re) * d$time)
d$yb <- rbinom(n, 1, plogis(0.4 * d$x + re))
d$ybin <- rbinom(n, d$nt, plogis(0.2 * d$x))
d$yo <- cut(d$y, c(-Inf, 0, 1, 2, Inf), labels = FALSE)
d$ycat <- factor(sample(c("p", "q", "r"), n, TRUE))
d$ypos <- exp(0.2 + 0.3 * d$x + rnorm(n, 0, 0.4))
d$ybeta <- plogis(0.3 * d$x + rnorm(n, 0, 0.5))
d$cc <- ifelse(d$y > 2.5, "right", "none")
d$yz <- ifelse(runif(n) < 0.3, 0, d$yc)
e <- as.vector(apply(matrix(rnorm(n), 10, 24), 2, function(z) {
  as.vector(stats::filter(z, 0.5, "recursive"))
}))
d$ya <- 1 + 0.5 * d$x + e
d$sub <- d$u > 0.3
d$xm <- d$x; d$xm[c(5, 17, 33)] <- NA
d$ym <- d$y; d$ym[c(8, 21)] <- NA

M <- list(
  gauss = list(bf(y ~ x + z), gaussian()),
  gauss_re = list(bf(y ~ x + (1 | g)), gaussian()),
  gauss_slope = list(bf(y ~ x + (1 + x | g)), gaussian()),
  gauss_f = list(bf(y ~ x * f), gaussian()),
  sigma_dpar = list(bf(y ~ x, sigma ~ z), gaussian()),
  student = list(bf(y ~ x + (1 | g)), student()),
  pois = list(bf(yc ~ x + (1 | g)), poisson()),
  pois_off = list(bf(yc ~ x + offset(log(time))), poisson()),
  pois_off_re = list(bf(yc ~ x + f + (1 | g) + offset(log(time))),
                     poisson()),
  pois_rate = list(bf(yc | rate(time) ~ x), poisson()),
  negbin = list(bf(yc ~ x + f), negbinomial()),
  bern = list(bf(yb ~ x + (1 | g)), bernoulli()),
  bern_probit = list(bf(yb ~ x + f), bernoulli(link = "probit")),
  binom = list(bf(ybin | trials(nt) ~ x), binomial()),
  weights = list(bf(y | weights(w) ~ x), gaussian()),
  se = list(bf(y | se(s) ~ x), gaussian()),
  se_sigma = list(bf(y | se(s, sigma = TRUE) ~ x), gaussian()),
  cens = list(bf(y | cens(cc) ~ x), gaussian()),
  trunc = list(bf(y | trunc(lb = -5) ~ x), gaussian()),
  cumul = list(bf(yo ~ x + (1 | g)), cumulative()),
  categ = list(bf(ycat ~ x), categorical()),
  lognorm = list(bf(ypos ~ x), lognormal()),
  gamma = list(bf(ypos ~ x + f), Gamma(link = "log")),
  beta = list(bf(ybeta ~ x), Beta()),
  zip = list(bf(yz ~ x, zi ~ z), zero_inflated_poisson()),
  hurdle = list(bf(yz ~ x), hurdle_poisson()),
  nl = list(bf(y ~ a + b * x, a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
            gaussian()),
  smooth = list(bf(y ~ s(x) + z), gaussian()),
  mo = list(bf(y ~ mo(yo) + x), gaussian()),
  ar_cond = list(bf(ya ~ x + ar(t, g)), gaussian()),
  ar_cov = list(bf(ya ~ x + ar(t, g, cov = TRUE)), gaussian()),
  arma_cond = list(bf(ya ~ x + arma(t, g)), gaussian()),
  mv = list(bf(y ~ x) + bf(y2 ~ x + z), gaussian()),
  mv_rescor = list(bf(mvbind(y, y2) ~ x) + set_rescor(TRUE), gaussian()),
  subset = list(bf(y | subset(sub) ~ x) + bf(y2 ~ x), gaussian()),
  mi = list(bf(y ~ mi(xm) + z) + bf(xm | mi() ~ z), gaussian()),
  mi_resp = list(bf(ym | mi() ~ x), gaussian()),
  mixture = list(bf(y ~ x), mixture(gaussian(), gaussian())),
  cs_intercept0 = list(bf(y ~ 0 + Intercept + x), gaussian()))

safe <- function(expr) tryCatch(expr, error = function(e) {
  structure(conditionMessage(e), class = "rev_err")
})
q <- function(expr) suppressWarnings(suppressMessages(expr))
out <- list()
for (nm in names(M)) {
  m <- M[[nm]]
  t0 <- proc.time()[["elapsed"]]
  fit <- safe(q(frm(m[[1]], data = d, family = m[[2]])))
  r <- list()
  if (inherits(fit, "rev_err")) {
    r$fit_error <- as.character(fit)
  } else {
    r$par <- fit$obj$par
    r$fn <- safe(fit$obj$fn(fit$obj$par))
    r$gr <- safe(fit$obj$gr(fit$obj$par))
    r$fn2 <- safe(fit$obj$fn(fit$obj$par * 0.9 + 0.01))
    r$gr2 <- safe(fit$obj$gr(fit$obj$par * 0.9 + 0.01))
    r$estimates <- fit$estimates
    r$logLik <- safe(q(logLik(fit)))
    r$fixef <- safe(q(fixef(fit)))
    r$fitted <- safe(q(fitted(fit)))
    r$fitted_nd <- safe(q(fitted(fit, newdata = d[1:12, ])))
    set.seed(11)
    r$predict <- safe(q(predict(fit, ndraws = 50)))
    set.seed(12)
    r$predict_nd <- safe(q(predict(fit, newdata = d[1:12, ], ndraws = 50)))
    r$residuals <- safe(q(residuals(fit)))
    r$simulate <- safe(q(simulate(fit, nsim = 2, seed = 13)))
    r$ce <- safe(q(lapply(conditional_effects(fit), as.data.frame)))
    r$emm <- safe(q(as.data.frame(summary(emmeans(fit, ~ 1)))))
    r$emm_x <- safe(q(as.data.frame(summary(emmeans(fit, ~ x)))))
    r$emm_ep <- safe(q(as.data.frame(summary(emmeans(fit, ~ 1,
                                                     epred = TRUE)))))
    r$data_names <- names(fit$data)
  }
  r$default_prior <- safe(q(as.data.frame(default_prior(m[[1]], data = d,
                                                        family = m[[2]]))))
  out[[nm]] <- r
  cat(sprintf("%-14s %s %.1fs\n", nm, if (is.null(r$fit_error)) "ok"
              else paste("FIT ERROR:", r$fit_error),
              proc.time()[["elapsed"]] - t0))
}
saveRDS(out, sprintf(
  "C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-log/regress-%s.rds",
  arm))
